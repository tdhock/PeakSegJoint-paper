works_with_R("3.2.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@9d1a4fc6751334513f2abb07689a187e04e4db76")

load("selected.by.set.RData")
load("train.sets.RData")

regions.file.list <- list()
regions.file.vec <- Sys.glob("PeakSegJoint-chunks/*/*/regions.RData")
for(regions.file in regions.file.vec){
  load(regions.file)
  regions.file.list[[regions.file]] <- regions
}

RData.vec <- Sys.glob("PeakSegJoint-chunks/*/*/problems.RData")
chunk.list <- list()
for(RData in RData.vec){
  objs <- load(RData)
  chunk.list[[RData]] <- step2.by.res
}

step1.error.list <- list()
for(set.name in names(chunk.problems)){
  regions.RData.vec <-
    Sys.glob(file.path("PeakSegJoint-chunks", set.name, "*", "regions.RData"))
  chunk.ids <- basename(dirname(regions.RData.vec))
  set.chunks <- paste0(set.name, "/", chunk.ids)
  set.error.df <- selected.by.set[[set.name]]
  for(split.i in 1:6){
    bases.per.problem <- set.error.df[split.i, "bases.per.problem"]
    res.str <- paste(bases.per.problem)
    train.validation <- train.sets[[set.name]][[split.i]]
    n.folds <- if(length(train.validation)==2) 2 else 3
    set.seed(1)
    fold.id <- sample(rep(1:n.folds, l=length(train.validation)))
    for(validation.fold in 1:n.folds){
      is.validation <- fold.id == validation.fold
      sets <- list(validation=train.validation[is.validation],
                   train=train.validation[!is.validation])
      train.list <- list()
      for(chunk.name in sets$train){
        problems.RData <-
          paste0("PeakSegJoint-chunks/", chunk.name, "/problems.RData")
        step2.by.res <- chunk.list[[problems.RData]]
        step2.by.problem <- step2.by.res[[res.str]]
        for(problem.name in names(step2.by.problem)){
          problem <- step2.by.problem[[problem.name]]
          train.list[[problem.name]] <-
            list(features=problem$features,
                 target=problem$error$target)
        }
      }
      fit <- IntervalRegressionProblems(train.list)
      stop("compute train/validation error curves")
    }
    is.train <- set.chunks %in% train.chunks
    test.chunks <- set.chunks[!is.train]
    split.name <- paste(set.name, "split", split.i)
    print(split.name)
    model.info <- step1[[split.name]]
    chunk.err.list <- list()
    for(chunk.name in test.chunks){
      data.by.problem <- data.by.chunk[[chunk.name]][[model.info$res.str]]
      peaks.by.problem <- list()
      for(problem.name in names(data.by.problem)){
        problem <- data.by.problem[[problem.name]]
        log.lambda.vec <- model.info$fit$predict(problem$features)
        log.lambda <- log.lambda.vec[[model.info$reg.str]]
        selected <- 
          subset(problem$modelSelection,
                 min.log.lambda < log.lambda & log.lambda < max.log.lambda)
        stopifnot(nrow(selected) == 1)
        peaks.by.problem[[problem.name]] <-
          problem$peaks[[paste(selected$peaks)]]
      }#problem.name
      chunk.peaks <- do.call(rbind, peaks.by.problem)
      pred.peaks <- if(nrow(chunk.peaks) == 0){
        Peaks()
      }else{
        clustered.peaks <- clusterPeaks(chunk.peaks)
        peaks.by.cluster <- split(clustered.peaks, clustered.peaks$cluster)
        pred.by.cluster <- list()
        for(cluster.name in names(peaks.by.cluster)){
          one.cluster <- peaks.by.cluster[[cluster.name]]
          pred.by.cluster[[cluster.name]] <- with(one.cluster, {
            data.frame(chromStart=min(chromStart),
                       chromEnd=max(chromEnd),
                       sample.id=unique(sample.id))
          })
        }
        do.call(rbind, pred.by.cluster)
      }
      regions.file <-
        sprintf("../chip-seq-paper/chunks/%s/regions.RData",
                chunk.name)
      regions <- regions.file.list[[regions.file]]
      stopifnot(is.data.frame(regions))
      stopifnot(nrow(regions) > 0)
      regions.by.sample <- split(regions, regions$sample.id)
      peaks.by.sample <-
        split(pred.peaks, pred.peaks$sample.id, drop=TRUE)
      for(sample.id in names(regions.by.sample)){
        sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
          peaks.by.sample[[sample.id]]
        }else{
          Peaks()
        }
        sample.regions <- regions.by.sample[[sample.id]]
        error.regions <- PeakErrorChrom(sample.peaks, sample.regions)
        chunk.err.list[[paste(chunk.name, sample.id)]] <- 
          data.frame(chunk.name, sample.id, error.regions)
      }#sample.id  
    }#chunk.name
    chunk.err <- do.call(rbind, chunk.err.list)
    step1.error.list[[split.name]] <-
      data.frame(set.name, split.i, chunk.err)
  }#split.i
}#set.name

step1.error <- do.call(rbind, step1.error.list)

save(step1.error, file="step1.error.RData")
