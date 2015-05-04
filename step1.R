works_with_R("3.2.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@c68c566a606aea75d646e81087d168e5e2f0531a")

load("train.sets.RData")
load("chunk.problems.RData")

### Compute the error of model fit with respect to chunk.name.
chunkError <- function(fit, chunk.name){
  data.by.problem <- data.by.chunk[[chunk.name]][[res.str]]
  peaks.by.regularization <- list()
  for(problem.name in names(data.by.problem)){
    problem <- data.by.problem[[problem.name]]
    log.lambda.vec <- fit$predict(problem$features)
    stopifnot(is.finite(log.lambda.vec))
    for(log.lambda.i in seq_along(log.lambda.vec)){
      regularization <- fit$regularization.vec[[log.lambda.i]]
      reg.str <- paste(regularization)
      log.lambda <- log.lambda.vec[[log.lambda.i]]
      selected.row <- 
        subset(problem$modelSelection,
               min.log.lambda < log.lambda &
                 log.lambda < max.log.lambda)
      peaks.by.regularization[[reg.str]][[problem.name]] <-
        problem$peaks[[paste(selected.row$peaks)]]
    }#log.lambda.i
  }#problem.name
  chunk.err.list <- list()
  peaks.out.list <- list()
  for(regularization.i in seq_along(peaks.by.regularization)){
    peaks.by.problem <- peaks.by.regularization[[regularization.i]]
    regularization <- fit$regularization.vec[[regularization.i]]
    reg.peaks <- do.call(rbind, peaks.by.problem)
    cluster.peaks <- if(nrow(reg.peaks) == 0){
      Peaks()
    }else{
      cl.peaks <- clusterPeaks(reg.peaks)
      peaks.by.cluster <- split(cl.peaks, cl.peaks$cluster)
      cluster.peaks.list <- list()
      for(cluster.name in names(peaks.by.cluster)){
        one.cluster <- peaks.by.cluster[[cluster.name]]
        cluster.peaks.list[[cluster.name]] <- 
          with(one.cluster, {
            data.frame(chromStart=min(chromStart),
                       chromEnd=max(chromEnd),
                       sample.id=unique(sample.id))
          })
      }#cluster.name
      do.call(rbind, cluster.peaks.list)
    }
    peaks.out.list[[paste(regularization)]] <- cluster.peaks
    regions.file <-
      sprintf("../chip-seq-paper/chunks/%s/regions.RData",
              chunk.name)
    regions <- regions.file.list[[regions.file]]
    stopifnot(is.data.frame(regions))
    stopifnot(nrow(regions) > 0)
    regions.by.sample <- split(regions, regions$sample.id)
    peaks.by.sample <-
      split(cluster.peaks, cluster.peaks$sample.id, drop=FALSE)
    for(sample.id in names(regions.by.sample)){
      sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
        peaks.by.sample[[sample.id]]
      }else{
        Peaks()
      }
      sample.regions <- regions.by.sample[[sample.id]]
      error.regions <- PeakErrorChrom(sample.peaks, sample.regions)
      chunk.err.list[[paste(sample.id, regularization)]] <- 
        data.frame(regularization, sample.id,
                   error.regions)
    }#sample.id  
  }#regularization
  list(errors=do.call(rbind, chunk.err.list),
       peak.list=peaks.out.list)
}

regions.file.list <- list()
regions.file.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/regions.RData")
for(regions.file in regions.file.vec){
  load(regions.file)
  regions.file.list[[regions.file]] <- regions
}

all.split.list <- list()

set.name <- "srf_NTNU_Gm12878"
split.i <- 1
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  data.by.chunk <- chunk.problems[[set.name]]
  for(split.i in seq_along(splits)){
    split.name <- paste(set.name, "split", split.i)
    under.name <- gsub(" ", "_", split.name)
    print(split.name)
    train.validation <- splits[[split.i]]
    set.seed(1)
    train.i <- sample(length(train.validation), length(train.validation)/2)
    sets <-
      list(train=train.validation[train.i],
           validation=train.validation[-train.i])
    model.error.list <- list() 
    fit.list <- list()
    pred.peak.list <- list()
    for(res.str in names(data.by.chunk[[1]])){
      RData.file <- sprintf("step1/%s_%s.RData", under.name, res.str)
      if(file.exists(RData.file)){
        load(RData.file)
      }else{
        train.data.list <- list()
        for(train.chunk in sets$train){
          data.by.problem <- data.by.chunk[[train.chunk]][[res.str]]
          for(problem.name in names(data.by.problem)){
            problem <- data.by.problem[[problem.name]]
            train.data.list[[paste(train.chunk, problem.name)]] <- problem
          }
        }
        fit <- IntervalRegressionProblems(train.data.list)
        tv.err.list <- list()
        tv.peak.list <- list()
        for(tv in names(sets)){
          set.chunks <- sets[[tv]]
          for(chunk.name in set.chunks){
            pred.list <- chunkError(fit, chunk.name)
            tv.peak.list[[chunk.name]] <- pred.list$peak.list
            tv.err.list[[chunk.name]] <-
              data.frame(tv, chunk.name, pred.list$errors)
          }#chunk.name
        }#tv
        tv.err.regions <- do.call(rbind, tv.err.list)
        tv.reg.err <- tv.err.regions %>%
          group_by(tv, regularization) %>%
          summarise(fp=sum(fp),
                    fn=sum(fn),
                    regions=n()) %>%
          mutate(errors=fp+fn)
        v.min.err <- tv.reg.err %>%
          filter(tv=="validation") %>%
          filter(errors == min(errors))
        resPlot <- 
          ggplot()+
            ggtitle(paste(res.str, "bases/problem"))+
            geom_point(aes(-log10(regularization), errors,
                           color=tv),
                       size=3,
                       pch=1,
                       data=v.min.err)+
            geom_line(aes(-log10(regularization), errors,
                          color=tv, group=tv),
                      data=tv.reg.err)+
            xlab("model complexity -log10(regularization)")
        print(resPlot)
        
        RData.dir <- dirname(RData.file)
        dir.create(RData.dir, showWarnings=FALSE, recursive=TRUE)
        save(fit, tv.reg.err, tv.peak.list, file=RData.file)
      }#if(file.exists(RData.file)) else
      model.error.list[[res.str]] <-
        data.frame(bases.per.problem=as.numeric(res.str),
                   tv.reg.err)
      pred.peak.list[[res.str]] <- tv.peak.list
      fit.list[[res.str]] <- fit
    }#res.str
    model.error <- do.call(rbind, model.error.list)
    glob.min.err <- model.error %>%
      filter(tv=="validation") %>%
      filter(errors == min(errors)) %>%
      tail(1)
    hline.min.err <- glob.min.err %>%
      select(-bases.per.problem)
    resPlot <- 
      ggplot()+
        geom_hline(aes(yintercept=errors, color=tv),
                   data=hline.min.err)+
        geom_point(aes(-log10(regularization), errors,
                       color=tv),
                   size=3,
                   pch=1,
                   data=glob.min.err)+
        geom_line(aes(-log10(regularization), errors,
                      color=tv, group=tv),
                  data=model.error)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(. ~ bases.per.problem)+
        xlab("model complexity -log10(regularization)")
    print(resPlot)
    res.str <- paste(glob.min.err$bases.per.problem)
    reg.str <- paste(glob.min.err$regularization)
    weight.vec <- fit.list[[res.str]]$weight.mat[, reg.str]
    nonzero.vec <- weight.vec[weight.vec != 0]
    print(nonzero.vec)

    out.peak.list <- list()
    res.peak.list <- pred.peak.list[[res.str]]
    for(chunk.name in names(res.peak.list)){
      out.peak.list[[chunk.name]] <- res.peak.list[[chunk.name]][[reg.str]]
    }
    split.data <-
      list(glob.min.err=glob.min.err,
           res.str=res.str,
           reg.str=reg.str,
           fit=fit.list[[res.str]],
           peak.list=out.peak.list)
    all.split.list[[split.name]] <- split.data
  }#split.i
}#set.name

hyper.list <- list()
seg.list <- list()
for(split.name in names(all.split.list)){
  hyper.row <- all.split.list[[split.name]]$glob.min.err
  set.name <- sub(" .*", "", split.name)
  data.by.chunk <- chunk.problems[[set.name]]
  bases.per.problem.vec <- as.numeric(names(data.by.chunk[[1]]))
  seg.list[[split.name]] <-
    data.frame(set.name, split.name,
               min.bases=min(bases.per.problem.vec),
               max.bases=max(bases.per.problem.vec))
  hyper.list[[split.name]] <- data.frame(set.name, split.name, hyper.row)
}
hyper.tall <- do.call(rbind, hyper.list)
seg.tall <- do.call(rbind, seg.list)
dcast(hyper.tall, set.name ~ bases.per.problem, fun.aggregate=length)
ggplot()+
  geom_segment(aes(min.bases, set.name,
                   xend=max.bases, yend=set.name),
               data=seg.tall,
               size=5)+
  scale_x_log10()+
  geom_point(aes(bases.per.problem, set.name),
             color="red",
             data=hyper.tall)

step1 <- all.split.list
save(step1, file="step1.RData")
