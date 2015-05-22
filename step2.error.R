works_with_R("3.2.0",
             data.table="1.9.4",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@b127ca0b27885f8058665284fd276c958e39f386")

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
  e <- new.env()
  objs <- load(RData, e)
  chunk.list[[RData]] <- e
}

problemData <- function(chunk.name.vec, FUN){
  data.by.problem <- list()
  for(chunk.name in chunk.name.vec){
    problems.RData <-
      paste0("PeakSegJoint-chunks/", chunk.name, "/problems.RData")
    chunk.env <- chunk.list[[problems.RData]]
    problems.dt <- chunk.env$step2.data.list[[res.str]]$problems
    for(problem.name in paste(problems.dt$problem.name)){
      problem <- chunk.env$step2.model.list[[problem.name]]
      if(problem.name %in% names(chunk.env$step2.error.list)){
        plist <- chunk.env$step2.error.list[[problem.name]]$problem
        problem[names(plist)] <- plist
      }
      data.by.problem[[problem.name]] <- FUN(problem)
    }
  }
  data.by.problem
}

labeledProblems <- function(chunk.name.vec, FUN){
  problemData(chunk.name.vec, function(problem){
    if("error.totals" %in% names(problem)){
      FUN(problem)
    }
  })
}

getTrain <- function(problem){
  if(is.numeric(problem$target)){
    list(features=problem$features,
         target=problem$target)
  }
}

estimate.regularization <- function(train.validation){
  n.folds <- if(length(train.validation)==2) 2 else 3
  set.seed(1)
  fold.id <- sample(rep(1:n.folds, l=length(train.validation)))
  picked.by.fold <- list()
  for(validation.fold in 1:n.folds){
    is.validation <- fold.id == validation.fold
    sets <- list(validation=train.validation[is.validation],
                 train=train.validation[!is.validation])
    train.list <- labeledProblems(sets$train, getTrain)
    fit <-
      IntervalRegressionProblems(train.list,
                                 initial.regularization=0.005,
                                 factor.regularization=1.1,
                                 verbose=0)
    error.by.tv <- list()
    for(tv in names(sets)){
      error.vec.list <- labeledProblems(sets[[tv]], function(problem){
        log.lambda.vec <- fit$predict(problem$features)
        error.vec <- rep(NA, length(log.lambda.vec))
        for(log.lambda.i in seq_along(log.lambda.vec)){
          log.lambda <- log.lambda.vec[[log.lambda.i]]
          selected <- 
            subset(problem$modelSelection,
                   min.log.lambda < log.lambda &
                     log.lambda < max.log.lambda)
          stopifnot(nrow(selected) == 1)
          error.vec[[log.lambda.i]] <- selected$error
        }#log.lambda
        error.vec
      })
      error.mat <- do.call(rbind, error.vec.list)
      error.by.tv[[tv]] <- 
        data.frame(tv,
                   errors=colSums(error.mat),
                   regularization=fit$regularization)
    }#tv
    tv.error <- do.call(rbind, error.by.tv)
    picked.i <- pick.best.index(error.by.tv$validation$errors)
    picked.error <- error.by.tv$validation[picked.i, ]
    tvPlot <- 
      ggplot()+
        ggtitle(paste(split.name, "validation fold", validation.fold))+
        geom_point(aes(-log10(regularization), errors, color=tv),
                   pch=1,
                   data=picked.error)+
        geom_line(aes(-log10(regularization), errors, color=tv),
                  data=tv.error)
    print(tvPlot)
    picked.by.fold[[validation.fold]] <- picked.error
  }#validation.fold
  picked <- do.call(rbind, picked.by.fold)
  mean(picked$regularization)
}  

step2.error.list <- list()
for(set.name in names(selected.by.set)){
  regions.RData.vec <-
    Sys.glob(file.path("PeakSegJoint-chunks", set.name, "*", "regions.RData"))
  chunk.ids <- basename(dirname(regions.RData.vec))
  set.chunks <- paste0(set.name, "/", chunk.ids)
  set.error.df <- selected.by.set[[set.name]]
  for(split.i in 1:6){
    split.name <- paste(set.name, "split", split.i)
    print(split.name)
    bases.per.problem <- set.error.df[split.i, "bases.per.problem"]
    res.str <- paste(bases.per.problem)
    train.chunk.vec <- train.sets[[set.name]][[split.i]]
    reg.est <- estimate.regularization(train.chunk.vec)
    train.list <- labeledProblems(train.chunk.vec, getTrain)
    fit <- IntervalRegressionProblems(train.list,
                                      initial.regularization=reg.est)
    
    is.train <- set.chunks %in% train.chunk.vec
    test.chunks <- set.chunks[!is.train]

    error.by.problem <- labeledProblems(test.chunks, function(problem){
      log.lambda <- fit$predict(problem$features)[1]
      selected <- 
        subset(problem$modelSelection,
               min.log.lambda < log.lambda &
                 log.lambda < max.log.lambda)
      stopifnot(nrow(selected) == 1)
      selected
    })

    split.error <- do.call(rbind, error.by.problem)
    
    step2.error.list[[split.name]] <-
      data.frame(set.name, split.i, split.error)
  }#split.i
}#set.name

step2.error <- do.call(rbind, step2.error.list)

save(step2.error, file="step2.error.RData")
