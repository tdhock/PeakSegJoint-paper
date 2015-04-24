load("train.sets.RData")
load("chunk.problems.RData")

matrixIntervalRegression <- function(problem.list){
  n.features <- ncol(problem.list$features)
  n.problems <- length(problem.list)

  features.list <- list()
  for(problem.i in seq_along(problem.list)){
    stopifnot(is.matrix(problem$features))
    stopifnot(is.numeric(problem$features))
    stopifnot(ncol(problem$features) == n.features)
    stopifnot(is.numeric(problem$target))
    stopifnot(length(problem$target) == 2)
    features.list[[problem.i]] <- colSums(problem$features)
  }
  all.features <- do.call(rbind, features.list)
  all.targets <- do.call(rbind, lapply(problem.list, "[[", "target"))
  bad.target <- apply(!is.finite(all.targets), 1, all)
  bad.feature <- stop("filter bad features")
  exclude <- bad.target | bad.feature
  stop("add constant feature = intercept")
  pred.penalty <- features %*% current.iterate
}

set.seed(1)
set.name <- "srf_NTNU_Gm12878"
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  data.by.chunk <- chunk.problems[[set.name]]
  for(split.i in seq_along(splits)){
    train.validation <- splits[[split.i]]
    train.i <- sample(length(train.validation), length(train.validation)/2)
    sets <-
      list(train=train.validation[train.i],
           validation=train.validation[-train.i])
    model.error.list <- list() # [[res.str, model.complexity]]
    for(res.str in names(data.by.chunk[[1]])){
      train.data.list <- list()
      for(train.chunk in sets$train){
        data.by.problem <- data.by.chunk[[train.chunk]][[res.str]]
        for(problem.name in names(data.by.problem)){
          problem <- data.by.problem[[problem.name]]
          train.data.list[[paste(train.chunk, problem.name)]] <- problem
        }
      }
      fit <- matrixIntervalRegression(train.data.list)
    }
  }
}
