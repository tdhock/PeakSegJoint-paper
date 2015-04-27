load("train.sets.RData")
load("chunk.problems.RData")

IntervalRegressionProblems <- function
### Compute a sequence of interval regression models for increasingly
### more L1 regularization, until we get to a regularization parameter
### that gives an optimal weight vector of zero.
(problem.list,
### List of problems with features (numeric matrix) and target
### (numeric vector of length 2).
 initial.regularization=0.001,
 factor.regularization=1.5,
 ...
 ){
  stopifnot(is.list(problem.list))
  stopifnot(is.numeric(initial.regularization))
  stopifnot(length(initial.regularization) == 1)
  stopifnot(initial.regularization > 0)
  stopifnot(is.numeric(factor.regularization))
  stopifnot(length(factor.regularization) == 1)
  stopifnot(factor.regularization > 0)
  
  n.input.features <- ncol(problem.list[[1]]$features)

  featureSum.list <- list()
  targets.list <- list()
  for(problem.i in seq_along(problem.list)){
    stopifnot(is.matrix(problem$features))
    stopifnot(is.numeric(problem$features))
    stopifnot(ncol(problem$features) == n.input.features)
    stopifnot(is.numeric(problem$target))
    stopifnot(length(problem$target) == 2)
    problem <- problem.list[[problem.i]]
    featureSum.list[[problem.i]] <- colSums(problem$features)
    targets.list[[problem.i]] <- problem$target
  }
  all.targets <- do.call(rbind, targets.list)
  is.trivial.target <- apply(!is.finite(all.targets), 1, all)
  nontrivial.i.vec <- which(!is.trivial.target)

  all.featureSum <- do.call(rbind, featureSum.list)
  is.finite.feature <- apply(is.finite(all.featureSum), 2, all)

  features.list <- list()
  for(problem.i in nontrivial.i.vec){
    problem <- problem.list[[problem.i]]
    features.list[[paste(problem.i)]] <-
      problem$features[, is.finite.feature, drop=FALSE]
  }
  finite.features <- do.call(rbind, features.list)
  all.mean.vec <- colMeans(finite.features)
  all.sd.vec <- apply(finite.features, 2, sd)
  is.invariant <- all.sd.vec == 0
  train.feature.i <- which(!is.invariant)
  train.feature.names <- colnames(finite.features)[train.feature.i]
  mean.vec <- all.mean.vec[train.feature.names]
  sd.vec <- all.sd.vec[train.feature.names]

  norm.featureSum.list <- list()
  for(problem.i in nontrivial.i.vec){
    problem <- problem.list[[problem.i]]
    raw.features.mat <- problem$features[, train.feature.names, drop=FALSE]
    train.mean.mat <-
      matrix(mean.vec, nrow(raw.features.mat), ncol(raw.features.mat),
             byrow=TRUE)
    train.sd.mat <-
      matrix(sd.vec, nrow(raw.features.mat), ncol(raw.features.mat),
             byrow=TRUE)
    norm.mat <- (raw.features.mat-train.mean.mat)/train.sd.mat
    intercept.mat <- cbind("(Intercept)"=1, norm.mat)
    norm.featureSum.list[[paste(problem.i)]] <- colSums(intercept.mat)
  }
  norm.featureSum.mat <- do.call(rbind, norm.featureSum.list)
  targets.mat <- do.call(rbind, targets.list[nontrivial.i.vec])
  rownames(norm.featureSum.mat) <- rownames(targets.mat) <-
    names(problem.list)[nontrivial.i.vec]

  ## Do we need to take the feature sd again and divide by that for
  ## optimization purposes?
  apply(norm.featureSum.mat, 2, sd)

  regularization <- initial.regularization
  n.features <- ncol(norm.featureSum.mat)
  weight.vec <- rep(0, n.features)
  n.nonzero <- n.features
  
  weight.vec.list <- list()
  while(n.nonzero > 1){
    weight.vec <-
      IntervalRegressionMatrix(norm.featureSum.mat, targets.mat,
                               weight.vec,
                               regularization)
    weight.vec.list[[paste(regularization)]] <- weight.vec
    regularization <- regularization * factor.regularization
  }
}

IntervalRegressionMatrix <- function
### Solve the squared hinge loss interval regression problem for one
### regularization parameter.
(features,
### Numeric feature matrix (problems x features).
 targets,
### Numeric target matrix (problems x 2).
 initial.weight.vec,
### initial guess for weight vector (features).
 regularization
### Degree of L1-regularization. 
 ){
  stopifnot(is.matrix(features))
  stopifnot(is.numeric(features))
  n.features <- ncol(features)
  n.problems <- nrow(features)

  stopifnot(is.matrix(targets))
  stopifnot(nrow(targets) == n.problems)
  stopifnot(ncol(targets) == 2)

  stopifnot(is.numeric(initial.weight.vec))
  stopifnot(length(initial.weight.vec) == n.features)

  weight.vec <- initial.weight.vec

  stop("optimization while loop")

  weight.vec
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
      fit <- IntervalRegressionProblems(train.data.list)
    }
  }
}
