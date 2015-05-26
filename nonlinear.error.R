works_with_R("3.2.0",
             data.table="1.9.4",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@b127ca0b27885f8058665284fd276c958e39f386")

load("selected.by.set.RData")
load("train.sets.RData")

IntervalRegressionMatrix <- function
### Solve the squared hinge loss interval regression problem for one
### regularization parameter: w* = argmin_w L(w) + regularization *
### ||w||_1 where L(w) is the average squared hinge loss with respect
### to the targets, and ||w||_1 is the L1-norm of the weight vector
### (excluding the first element, which is the un-regularized
### intercept or bias term).
(features,
### Scaled numeric feature matrix (problems x features). The first
### column/feature should be all ones and will not be regularized.
 targets,
### Numeric target matrix (problems x 2).
 initial.weight.vec,
### initial guess for weight vector (features).
 regularization,
### Degree of L1-regularization.
 threshold=1e-2,
### When the stopping criterion gets below this threshold, the
### algorithm stops and declares the solution as optimal.
 max.iterations=1e4,
### Error if the algorithm has not found an optimal solution after
### this many iterations.
 Lipschitz=NULL,
### A numeric scalar or NULL, which means to compute Lipschitz as the
### mean of the squared L2-norms of the rows of the feature matrix.
 verbose=2
### Print messages if >= 2.
 ){
  n.features <- ncol(features)
  n.problems <- nrow(features)

  stopifnot(nrow(targets) == n.problems)
  stopifnot(ncol(targets) == 2)

  if(is.null(Lipschitz)){
    Lipschitz <- mean(rowSums(features * features))
  }
  stopifnot(is.numeric(Lipschitz))
  stopifnot(length(Lipschitz) == 1)

  stopifnot(is.numeric(max.iterations))
  stopifnot(length(max.iterations) == 1)

  stopifnot(is.numeric(threshold))
  stopifnot(length(threshold) == 1)

  stopifnot(length(initial.weight.vec) == n.features)

  ## Return 0 for a negative number and the same value otherwise.
  positive.part <- function(x){
    ifelse(x<0, 0, x)
  }
  squared.hinge <- function(x){
    ifelse(x<1,(x-1)^2,0)
  }
  squared.hinge.deriv <- function(x){
    ifelse(x<1,2*(x-1),0)
  }  
  calc.loss <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge(linear.predictor-targets[,1])
    right.term <- squared.hinge(targets[,2]-linear.predictor)
    mean(left.term+right.term)
  }
  calc.grad <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge.deriv(linear.predictor-targets[,1])
    right.term <- squared.hinge.deriv(targets[,2]-linear.predictor)
    full.grad <- features * (left.term-right.term)
    colSums(full.grad)/nrow(full.grad)
  }    
  calc.penalty <- function(x){
    regularization * sum(abs(x[-1]))
  }
  calc.cost <- function(x){
    calc.loss(x) + calc.penalty(x)
  }
  soft.threshold <- function(x,thresh){
    ifelse(abs(x) < thresh, 0, x-thresh*sign(x))
  }
  ## do not threshold the intercept.
  prox <- function(x,thresh){
    x[-1] <- soft.threshold(x[-1],thresh)
    x
  }
  ## p_L from the fista paper.
  pL <- function(x,L){
    grad <- calc.grad(x)
    prox(x - grad/L, regularization/L)
  }
  dist2subdiff.opt <- function(w,g){
    ifelse(w==0,positive.part(abs(g)-regularization),
           ifelse(w<0,abs(-regularization+g),abs(regularization+g)))
  }

  iterate.count <- 1
  stopping.crit <- threshold
  last.iterate <- this.iterate <- y <- initial.weight.vec
  this.t <- 1
  while(stopping.crit >= threshold){
    ## here we implement the FISTA method with constant step size, as
    ## described by in the Beck and Tebolle paper.
    last.iterate <- this.iterate
    this.iterate <- pL(y, Lipschitz)
    last.t <- this.t
    this.t <- (1+sqrt(1+4*last.t^2))/2
    y <- this.iterate + (last.t - 1)/this.t*(this.iterate-last.iterate)
    ## here we calculate the subgradient optimality condition, which
    ## requires 1 more gradient evaluation per iteration.
    after.grad <- calc.grad(this.iterate)
    w.dist <- dist2subdiff.opt(this.iterate[-1],after.grad[-1])
    zero.at.optimum <- c(abs(after.grad[1]),w.dist)
    stopping.crit <- max(zero.at.optimum)

    if(verbose >= 2){
      cost <- calc.cost(this.iterate)
      cat(sprintf("%10d cost %10f crit %10.7f\n",
                  iterate.count,
                  cost,
                  stopping.crit))
    }
    iterate.count <- iterate.count + 1
    if(iterate.count > max.iterations){
      Lipschitz <- Lipschitz * 1.5
      iterate.count <- 1
      cat(max.iterations," iterations, increasing Lipschitz.\n")
    }
  }
  this.iterate
### Numeric vector of scaled weights w of the affine function f_w(X) =
### X %*% w for a scaled feature matrix X with the first row entirely
### ones.
}

NonlinearProblems <- structure(function
### Compute a sequence of interval regression models for increasingly
### more L1 regularization, until we get to a regularization parameter
### that gives an optimal weight vector of zero. The problem is w* =
### argmin_w L(w) + regularization * ||w||_1 where L(w) is the mean
### squared hinge loss and ||w||_1 is the L1-norm of the non-intercept
### coefficients. We first scale the input features and then
### repeatedly call IntervalRegressionMatrix, using warm restarts.
(problem.list,
### List of problems with features (numeric matrix) and target
### (numeric vector of length 2).
 initial.regularization=0.001,
### Initial regularization parameter.
 threshold=initial.regularization,
 factor.regularization=1.5,
### Increase regularization by this factor after finding an optimal
### solution.
 verbose=1,
### Print messages if >= 1.
 featureVec=function(mat){
   x <- 
     c(sum=colSums(mat, na.rm=TRUE),
       quartile.0=apply(mat, 2, quantile, na.rm=TRUE, 0),
       quartile.25=apply(mat, 2, quantile, na.rm=TRUE, 0.25),
       quartile.50=apply(mat, 2, quantile, na.rm=TRUE, 0.50),
       quartile.75=apply(mat, 2, quantile, na.rm=TRUE, 0.75),
       quartile.1=apply(mat, 2, quantile, na.rm=TRUE, 1),
       mean=apply(mat, 2, mean, na.rm=TRUE),
       sd=apply(mat, 2, sd, na.rm=TRUE),
       mad=apply(mat, 2, mad, na.rm=TRUE))
   ## suppressWarnings({
   ##   c(x,
   ##   log=log(x),
   ##   `log+1`=log(x+1),
   ##   log.log=log(log(x)))
   ## })
 },
 ...
### Other parameters to pass to IntervalRegressionMatrix.
 ){
  require(Matrix)
  stopifnot(is.list(problem.list))
  stopifnot(is.numeric(initial.regularization))
  stopifnot(length(initial.regularization) == 1)
  stopifnot(initial.regularization > 0)
  stopifnot(is.numeric(factor.regularization))
  stopifnot(length(factor.regularization) == 1)
  stopifnot(factor.regularization > 1)
  
  n.input.features <- ncol(problem.list[[1]]$features)

  featureVec.list <- list()
  targets.list <- list()
  for(problem.i in seq_along(problem.list)){
    problem <- problem.list[[problem.i]]
    stopifnot(is.matrix(problem$features))
    stopifnot(is.numeric(problem$features))
    stopifnot(ncol(problem$features) == n.input.features)
    stopifnot(is.numeric(problem$target))
    if(length(problem$target) != 2){
      print(problem)
      stop("target should be numeric vector of length 2")
    }
    problem <- problem.list[[problem.i]]
    fvec <- featureVec(problem$features)
    names.tab <- table(names(fvec))
    stopifnot(names.tab == 1)
    featureVec.list[[problem.i]] <- fvec
    targets.list[[problem.i]] <- problem$target
  }
  all.targets <- do.call(rbind, targets.list)
  is.trivial.target <- apply(!is.finite(all.targets), 1, all)
  nontrivial.i.vec <- which(!is.trivial.target)

  all.featureVec.mat <- do.call(rbind, featureVec.list)
  is.finite.feature <- apply(is.finite(all.featureVec.mat), 2, all)
  finite.featureVec.mat <- all.featureVec.mat[, is.finite.feature]
  all.mean.vec <- colMeans(finite.featureVec.mat)
  all.sd.vec <- apply(finite.featureVec.mat, 2, sd)
  is.invariant <- all.sd.vec == 0
  train.feature.i <- which(!is.invariant)
  train.feature.names <- colnames(finite.featureVec.mat)[train.feature.i]
  mean.vec <- all.mean.vec[train.feature.names]
  sd.vec <- all.sd.vec[train.feature.names]

  raw.featureVec.mat <- finite.featureVec.mat[, train.feature.names]
  train.mean.mat <-
    matrix(mean.vec, nrow(raw.featureVec.mat), ncol(raw.featureVec.mat),
           byrow=TRUE)
  train.sd.mat <-
    matrix(sd.vec, nrow(raw.featureVec.mat), ncol(raw.featureVec.mat),
           byrow=TRUE)
  norm.featureVec.mat <- (raw.featureVec.mat-train.mean.mat)/train.sd.mat
  intercept.mat <-
    Matrix(cbind("(Intercept)"=1, norm.featureVec.mat))[nontrivial.i.vec,]

  targets.mat <- Matrix(do.call(rbind, targets.list[nontrivial.i.vec]))
  rownames(intercept.mat) <- rownames(targets.mat) <-
    names(problem.list)[nontrivial.i.vec]

  regularization <- initial.regularization
  n.weights <- ncol(intercept.mat)
  weight.vec <- Matrix(rep(0, n.weights))
  n.nonzero <- n.weights
  
  weight.vec.list <- list()
  regularization.vec.list <- list()
  while(n.nonzero > 1){
    weight.vec <-
      IntervalRegressionMatrix(intercept.mat, targets.mat,
                               weight.vec,
                               regularization,
                               verbose=verbose,
                               threshold=threshold,
                               ...)
    n.zero <- sum(weight.vec == 0)
    n.nonzero <- sum(weight.vec != 0)
    l1.norm <- sum(abs(weight.vec[-1]))
    if(verbose >= 1){
      cat(sprintf("regularization=%8.4f L1norm=%8.4f zeros=%d\n",
                  regularization, l1.norm, n.zero))
    }
    weight.vec.list[[paste(regularization)]] <- weight.vec
    regularization.vec.list[[paste(regularization)]] <- regularization
    regularization <- regularization * factor.regularization
  }
  train.weight.mat <- Matrix(do.call(cbind, weight.vec.list))
  regularization.vec <- do.call(c, regularization.vec.list)
  dimnames(train.weight.mat) <-
    list(feature=colnames(intercept.mat),
         regularization=regularization.vec)
  can.ignore <- apply(train.weight.mat==0, 1, all)
  pred.weight.mat <- t(train.weight.mat[!can.ignore, ])
  pred.feature.names <- colnames(pred.weight.mat)[-1]
  pred.mean.vec <- mean.vec[pred.feature.names]
  pred.sd.vec <- sd.vec[pred.feature.names]
  list(pred.weight.mat=pred.weight.mat,
       featureVec=featureVec,
       regularization.vec=regularization.vec,
       pred.mean.vec=pred.mean.vec,
       pred.sd.vec=pred.sd.vec,
       pred.feature.names=pred.feature.names,
       predict=function(mat){
         stopifnot(is.matrix(mat))
         stopifnot(is.numeric(mat))
         fvec <- featureVec(mat)
         stopifnot(pred.feature.names %in% names(fvec))
         raw.vec <- fvec[pred.feature.names]
         raw.vec[!is.finite(raw.vec)] <- 0 
         norm.vec <- (raw.vec-pred.mean.vec)/pred.sd.vec
         intercept.vec <- c("(Intercept)"=1, norm.vec)
         pred.weight.mat %*% intercept.vec
       })
### List representing fit model. You can do
### fit$predict(feature.matrix) to get a predicted log penalty
### value. The pred.mean.vec and pred.sd.vec were used for scaling the
### training data matrices. The pred.weight.mat is the
### n.regularization * n.features numeric matrix of optimal
### coefficients.
}, ex=function(){
  library(PeakSegJoint)
  data(H3K4me3.PGP.immune.4608)
  chrom.vec <- sub(":.*", "", names(H3K4me3.PGP.immune.4608))
  table(chrom.vec)
  train.chroms <- c("chr1", "chr9")
  sets <-
    list(train=chrom.vec %in% train.chroms,
         validation=! chrom.vec %in% train.chroms)
  train.problems <- H3K4me3.PGP.immune.4608[sets$train]
  fit.list <-
    list(nonlinear=NonlinearProblems(train.problems),
         linear=IntervalRegressionProblems(train.problems))
  set.error.list <- list()
  for(set.name in names(sets)){
    in.set <- sets[[set.name]]
    problem.list <- H3K4me3.PGP.immune.4608[in.set]
    for(algorithm in names(fit.list)){
      fit <- fit.list[[algorithm]]
      error.mat.list <- list()
      for(problem.name in names(problem.list)){
        problem <- problem.list[[problem.name]]
        pred.log.lambda <- fit$predict(problem$features)
        too.hi <- problem$target[2] < pred.log.lambda
        too.lo <- pred.log.lambda < problem$target[1]
        is.error <- too.hi | too.lo
        error.mat.list[[problem.name]] <- is.error
      }
      error.mat <- do.call(cbind, error.mat.list)
      percent.error <- rowMeans(error.mat) * 100
      set.error.list[[paste(set.name, algorithm)]] <-
        data.frame(set.name,
                   algorithm,
                   regularization=fit$regularization,
                   percent.error)
    }
  }
  set.error <- do.call(rbind, set.error.list)
  library(ggplot2)
  ggplot()+
    geom_line(aes(-log10(regularization), percent.error,
                  group=interaction(algorithm, set.name),
                  color=algorithm,
                  linetype=set.name),
              data=set.error)
})

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
      if(!is.null(problem)){
        if(problem.name %in% names(chunk.env$step2.error.list)){
          plist <- chunk.env$step2.error.list[[problem.name]]$problem
          not.null <- !sapply(plist, is.null)
          data.names <- names(plist)[not.null]
          problem[data.names] <- plist
        }
        data.by.problem[[problem.name]] <- FUN(problem)
      }#!is.null(problem.name
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
      NonlinearProblems(train.list,
                        initial.regularization=0.005,
                        factor.regularization=1.1,
                        verbose=0)
    error.by.tv <- list()
    for(tv in names(sets)){
      error.vec.list <- labeledProblems(sets[[tv]], function(problem){
        log.lambda.mat <- fit$predict(problem$features)
        error.vec <- rep(NA, length(log.lambda.mat))
        for(log.lambda.i in seq_along(error.vec)){
          log.lambda <- log.lambda.mat[log.lambda.i, ]
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

nonlinear.error.all.list <- list()
nonlinear.error.list <- list()
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
    fit <- NonlinearProblems(train.list,
                             initial.regularization=reg.est)
    
    is.train <- set.chunks %in% train.chunk.vec
    test.chunks <- set.chunks[!is.train]

    ## This method of computing the error under-estimates the false
    ## positive rate, since it does not include any of the unlabeled
    ## problems.
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

    ## Compute the error on all the problems.
    for(test.chunk in test.chunks){
      peaks.by.problem <- problemData(test.chunk, function(problem){
        if(is.null(problem$peaks))return(NULL)
        log.lambda <- fit$predict(problem$features)[1]
        selected <- 
          subset(problem$modelSelection,
                 min.log.lambda < log.lambda &
                   log.lambda < max.log.lambda)
        stopifnot(nrow(selected) == 1)
        subset(problem$peaks, peaks == selected$peaks)
      })
      pred.peaks <- if(length(peaks.by.problem)==0){
        Peaks()
      }else{
        do.call(rbind, peaks.by.problem)
      }
      regions.RData <-
        sprintf("PeakSegJoint-chunks/%s/regions.RData", test.chunk)
      test.regions <- regions.file.list[[regions.RData]]
      chunk.error <- PeakErrorSamples(pred.peaks, test.regions)
      nonlinear.error.all.list[[paste(split.name, test.chunk)]] <- 
        data.frame(set.name, split.i, test.chunk, chunk.error)
    }
    
    nonlinear.error.list[[split.name]] <-
      data.frame(set.name, split.i, split.error)
  }#split.i
}#set.name

nonlinear.error <- do.call(rbind, nonlinear.error.list)
nonlinear.error.all <- do.call(rbind, nonlinear.error.all.list)

save(nonlinear.error,
     nonlinear.error.all,
     file="nonlinear.error.RData")
