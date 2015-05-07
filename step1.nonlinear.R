works_with_R("3.2.0",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@ff5a7c58e297b54b328047f4e02285f0cb5d2838")

load("train.sets.RData")
load("chunk.problems.RData")

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
       quartile=apply(mat, 2, quantile, na.rm=TRUE),
       mean=apply(mat, 2, mean, na.rm=TRUE),
       sd=apply(mat, 2, sd, na.rm=TRUE),
       mad=apply(mat, 2, mad, na.rm=TRUE))
   suppressWarnings({
     c(x,
     log=log(x),
     `log+1`=log(x+1),
     log.log=log(log(x)))
   })
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
    stopifnot(table(names(fvec)) == 1)
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
  intercept.mat <- Matrix(cbind("(Intercept)"=1, norm.featureVec.mat))

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
    list(nonlinear=NonlinearProblems(train.problems)
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

step1 <- list()
set.name <- "srf_NTNU_Gm12878"
split.i <- 1
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  data.by.chunk <- chunk.problems[[set.name]]
  for(split.i in seq_along(splits)){
    split.name <- paste(set.name, "split", split.i)
    under.name <- gsub(" ", "_", split.name)

    RData.file <- sprintf("step1.nonlinear/%s.RData", under.name)
    if(file.exists(RData.file)){
      load(RData.file)
    }else{
      print(split.name)
      train.validation <- splits[[split.i]]

      ## Use minimum incorrect regions to select the res.str bases per
      ## problem hyper-parameter.
      err.by.res <- list()
      for(chunk.name in train.validation){
        data.by.res <- data.by.chunk[[chunk.name]]
        for(res.str in names(data.by.res)){
          data.by.problem <- data.by.res[[res.str]]
          peaks.by.problem <- list()
          for(problem.name in names(data.by.problem)){
            problem <- data.by.problem[[problem.name]]
            error.range <- range(problem$modelSelection$weighted.error)
            selected <- if(diff(error.range) == 0){
              subset(problem$modelSelection, peaks==0)
            }else{
              min.err <-
                subset(problem$modelSelection,
                       weighted.error==min(weighted.error))
              subset(min.err, peaks == min(peaks))
            }
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
            err.by.res[[res.str]][[paste(chunk.name, sample.id)]] <- 
              data.frame(chunk.name, sample.id, error.regions)
          }#sample.id
        }#res.str
      }#chunk.name
      stats.by.res <- list()
      for(res.str in names(err.by.res)){
        error.regions <- do.call(rbind, err.by.res[[res.str]])
        fp <- sum(error.regions$fp)
        fn <- sum(error.regions$fn)
        errors <- fp+fn
        stats.by.res[[res.str]] <- 
          data.frame(res.str, fp, fn, errors, regions=nrow(error.regions))
      }
      res.stats <- do.call(rbind, stats.by.res)
      print(res.stats)
      stop(1)
      best.res <- res.stats[which.min(res.stats$errors), ]
      res.str <- paste(best.res$res.str)

      ## Determine train/validation splits.
      n.folds <- min(length(train.validation), 4)
      set.seed(1)
      fold.id <- sample(rep(1:n.folds, l=length(train.validation)))
      fold.regularization.list <- list()
      for(fold.i in 1:n.folds){
        is.validation <- fold.i == fold.id
        sets <-
          list(train=train.validation[!is.validation],
               validation=train.validation[is.validation])
        train.data.list <- list()
        for(train.chunk in sets$train){
          data.by.problem <- data.by.chunk[[train.chunk]][[res.str]]
          for(problem.name in names(data.by.problem)){
            problem <- data.by.problem[[problem.name]]
            train.data.list[[paste(train.chunk, problem.name)]] <- problem
          }
        }
        fit <-
          NonlinearProblems(train.data.list)
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
            ggtitle(paste("fold", fold.i))+
            geom_point(aes(-log10(regularization), errors,
                           color=tv),
                       size=3,
                       pch=1,
                       data=v.min.err)+
            geom_line(aes(-log10(regularization), errors,
                          color=tv, group=tv),
                      data=tv.reg.err)+
            geom_point(aes(-log10(regularization), errors,
                           color=tv, group=tv),
                       size=1,
                       data=tv.reg.err)+
            xlab("model complexity -log10(regularization)")
        print(resPlot)

        fold.regularization.list[[fold.i]] <-
          mean(v.min.err$regularization)
      }#fold.i

      reg.num <- mean(unlist(fold.regularization.list))

      train.data.list <- list()
      for(chunk.name in train.validation){
        data.by.problem <- data.by.chunk[[chunk.name]][[res.str]]
        for(problem.name in names(data.by.problem)){
          problem <- data.by.problem[[problem.name]]
          train.data.list[[paste(chunk.name, problem.name)]] <- problem
        }
      }
      fit <-
        IntervalRegressionProblems(train.data.list,
                                   initial.regularization=reg.num)
      reg.str <- colnames(fit$weight.mat)[1]
      
      split.data <-
        list(res.str=res.str,
             reg.str=reg.str,
             res.stats=res.stats,
             best.res=best.res,
             fit=fit)
      
      RData.dir <- dirname(RData.file)
      dir.create(RData.dir, showWarnings=FALSE, recursive=TRUE)
      save(split.data, file=RData.file)
    }#if(file.exists(RData.file)) else
    step1[[split.name]] <- split.data
  }#split.i
}#set.name

save(step1, file="step1.RData")
