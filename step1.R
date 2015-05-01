works_with_R("3.2.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@c68c566a606aea75d646e81087d168e5e2f0531a")

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
 verbose=1,
 ...
### Other parameters to pass to IntervalRegressionMatrix.
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
    problem <- problem.list[[problem.i]]
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
  regularization.vec.list <- list()
  while(n.nonzero > 1){
    weight.vec <-
      IntervalRegressionMatrix(norm.featureSum.mat, targets.mat,
                               weight.vec,
                               regularization,
                               verbose=verbose,
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
  weight.mat <- do.call(cbind, weight.vec.list)
  list(weight.mat=weight.mat,
       regularization.vec=do.call(c, regularization.vec.list),
       mean.vec=mean.vec,
       sd.vec=sd.vec,
       train.feature.names=train.feature.names,
       predict=function(mat){
         stopifnot(is.matrix(mat))
         stopifnot(is.numeric(mat))
         stopifnot(train.feature.names %in% colnames(mat))
         raw.mat <- mat[, train.feature.names, drop=FALSE]
         raw.mat[!is.finite(raw.mat)] <- 0 
         mean.mat <- matrix(mean.vec, nrow(raw.mat), ncol(raw.mat), byrow=TRUE)
         sd.mat <- matrix(sd.vec, nrow(raw.mat), ncol(raw.mat), byrow=TRUE)
         norm.mat <- (raw.mat-mean.mat)/sd.mat
         intercept.mat <- cbind("(Intercept)"=1, norm.mat)
         colSums(intercept.mat %*% weight.mat)
       })
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
### Print messages if > 0.
 ){
  stopifnot(is.matrix(features))
  stopifnot(is.numeric(features))
  n.features <- ncol(features)
  n.problems <- nrow(features)

  stopifnot(is.matrix(targets))
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

  stopifnot(is.numeric(initial.weight.vec))
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
}

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

set.seed(1)

all.split.list <- list()

set.name <- "srf_NTNU_Gm12878"
split.i <- 1
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  data.by.chunk <- chunk.problems[[set.name]]
  for(split.i in seq_along(splits)){
    split.name <- paste(set.name, "split", split.i)
    train.validation <- splits[[split.i]]
    train.i <- sample(length(train.validation), length(train.validation)/2)
    sets <-
      list(train=train.validation[train.i],
           validation=train.validation[-train.i])
    model.error.list <- list() 
    fit.list <- list()
    pred.peak.list <- list()
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
      fit.list[[res.str]] <- fit
      tv.err.list <- list()
      for(tv in names(sets)){
        set.chunks <- sets[[tv]]
        for(chunk.name in set.chunks){
          pred.list <- chunkError(fit, chunk.name)
          pred.peak.list[[res.str]][[chunk.name]] <- pred.list$peak.list
          tv.err.list[[paste(tv, chunk.name)]] <-
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
      model.error.list[[res.str]] <-
        data.frame(bases.per.problem=as.numeric(res.str),
                   tv.reg.err)
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
    all.split.list[[split.name]] <- 
      list(glob.min.err=glob.min.err,
           res.str=res.str,
           reg.str=reg.str,
           fit=fit.list[[res.str]],
           peak.list=out.peak.list)
  }#split.i
}#set.name
