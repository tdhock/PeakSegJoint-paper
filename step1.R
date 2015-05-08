works_with_R("3.2.0",
             data.table="1.9.4",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@73ef55d75ba229e3cfac88b2d7c3d954cb1c0f84")

load("train.sets.RData")
load("chunk.problems.RData")

regions.file.list <- list()
regions.file.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/regions.RData")
for(regions.file in regions.file.vec){
  load(regions.file)
  regions.dt <- data.table(regions)
  setkey(regions.dt, chromStart, chromEnd)
  regions.file.list[[regions.file]] <- regions.dt
}

counts.file.list <- list()
counts.file.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/counts.RData")
for(counts.file in counts.file.vec){
  load(counts.file)
  counts.dt <- data.table(counts)
  counts.dt[, count := coverage ]
  setkey(counts.dt, chromStart, chromEnd)
  counts.file.list[[counts.file]] <- counts.dt
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

step1 <- list()
set.name <- "srf_NTNU_Gm12878"
split.i <- 1
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  data.by.chunk <- chunk.problems[[set.name]]
  for(split.i in seq_along(splits)){
    split.name <- paste(set.name, "split", split.i)
    under.name <- gsub(" ", "_", split.name)

    RData.file <- sprintf("step1/%s.RData", under.name)
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
      best.res <- res.stats[which.min(res.stats$errors), ]
      res.str <- paste(best.res$res.str)

      ## Cluster peaks, fit another segmentation model, and take its
      ## ann error.
      err.by.res <- list()
      for(chunk.i in seq_along(train.validation)){
        chunk.name <- train.validation[[chunk.i]]
        cat(sprintf("%4d / %4d chunks\n", chunk.i, length(train.validation)))
        data.by.res <- data.by.chunk[[chunk.name]]

        counts.file <-
          sprintf("../chip-seq-paper/chunks/%s/counts.RData",
                  chunk.name)
        counts <- counts.file.list[[counts.file]]
        stopifnot(is.data.frame(counts))
        stopifnot(nrow(counts) > 0)

        regions.file <-
          sprintf("../chip-seq-paper/chunks/%s/regions.RData",
                  chunk.name)
        regions <- regions.file.list[[regions.file]]
        stopifnot(is.data.frame(regions))
        stopifnot(nrow(regions) > 0)
        regions.by.sample <- split(regions, regions$sample.id)
        regions.by.chromStart <- split(regions, regions$chromStart)
        chrom <- paste(regions$chrom[1])
        
        for(res.str in names(data.by.res)){
          bases.per.problem <- as.integer(res.str)
          data.by.problem <- data.by.res[[res.str]]
          peaks.by.problem <- list()
          for(problem.name in names(data.by.problem)){
            problem <- data.by.problem[[problem.name]]
            selected <- subset(problem$modelSelection, peaks==max(peaks))
            stopifnot(nrow(selected) == 1)
            if(selected$peaks > 0){
              peaks.str <- paste(selected$peaks)
              peaks.df <- problem$peaks[[peaks.str]]
              peaks.by.problem[[problem.name]] <- peaks.df[1,]
            }
          }#problem.name
          step1.peaks <- do.call(rbind, peaks.by.problem)

          problem.list <- list()
          if(nrow(step1.peaks) > 0){
            clustered.peaks <- clusterPeaks(step1.peaks)
            peaks.by.cluster <- split(clustered.peaks, clustered.peaks$cluster)
            pred.by.cluster <- list()
            for(cluster.name in names(peaks.by.cluster)){
              cluster <- peaks.by.cluster[[cluster.name]]
              merged.peak <- with(cluster, {
                data.frame(chromStart=min(chromStart),
                           chromEnd=max(chromEnd))
              })
              pred.by.cluster[[cluster.name]] <-
                data.frame(merged.peak,
                           sample.id=unique(cluster$sample.id))
              cluster.chromStart <- min(cluster$chromStart)
              cluster.chromEnd <- max(cluster$chromEnd)
              cluster.mid <-
                as.integer((cluster.chromEnd + cluster.chromStart)/2)
              half.bases <- as.integer(bases.per.problem/2)
              cluster.num <- as.numeric(cluster.name)
              before.name <- paste(cluster.num-1)
              chromEnd.before <- if(before.name %in% names(peaks.by.cluster)){
                max(peaks.by.cluster[[before.name]]$chromEnd)
              }else{
                0
              }
              after.name <- paste(cluster.num+1)
              chromStart.after <- if(after.name %in% names(peaks.by.cluster)){
                min(peaks.by.cluster[[after.name]]$chromStart)
              }else{
                Inf
              }
              problemStart <- as.integer(cluster.chromStart - half.bases)
              if(problemStart < chromEnd.before){
                problemStart <-
                  as.integer((chromEnd.before+cluster.chromStart)/2)
              }
              problemEnd <- as.integer(cluster.chromEnd + half.bases)
              if(chromStart.after < problemEnd){
                problemEnd <- as.integer((chromStart.after+cluster.chromEnd)/2)
              }
              stopifnot(problemStart <= cluster.chromStart)
              stopifnot(cluster.chromEnd <= problemEnd)
              problem.i <- as.numeric(cluster.name)+1
              problem.list[[problem.i]] <-
                data.frame(problem.i,
                           problem.name=sprintf("%s:%d-%d",
                             chrom, problemStart, problemEnd),
                           problemStart, problemEnd,
                           peakStart=merged.peak$chromStart,
                           peakEnd=merged.peak$chromEnd)
            }
          }
          step2.overlap <- do.call(rbind, problem.list)
          step2.problems <- with(step2.overlap, {
            prev.problemEnd <- problemEnd[-length(problemEnd)]
            next.problemStart <- problemStart[-1]
            overlaps.next <- which(next.problemStart <= prev.problemEnd)
            mid <- as.integer((prev.problemEnd+next.problemStart)/2)
            problemEnd[overlaps.next] <- mid[overlaps.next]
            problemStart[overlaps.next+1] <- mid[overlaps.next]+1L
            data.frame(problem.i=seq_along(problemStart),
                       problem.name=sprintf("%s:%d-%d",
                         chrom, problemStart, problemEnd),
                       problemStart, problemEnd)
          })
          stopifnot(with(step2.problems, {
            problemEnd[-length(problemEnd)] < problemStart[-1]
          }))
          
          ggplot()+
            geom_segment(aes(problemStart/1e3, problem.i,
                             color=what,
                             xend=problemEnd/1e3, yend=problem.i),
                         data=data.frame(step2.overlap, what="overlap"))+
            geom_segment(aes(problemStart/1e3, problem.i,
                             color=what,
                             xend=problemEnd/1e3, yend=problem.i),
                         data=data.frame(step2.problems, what="corrected"))

          problems.dt <- data.table(step2.problems)
          setkey(problems.dt, problemStart, problemEnd)
          over.regions <- foverlaps(regions, problems.dt, nomatch=0L)
          regions.by.problem <-
            split(over.regions, over.regions$problem.name, drop=TRUE)
          peaks.by.problem <- list()
          for(problem.name in names(regions.by.problem)){
            problem.regions <- regions.by.problem[[problem.name]]
            problem.i <- problem.regions$problem.i[1]
            problem <- problems.dt[problem.i, ]
            problem.counts <-
              foverlaps(counts, problem, nomatch=0L, type="within")
            
            tryCatch({
              fit <- PeakSegJointHeuristic(problem.counts)
              converted <- ConvertModelList(fit)
              prob.err.list <- PeakSegJointError(converted, problem.regions)
              best.models <-
                subset(prob.err.list$modelSelection, errors==min(errors))
              peaks.num <- min(best.models$peaks)
              if(peaks.num > 0){
                show.peaks <- subset(converted$peaks, peaks == peaks.num)
                peaks.by.problem[[problem.name]] <- show.peaks
              }
            }, error=function(e){
              print(e)
            })
          }#problem.name
          pred.peaks <- do.call(rbind, peaks.by.problem)
          ggplot()+
            geom_point(aes(chromStart/1e3, sample.id),
                       data=pred.peaks,
                       pch=1)+
            geom_segment(aes(chromStart/1e3, sample.id,
                             xend=chromEnd/1e3, yend=sample.id),
                         data=pred.peaks)
          error.regions <- PeakErrorSamples(pred.peaks, regions)
          err.by.res[[res.str]][[chunk.name]] <- 
            with(error.regions, {
              data.frame(fp=sum(fp),
                         fn=sum(fn),
                         errors=sum(fp+fn),
                         regions=length(fp))
            })
        }#res.str
      }#chunk.i
      stats.by.res <- list()
      for(res.str in names(err.by.res)){
        res.totals <- do.call(rbind, err.by.res[[res.str]])
        stats.by.res[[res.str]] <- with(res.totals, {
          data.frame(res.str,
                     fp=sum(fp),
                     fn=sum(fn),
                     errors=sum(fp+fn),
                     regions=sum(regions))
        })
      }
      step2.stats <- do.call(rbind, stats.by.res)
      print(step2.stats)
      print(res.stats)
      ##best.res <- res.stats[which.min(res.stats$errors), ]
      ##res.str <- paste(best.res$res.str)
      
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
          IntervalRegressionProblems(train.data.list)
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
