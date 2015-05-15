works_with_R("3.2.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@c91c0deafc39321023307926a289c6acf5d0ced4",
             data.table="1.9.4",
             ggplot2="1.0")

## The script computes a list of PeakSegJoint problem regions for each
## data set and bases.per.problem size.

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

genome.pos.pattern <-
  paste0("(?<chrom>chr.*?)",
         ":",
         "(?<chromStart>[0-9]+)",
         "-",
         "(?<chromEnd>[0-9]+)")

tf.size <- 4.5 * 2^(4:13)

target.sizes <-
  list(H3K36me3=4.5 * 2^(12:20),
       H3K4me3=4.5 * 2^(6:17),
       nrsf=tf.size,
       srf=tf.size,
       max=tf.size)

bases.per.problem.all <- 4.5 * 2^(4:20)

set.dir.i <- 5
chunk.id <- "13"
bases.per.problem <- 144
bases.per.problem <- 2304
bases.per.problem <- 36864

set.dirs <- Sys.glob("../chip-seq-paper/chunks/*_*_*") #include TF data!
for(set.dir.i in seq_along(set.dirs)){
  set.dir <- set.dirs[[set.dir.i]]
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  bases.per.problem.vec <- target.sizes[[experiment]]
  chunk.ids <- dir(set.dir)
  for(chunk.id in chunk.ids){
    chunk.path <- file.path(set.dir, chunk.id)
    chunk.name <- paste0(set.name, "/", chunk.id)
    counts.RData <- file.path(chunk.path, "counts.RData")
    load(counts.RData)
    regions.RData <- file.path(chunk.path, "regions.RData")
    load(regions.RData)
    regions.dt <- data.table(regions)
    counts.dt <- data.table(counts)
    counts.dt[, count := coverage]
    setkey(regions.dt, chromStart, chromEnd)
    setkey(counts.dt, chromStart, chromEnd)
    regions.dt$region.i <- 1:nrow(regions.dt)
    chrom <- paste(regions$chrom[1])
    max.chromEnd <- max(counts$chromEnd)
    min.chromStart <- min(counts$chromStart)
    for(bases.per.problem in bases.per.problem.vec){
      res.str <- paste(bases.per.problem)
      RData.file <-
        sprintf("chunk.problems/%s/%s/%s.RData",
                set.name, chunk.id, bases.per.problem)
      if(file.exists(RData.file)){
        obj.names <- load(RData.file)
      }else{
        cat(sprintf("%s %s\n", chunk.name, bases.per.problem))
        problemSeq <- seq(0, max.chromEnd, by=bases.per.problem)
        problemStart <-
          as.integer(sort(c(problemSeq,
                            problemSeq+bases.per.problem/2)))
        problemEnd <- problemStart+bases.per.problem
        is.overlap <- min.chromStart < problemEnd &
          problemStart < max.chromEnd
        problem.name <- sprintf("%s:%d-%d", chrom, problemStart, problemEnd)
        problems <- 
          data.frame(problem.name, 
                     bases.per.problem, problemStart, problemEnd)[is.overlap,]
        all.problems.dt <- data.table(problems)
        setkey(all.problems.dt, problemStart, problemEnd)

        overlap.regions <- foverlaps(regions.dt, all.problems.dt, nomatch=0L)
        region.counts <- table(overlap.regions$region.i)
        region.weights <- 1/region.counts
        overlap.regions[, weight := region.weights[paste(region.i)]]
        table(paste(overlap.regions$problem.name))

        problems.with.regions <- paste(unique(overlap.regions$problem.name))
        setkey(overlap.regions, problem.name)

        step1.by.problem <- list() #[[problem.name]]
        problem.stats.list <- list() #[[problem.name]]
        for(problem.i in 1:nrow(all.problems.dt)){
          problem <- all.problems.dt[problem.i, ]
          problem.name <- paste(problem$problem.name)
          problem.regions <- overlap.regions[problem.name]
          regions.by.sample <-
            split(problem.regions, problem.regions$sample.id, drop=TRUE)
          ## microbenchmark(foverlaps={
          ##   problem.counts <- foverlaps(counts.dt, problem, nomatch=0L)
          ##   print(nrow(problem.counts))
          ## },vector.ops={
          ##   problem.counts <-
          ##     counts.dt[! (chromEnd < problem$problemStart |
          ##                 problem$problemEnd < chromStart), ]
          ##   print(nrow(problem.counts))
          ## }, times=10)
          problem.counts <-
            counts.dt[! (chromEnd < problem$problemStart |
                           problem$problemEnd < chromStart), ]
          profile.list <- ProfileList(problem.counts)
          converted <- tryCatch({
            fit <- PeakSegJointHeuristic(profile.list)
            ConvertModelList(fit)
          }, error=function(e){
            list(loss=data.frame(peaks=0, loss=0))
          })
          peaks.by.peaks <- list("0"=Peaks())
          if(!is.null(converted$peaks)){
            peaks.list <- with(converted, split(peaks, peaks$peaks))
            peaks.by.peaks[names(peaks.list)] <- peaks.list
          }
          error.by.peaks <- list()
          if(length(regions.by.sample))for(peaks.str in names(peaks.by.peaks)){
            model.peaks <- peaks.by.peaks[[peaks.str]]
            peaks.by.sample <- split(model.peaks, model.peaks$sample.id)
            error.by.sample <- list()
            for(sample.id in names(regions.by.sample)){
              sample.regions <- regions.by.sample[[sample.id]]
              sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
                peaks.by.sample[[sample.id]]
              }else{
                Peaks()
              }
              sample.error <- PeakErrorChrom(sample.peaks, sample.regions)
              sample.error$weight <- sample.regions$weight
              sample.error$weighted.error <- with(sample.error, (fp+fn)*weight)
              error.by.sample[[sample.id]] <-
                data.frame(sample.id, peaks=as.numeric(peaks.str), sample.error)
            }#sample.id
            error.by.peaks[[peaks.str]] <- do.call(rbind, error.by.sample)
          }#peaks.str
          loss.df <- converted$loss
          rownames(loss.df) <- loss.df$peaks
          if(length(regions.by.sample)){
            loss.df[names(error.by.peaks), "weighted.error"] <-
              sapply(error.by.peaks, with, sum(weighted.error))
            loss.df[names(error.by.peaks), "total.weight"] <-
              sapply(error.by.peaks, with, sum(weight))
          }else{
            loss.df$weighted.error <- 0
            loss.df$total.weight <- 0
          }
          loss.df$cummin <- cummin(loss.df$loss)
          some.loss <- subset(loss.df, loss==cummin & !is.na(total.weight))
          exact <-
            with(some.loss, exactModelSelection(loss, peaks, peaks))
          exact$weighted.error <- loss.df[paste(exact$peaks), "weighted.error"]
          stopifnot(!is.na(exact$weighted.error))
          indices <- with(exact, {
            largestContinuousMinimum(weighted.error,
                                     max.log.lambda-min.log.lambda)
          })
          total.weight <- some.loss$total.weight[1]
          stopifnot(total.weight == some.loss$total.weight)
          problem.stats.list[[problem.name]] <-
            data.frame(problem.name,
                       weighted.error=min(some.loss$weighted.error),
                       total.weight)
          step1.by.problem[[problem.name]] <- 
            list(features=featureMatrix(profile.list),
                 modelSelection=exact,
                 peaks=peaks.by.peaks,
                 target=c(min.log.lambda=exact$min.log.lambda[indices$start],
                   max.log.lambda=exact$max.log.lambda[indices$end]))
        }#problem.i
        problem.stats <- do.call(rbind, problem.stats.list)
        total.weight <- sum(problem.stats$total.weight)
        stopifnot(all.equal(nrow(regions), total.weight))
        error.row <- 
          data.frame(set.name, chunk.name, bases.per.problem,
                     weighted.error=sum(problem.stats$weighted.error),
                     total.weight)
        ## Also compute step2 data.
        peaks.by.problem <- list()
        for(problem.i in seq_along(step1.by.problem)){
          problem.name <- names(step1.by.problem)[[problem.i]]
          pos.row <- str_match_perl(problem.name, genome.pos.pattern)
          problem <- step1.by.problem[[problem.name]]
          selected <- subset(problem$modelSelection, peaks==max(peaks))
          stopifnot(nrow(selected) == 1)
          if(selected$peaks > 0){
            peaks.str <- paste(selected$peaks)
            peaks.df <- problem$peaks[[peaks.str]]
            peaks.row <-
              data.frame(problem.i, problem.name,
                         problemStart=as.integer(pos.row[, "chromStart"]),
                         problemEnd=as.integer(pos.row[, "chromEnd"]),
                         peaks.df[1,])
            peaks.row$sample.id <- "step 1"
            peaks.by.problem[[problem.name]] <- peaks.row
          }
        }#problem.name

        probPlot <- 
        ggplot()+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3),
                        data=data.frame(chromStart=min.chromStart,
                          chromEnd=max.chromEnd),
                        fill="grey")+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                            fill=annotation),
                        data=regions)+
          scale_fill_manual(values=ann.colors)+
          geom_segment(aes(problemStart/1e3, problem.name,
                           xend=problemEnd/1e3, yend=problem.name),
                       data=problems)

        step1.peaks <- do.call(rbind, peaks.by.problem)
        
        peakPlot <-
          probPlot+
            geom_segment(aes(chromStart/1e3, problem.name,
                             xend=chromEnd/1e3, yend=problem.name),
                         size=2,
                         data=step1.peaks)

        problemPlot <- 
          ggplot()+
            scale_y_continuous("aligned read coverage",
                               breaks=function(limits){
                                 floor(limits[2])
                               })+
            scale_linetype_manual("error type",
                                  limits=c("correct", 
                                    "false negative",
                                    "false positive"
                                           ),
                                  values=c(correct=0,
                                    "false negative"=3,
                                    "false positive"=1))+
            xlab("position on chromosome (kilo bases = kb)")+
            geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                              fill=annotation),
                          alpha=0.5,
                          color="grey",
                          data=regions)+
            scale_fill_manual(values=ann.colors)+
            theme_bw()+
            theme(panel.margin=grid::unit(0, "cm"))+
            facet_grid(sample.id ~ ., labeller=function(var, val){
              sub("McGill0", "", sub(" ", "\n", val))
            }, scales="free")+
            geom_step(aes(chromStart/1e3, coverage),
                      data=counts,
                      color="grey50")+
            geom_segment(aes(problemStart/1e3, problem.i,
                             xend=problemEnd/1e3, yend=problem.i),
                         data=step1.peaks)+
            geom_segment(aes(chromStart/1e3, problem.i,
                             xend=chromEnd/1e3, yend=problem.i),
                         color="deepskyblue",
                         size=2,
                         data=step1.peaks)

        if(length(peaks.by.problem) == 0){
          pred.peaks <- Peaks()
          step2.by.problem <- list()
        }else{
          step2.overlap.list <- list()
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
            step2.overlap.list[[problem.i]] <-
              data.frame(problem.i,
                         problem.name=sprintf("%s:%d-%d",
                           chrom, problemStart, problemEnd),
                         problemStart, problemEnd,
                         peakStart=merged.peak$chromStart,
                         peakEnd=merged.peak$chromEnd)
          }
          
          step2.overlap <- do.call(rbind, step2.overlap.list)
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

          step2.prob.plot <- 
          problemPlot+
            geom_segment(aes(problemStart/1e3, problem.i,
                             xend=problemEnd/1e3, yend=problem.i),
                         data=data.frame(step2.problems, sample.id="step 2"))

          problems.dt <- data.table(step2.problems)
          setkey(problems.dt, problemStart, problemEnd)
          setkey(regions.dt, chromStart, chromEnd)
          over.regions <- foverlaps(regions.dt, problems.dt, nomatch=0L)
          over.regions[,
            `:=`(overlapStart=ifelse(problemStart < chromStart,
                   chromStart, problemStart),
                 overlapEnd=ifelse(problemEnd < chromEnd,
                   problemEnd, chromEnd))]
          over.regions[, overlapBases := overlapEnd-overlapStart]
          region.i.problems <- over.regions[,
            .(problem.name=problem.name[which.max(overlapBases)]),
                       by=region.i]
          stopifnot(nrow(region.i.problems) == nrow(regions.dt))
          setkey(regions.dt, region.i)
          setkey(region.i.problems, region.i)
          assigned.regions <- regions.dt[region.i.problems,]
          stopifnot(nrow(assigned.regions) == nrow(regions.dt))
          regions.by.problem <-
            split(assigned.regions, assigned.regions$problem.name, drop=TRUE)
          setkey(problems.dt, problem.name)
          peaks.by.problem <- list()
          step2.by.problem <- list()
          step2.peak.list <- list()
          saved.problem.list <- list()
          for(problem.name in names(regions.by.problem)){
            problem.regions <- regions.by.problem[[problem.name]]
            problem <- problems.dt[problem.name]
            setkey(problem, problemStart, problemEnd)
            problem.i <- problem$problem.i
            problem.counts <-
              foverlaps(counts.dt, problem, nomatch=0L, type="within")
            
            tryCatch({
              profile.list <- ProfileList(problem.counts)
              fit <- PeakSegJointHeuristic(profile.list)
              converted <- ConvertModelList(fit)
              prob.err.list <- PeakSegJointError(converted, problem.regions)
              step2.by.problem[[problem.name]] <-
                list(converted=converted,
                     error=prob.err.list,
                     features=featureMatrix(profile.list))
              saved.problem.list[[problem.name]] <-
                data.frame(problem, sample.id="step 2")
              best.models <-
                subset(prob.err.list$modelSelection, errors==min(errors))
              peaks.num <- min(best.models$peaks)
              if(peaks.num > 0){
                show.peaks <- subset(converted$peaks, peaks == peaks.num)
                peaks.by.problem[[problem.name]] <- show.peaks
                peak.row <- show.peaks[1,]
                peak.row$sample.id <- "step 2"
                step2.peak.list[[problem.name]] <-
                  data.frame(problem.i, peak.row)
              }
            }, error=function(e){
              paste(e)
            })
          }#problem.name
          pred.peaks <- do.call(rbind, peaks.by.problem)
          step2.peaks <- do.call(rbind, step2.peak.list)
          saved.problems <- do.call(rbind, saved.problem.list)

          ggplot()+
            geom_point(aes(chromStart/1e3, sample.id),
                       data=pred.peaks,
                       pch=1)+
            geom_segment(aes(chromStart/1e3, sample.id,
                             xend=chromEnd/1e3, yend=sample.id),
                         data=pred.peaks)
        }
        if(is.null(pred.peaks))pred.peaks <- Peaks()
        error.regions <- PeakErrorSamples(pred.peaks, regions.dt)

        step2.peak.plot <-
          problemPlot+
            geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                             linetype=status),
                         data=error.regions,
                         fill=NA,
                         size=2,
                         color="black")+
            geom_segment(aes(chromStart/1e3, 0,
                             xend=chromEnd/1e3, yend=0),
                         data=pred.peaks,
                         size=2,
                         color="deepskyblue")+
            geom_segment(aes(problemStart/1e3, problem.i,
                             xend=problemEnd/1e3, yend=problem.i),
                         data=data.frame(saved.problems, sample.id="step 2"))+
            geom_segment(aes(chromStart/1e3, problem.i,
                             xend=chromEnd/1e3, yend=problem.i),
                         data=step2.peaks,
                         size=2,
                         color="deepskyblue")
        
        best.step2.error <- 
          with(error.regions, {
            data.frame(set.name,
                       chunk.name,
                       bases.per.problem,
                       fp=sum(fp),
                       fn=sum(fn),
                       errors=sum(fp+fn),
                       regions=length(fp))
          })
        ## Save results for this chunk/resolution.
        RData.dir <- dirname(RData.file)
        dir.create(RData.dir, showWarnings=FALSE, recursive=TRUE)
        save(step1.by.problem, error.row,
             step2.by.problem, best.step2.error,
             file=RData.file)
      }
    }#bases.per.problem
  }#chunk.id
}#set.dir.i

system("touch chunk.problems.RData")
