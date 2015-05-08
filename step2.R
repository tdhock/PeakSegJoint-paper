works_with_R("3.2.0",
             data.table="1.9.4",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@d7867c03acfbdbb2f19b358bae67bdc374d2432b")

load("step1.RData")
load("chunk.problems.RData")
load("train.sets.RData")

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

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

genome.pos.pattern <-
  paste0("(?<chrom>chr.*?)",
         ":",
         "(?<chromStart>[0-9]+)",
         "-",
         "(?<chromEnd>[0-9]+)")

step2 <- list()
for(set.name in names(chunk.problems)){
  data.by.chunk <- chunk.problems[[set.name]]
  set.chunks <- names(data.by.chunk)
  for(split.i in 1:6){
    train.chunks <- train.sets[[set.name]][[split.i]]
    is.train <- set.chunks %in% train.chunks
    train.validation <- set.chunks[is.train]
    split.name <- paste(set.name, "split", split.i)
    model.info <- step1[[split.name]]
    bases.per.problem <- as.numeric(model.info$res.str)
    chunk.err.list <- list()
    for(chunk.name in train.validation){
      data.by.problem <- data.by.chunk[[chunk.name]][[model.info$res.str]]
      peaks.by.problem <- list()
      step1.peak.list <- list()
      for(problem.i in seq_along(data.by.problem)){
        problem.name <- names(data.by.problem)[[problem.i]]
        problem <- data.by.problem[[problem.name]]
        ## Make a prediction using the fitted model, and show that on
        ## the plot.
        log.lambda.vec <- model.info$fit$predict(problem$features)
        log.lambda <- log.lambda.vec[[model.info$reg.str]]
        selected <- 
          subset(problem$modelSelection,
                 min.log.lambda < log.lambda & log.lambda < max.log.lambda)
        stopifnot(nrow(selected) == 1)
        selected.peaks <- problem$peaks[[paste(selected$peaks)]]
        peaks.by.problem[[problem.name]] <- selected.peaks
        ## For step2, even if there are no predicted peaks, try to get
        ## a peak so we can get problems on the noPeaks regions.
        peak.models <- subset(problem$modelSelection, peaks != 0)
        if(selected$peaks == 0 && nrow(peak.models) > 0){
          selected <- peak.models[which.max(peak.models$peaks),]
        }
        selected.peaks <- problem$peaks[[paste(selected$peaks)]]
        if(nrow(selected.peaks) > 0){
          step1.peak.list[[problem.name]] <-
            data.frame(problem.i, problem.name,
                       sample.id="step 1",
                       selected.peaks[1,])
        }
      }#problem.name
      step1.mat <-
        str_match_perl(names(data.by.problem),
                       genome.pos.pattern)
      chrom <- step1.mat[1, "chrom"]
      step1.problems <-
        data.frame(problemStart=as.integer(step1.mat[,"chromStart"]),
                   problemEnd=as.integer(step1.mat[,"chromEnd"]),
                   sample.id="step 1",
                   problem.i=1:nrow(step1.mat))
      step1.peaks <- do.call(rbind, step1.peak.list)
      problem.list <- list()
      pred.peaks <- if(nrow(step1.peaks) == 0){
        Peaks()
      }else{
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
          cluster.mid <- as.integer((cluster.chromEnd + cluster.chromStart)/2)
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
          ## cat("cluster", cluster.name,
          ##     "before", chromEnd.before,
          ##     "after", chromStart.after,
          ##     after.name, "\n")
          ##problemStart <- cluster.mid - half.bases
          problemStart <- as.integer(cluster.chromStart - half.bases)
          if(problemStart < chromEnd.before){
            ##cat("overlap before cluster", cluster.name, "\n")
            problemStart <- as.integer((chromEnd.before+cluster.chromStart)/2)
          }
          ##problemEnd <- cluster.mid + half.bases
          problemEnd <- as.integer(cluster.chromEnd + half.bases)
          if(chromStart.after < problemEnd){
            ##cat("overlap after cluster", cluster.name, "\n")
            problemEnd <- as.integer((chromStart.after+cluster.chromEnd)/2)
          }
          stopifnot(problemStart < cluster.chromStart)
          stopifnot(cluster.chromEnd < problemEnd)
          problem.i <- as.numeric(cluster.name)+1
          problem.list[[problem.i]] <-
            data.frame(problem.i,
                       problem.name=sprintf("%s:%d-%d",
                         chrom, problemStart, problemEnd),
                       problemStart, problemEnd,
                       peakStart=merged.peak$chromStart,
                       peakEnd=merged.peak$chromEnd)
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
      regions.by.chromStart <- split(regions, regions$chromStart)
      step2.problems <- do.call(rbind, problem.list)
      counts.file <-
        sprintf("../chip-seq-paper/chunks/%s/counts.RData",
                chunk.name)
      counts <- counts.file.list[[counts.file]]
      stopifnot(is.data.frame(counts))
      stopifnot(nrow(counts) > 0)

      show.peaks <- do.call(rbind, peaks.by.problem)

      chunkPlot <- 
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
        geom_step(aes(chromStart/1e3, coverage),
                  data=counts,
                  color="grey50")
      if(nrow(show.peaks) > 0){
        chunkPlot <- chunkPlot+
        geom_segment(aes(chromStart/1e3, 0,
                         xend=chromEnd/1e3, yend=0),
                     size=2,
                     data=show.peaks)
      }
      chunkPlot <- chunkPlot+
        theme_bw()+
        geom_segment(aes(chromStart/1e3, problem.i,
                         xend=chromEnd/1e3, yend=problem.i),
                     size=2,
                     data=step1.peaks)+
        geom_segment(aes(problemStart/1e3, problem.i,
                         xend=problemEnd/1e3, yend=problem.i),
                     data=step1.problems)+
        geom_segment(aes(problemStart/1e3, problem.i,
                         xend=problemEnd/1e3, yend=problem.i),
                     data=data.frame(step2.problems, sample.id="step 2"))+
        geom_text(aes(max(step1.problems$problemEnd)/1e3, 1,
                      label=paste(bases.per.problem, "bases/problem")),
                  vjust=0,
                  hjust=1,
                  data=data.frame(sample.id="step 1", model.info$best.res))+
        theme(panel.margin=grid::unit(0, "cm"))+
        ##coord_cartesian(xlim=(chrom.range+c(-1,1)*chrom.bases*0/8)/1e3)+
        facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
          sub("McGill0", "", sub(" ", "\n", val))
        })+
        ggtitle(paste(split.name, chunk.name))
      chunk.png <- sprintf("step2/%s_split%d.png", chunk.name, split.i)
      png.dir <- dirname(chunk.png)
      dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)
      print(chunk.png)
      png(chunk.png, width=9, h=7, res=200, units="in")
      print(chunkPlot)
      dev.off()
      
      problems.dt <- data.table(step2.problems)
      setkey(problems.dt, problemStart, problemEnd)
      over.regions <- foverlaps(regions, problems.dt, nomatch=0L)
      regions.by.problem <-
        split(over.regions, over.regions$problem.name, drop=TRUE)
      for(problem.name in names(regions.by.problem)){
        problem.regions <- regions.by.problem[[problem.name]]
        problem.i <- problem.regions$problem.i[1]
        problem <- problems.dt[problem.i, ]
        problem.counts <- foverlaps(counts, problem, nomatch=0L)
        fit <- PeakSegJointHeuristic(problem.counts)
        converted <- ConvertModelList(fit)
        prob.err.list <- PeakSegJointError(converted, problem.regions)
      }
    }#chunk.name
  }#split.i
}#set.name

save(step2, file="step2.RData")
