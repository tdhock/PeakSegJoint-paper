works_with_R("3.2.0",
             data.table="1.9.4",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@dce1572ae2067c0dd6b0e6c6efef99e0de2e3274")

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

for(set.name in dir("PeakSegJoint-chunks")){
  set.dir <- file.path("PeakSegJoint-chunks", set.name)
  problems.RData.vec <- Sys.glob(file.path(set.dir, "*", "problems.RData"))
  for(problems.RData in problems.RData.vec){
    objs <- load(problems.RData)
    if(min(step2.error$errors) > 0){
      regions.RData <- sub("problems", "regions", problems.RData)
      robjs <- load(regions.RData)
      counts.list <- list()
      chunk.dir <- dirname(problems.RData)
      counts.RData.vec <- Sys.glob(file.path(chunk.dir, "*", "*.RData"))
      for(counts.RData in counts.RData.vec){
        cell.type <- basename(dirname(counts.RData))
        sample.id <- sub(".RData$", "", basename(counts.RData))
        load(counts.RData)
        counts.list[[counts.RData]] <-
          data.table(cell.type, sample.id, counts)
      }
      counts <- do.call(rbind, counts.list)
      best.res <- with(step2.error, {
        paste(bases.per.problem[which.min(errors)])
      })
      step2.by.problem <- step2.by.res[[best.res]]
      step2.peak.list <- list()
      sample.peak.list <- list()
      pos.mat <- str_match_perl(names(step2.by.problem), genome.pos.pattern)
      for(problem.i in seq_along(step2.by.problem)){
        pos.row <- pos.mat[problem.i, ]
        problem.name <- pos.row[[1]]
        problem <- step2.by.problem[[problem.i]]
        min.err <- subset(problem$error$modelSelection, errors==min(errors))
        best.peaks <- min(min.err$peaks)
        if(is.data.frame(problem$converted$peaks)){
          selected.peaks <- subset(problem$converted$peaks, peaks==best.peaks)
          if(nrow(selected.peaks) > 0){
            sample.peak.list[[problem.name]] <- selected.peaks
          }
        }
        problemStart <- as.integer(pos.row[["chromStart"]])
        problemEnd <- as.integer(pos.row[["chromEnd"]])
        step2.peaks <-
          data.frame(problem.i,
                     problemStart,
                     problemEnd,
                     selected.peaks[1,])
        step2.peaks$sample.id <- "step 2"
        step2.peak.list[[problem.i]] <- step2.peaks
      }
      step2.peaks <- do.call(rbind, step2.peak.list)
      sample.peaks <- if(length(sample.peak.list) == 0){
        Peaks()
      }else{
        do.call(rbind, sample.peak.list)
      }

      error.regions <- PeakErrorSamples(sample.peaks, regions)

      limits <- regions[, .(chromStart=min(chromStart),
                            chromEnd=max(chromEnd))]
      limits[, bases := chromEnd - chromStart ]
      limits[, expand := as.integer(bases/10) ]
      limits[, min := chromStart - expand]
      limits[, max := chromEnd + expand]
      lim.vec <- with(limits, c(min, max))/1e3
      some.counts <- counts[limits$min < chromEnd &
                              chromStart < limits$max,]

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
          scale_x_continuous("position on chromosome (kilo bases = kb)")+
          coord_cartesian(xlim=lim.vec)+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                            fill=annotation),
                        alpha=0.5,
                        color="grey",
                        data=regions)+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                            linetype=status),
                        color="black",
                        fill=NA,
                        size=1.5,
                        data=error.regions)+
          scale_fill_manual(values=ann.colors)+
          theme_bw()+
          theme(panel.margin=grid::unit(0, "cm"))+
          facet_grid(sample.id ~ ., labeller=function(var, val){
            sub("McGill0", "", sub(" ", "\n", val))
          }, scales="free")+
          geom_step(aes(chromStart/1e3, count),
                    data=some.counts,
                    color="grey50")+
          geom_segment(aes(problemStart/1e3, problem.i,
                           xend=problemEnd/1e3, yend=problem.i),
                       data=step2.peaks)+
          geom_segment(aes(chromStart/1e3, problem.i,
                           xend=chromEnd/1e3, yend=problem.i),
                       size=2,
                       color="deepskyblue",
                       data=step2.peaks)
      if(is.data.frame(sample.peaks) && nrow(sample.peaks) > 0){
        problemPlot <- problemPlot+
          geom_segment(aes(chromStart/1e3, 0,
                           xend=chromEnd/1e3, yend=0),
                       size=2,
                       color="deepskyblue",
                       data=sample.peaks)
      }

      png.base <- gsub("/", "_", sub(".*?/", "", chunk.dir))
      png.name <-
        sprintf("figure-train-errors/%s.png",
                png.base)
      png.dir <- dirname(png.name)
      dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)
      print(png.name)
      png(png.name, width=15, h=10, res=100, units="in")
      print(problemPlot)
      dev.off()
    }
  }#chunk.name
}#set.name

