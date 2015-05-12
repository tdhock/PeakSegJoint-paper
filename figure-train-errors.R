works_with_R("3.2.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@f1682b7e9a3d11c410c208325610ea1ede13edfa")

load("selected.by.set.RData")

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

regions.file.list <- list()
regions.file.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/regions.RData")
for(regions.file in regions.file.vec){
  load(regions.file)
  regions.file.list[[regions.file]] <- regions
}

RData.vec <- Sys.glob("chunk.problems/*/*/*.RData")
chunk.list <- list()
for(RData in RData.vec){
  objs <- load(RData)
  chunk.list[[RData]] <- best.step2.error
}

genome.pos.pattern <-
  paste0("(?<chrom>chr.*?)",
         ":",
         "(?<chromStart>[0-9]+)",
         "-",
         "(?<chromEnd>[0-9]+)")

for(set.name in names(selected.by.set)){
  selected.df <- selected.by.set[[set.name]]
  chunk.ids <- dir(file.path("chunk.problems", set.name))
  set.chunks <- paste0(set.name, "/", chunk.ids)
  for(chunk.name in set.chunks){
    RData.vec <- Sys.glob(sprintf("chunk.problems/%s/*.RData", chunk.name))
    error.list <- list()
    for(RData in RData.vec){
      objs <- load(RData)
      error.list[[RData]] <- best.step2.error
    }
    error.df <- do.call(rbind, error.list)
    error.ordered <- error.df[order(error.df$bases.per.problem),]
    if(min(error.ordered$errors) > 0){
      best.res <- with(error.ordered, {
        paste(bases.per.problem[which.min(errors)])
      })
      RData <- grep(best.res, RData.vec, value=TRUE)
      load(RData)
      step1.peak.list <- list()
      pos.mat <- str_match_perl(names(step1.by.problem), genome.pos.pattern)
      for(problem.i in seq_along(step1.by.problem)){
        problem <- step1.by.problem[[problem.i]]
        step1.peaks.all <- problem$peaks[[length(problem$peaks)]]
        pos.row <- pos.mat[problem.i, ]
        problemStart <- as.integer(pos.row[["chromStart"]])
        problemEnd <- as.integer(pos.row[["chromEnd"]])
        step1.row <- step1.peaks.all[1,]
        step1.peaks <-
          data.frame(problem.i,
                     problemStart,
                     problemEnd,
                     chromStart=step1.row$chromStart,
                     chromEnd=step1.row$chromEnd)
        step1.peaks$sample.id <- "step 1"
        step1.peak.list[[problem.i]] <- step1.peaks
      }
      step1.peaks <- do.call(rbind, step1.peak.list)

      step2.peak.list <- list()
      sample.peak.list <- list()
      pos.mat <- str_match_perl(names(step2.by.problem), genome.pos.pattern)
      for(problem.i in seq_along(step2.by.problem)){
        pos.row <- pos.mat[problem.i, ]
        problem.name <- pos.row[[1]]
        problem <- step2.by.problem[[problem.i]]
        min.err <- subset(problem$error$error.totals, errors==min(errors))
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
      sample.peaks <- do.call(rbind, sample.peak.list)

      chunk.dir <- paste0("../chip-seq-paper/chunks/", chunk.name)
      load(file.path(chunk.dir, "regions.RData"))
      load(file.path(chunk.dir, "counts.RData"))

      error.regions <- PeakErrorSamples(sample.peaks, regions)

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
          geom_step(aes(chromStart/1e3, coverage),
                    data=counts,
                    color="grey50")+
          geom_segment(aes(chromStart/1e3, 0,
                           xend=chromEnd/1e3, yend=0),
                       size=2,
                       color="deepskyblue",
                       data=sample.peaks)+
          geom_segment(aes(problemStart/1e3, problem.i,
                           xend=problemEnd/1e3, yend=problem.i),
                       color="deepskyblue",
                       data=step1.peaks)+
          geom_segment(aes(chromStart/1e3, problem.i,
                           xend=chromEnd/1e3, yend=problem.i),
                       color="deepskyblue",
                       size=2,
                       data=step1.peaks)+
          geom_segment(aes(problemStart/1e3, problem.i,
                           xend=problemEnd/1e3, yend=problem.i),
                       color="deepskyblue",
                       data=step2.peaks)+
          geom_segment(aes(chromStart/1e3, problem.i,
                           xend=chromEnd/1e3, yend=problem.i),
                       size=2,
                       color="deepskyblue",
                       data=step2.peaks)
      png.name <- sprintf("figure-train-errors/%s.png", chunk.name)
      png.dir <- dirname(png.name)
      dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)
      if(!file.exists(png.name)){
        print(png.name)
        png(png.name, width=9, h=7, res=200, units="in")
        print(problemPlot)
        dev.off()
      }
    }
  }#chunk.name
}#set.name

