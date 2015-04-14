works_with_R("3.1.3",
             data.table="1.9.5",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakSegDP@fe06a5b91d68c5d1ec471cb15c3ec3935dc2624d",
             "tdhock/PeakSegJoint@547ca81ce92c38b7c7f78cd071efc5afa96cb289",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732")

data(H3K36me3.TDH.other.chunk1)

profile.list <- with(H3K36me3.TDH.other.chunk1, split(counts, counts$sample.id))
region.list <-
  with(H3K36me3.TDH.other.chunk1, split(regions, regions$sample.id))
chrom.range <- with(H3K36me3.TDH.other.chunk1$counts, {
  c(min(chromStart), max(chromEnd))
})
chrom.bases <- chrom.range[2]-chrom.range[1]
n.bins <- 2000L
bases.per.bin <- as.integer(chrom.bases/n.bins)
fit.list <- list()
timing.list <- list()
error.region.list <- list()
peak.list <- list()
for(sample.id in names(profile.list)){
  profile <- profile.list[[sample.id]]
  regions <- region.list[[sample.id]]
  seconds <- system.time({
    fit.list[[sample.id]] <- fit <- 
      segmentBins(profile,
                  chrom.range[1], bases.per.bin)
  })[["elapsed"]]
  sample.error.list <- list()
  sample.region.list <- list()
  for(peaks.str in names(fit$fit$peaks)){
    peaks.df <- fit$fit$peaks[[peaks.str]]
    peaks <- as.numeric(peaks.str)
    error <- PeakErrorChrom(peaks.df, regions)
    sample.region.list[[peaks.str]] <-
      data.frame(sample.id, peaks, error)
    sample.error.list[[peaks.str]] <-
      data.frame(sample.id, peaks,
                 errors=with(error, sum(fp+fn)))
  }
  errors <- do.call(rbind, sample.error.list)
  best.error <- errors[which.min(errors$errors),]
  peaks.str <- paste(best.error$peaks)
  y <- max(profile$count)
  h <- y/10
  error.region.list[[sample.id]] <-
    data.frame(ymin=y-h, ymax=y+h, sample.region.list[[peaks.str]])
  peak.list[[sample.id]] <-
    data.frame(sample.id, fit$fit$peaks[[peaks.str]], y)
  timing.list[[sample.id]] <- data.frame(seconds, sample.id)
}
error.regions <- do.call(rbind, error.region.list)
peaks.df <- do.call(rbind, peak.list)

bases.per.problem.vec <- (4.5) * (2^(12:17))

problem.list <- list()
for(bases.per.problem in bases.per.problem.vec){
  problemStart <- seq(0, chrom.range[2], by=bases.per.problem)
  problemStart <- c(problemStart, problemStart+bases.per.problem/2)
  problemEnd <- problemStart+bases.per.problem
  is.overlap <- chrom.range[1] < problemEnd
  problem.name <- sprintf("chr21:%d-%d", problemStart, problemEnd)
  ## TODO: weight regions based on number of overlapping problems,
  ## compute weighted error.
  problem.list[[paste(bases.per.problem)]] <- 
    data.frame(problem.name,
               bases.per.problem, problemStart, problemEnd)[is.overlap,]
}
problems <- do.call(rbind, problem.list)
ggplot()+
  geom_point(aes(problemStart/1e3, bases.per.problem),
               data=problems)+
  geom_segment(aes(problemStart/1e3, bases.per.problem,
                   xend=problemEnd/1e3, yend=bases.per.problem),
               data=problems)+
  scale_y_log10()

problem.dt <- data.table(problems)
setkey(problem.dt, problemStart, problemEnd)
count.dt <- data.table(H3K36me3.TDH.other.chunk1$counts)
setkey(count.dt, chromStart, chromEnd)
overlap.dt <- foverlaps(count.dt, problem.dt)

joint.time.list <- list()
overlap.list <- split(overlap.dt, overlap.dt$problem.name, drop=TRUE)
for(problem.name in names(overlap.list)){
  problem.coverage <- overlap.list[[problem.name]]
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    xlab("position on chromosome (kilo bases = kb)")+
    geom_line(aes(chromStart/1e3, count),
              data=problem.coverage,
              color="grey50")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  table(problem.coverage$sample.id)
  seconds <- system.time({
    fit <- PeakSegJointHeuristic(problem.coverage)
    converted <- ConvertModelList(fit)
  })[["elapsed"]]
  joint.time.list[[problem.name]] <-
    data.frame(problem.name, seconds,
               bases.per.problem=problem.coverage$bases.per.problem[1],
               data=nrow(problem.coverage))
  peaks.by.sample <- with(converted, split(peaks, peaks$sample.id))
  for(sample.id in names(peaks.by.sample)){
    sample.peaks <- peaks.by.sample[[sample.id]]
    sample.regions <- region.list[[sample.id]]
    peaks.by.peaks <- split(sample.peaks, sample.peaks$peaks)
  }
}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

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
  geom_blank(aes(linetype=status),
             data=data.frame(status=c("correct", "false positive",
                               "false negative")))+
  xlab("position on chromosome (kilo bases = kb)")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                alpha=0.5,
                color="grey",
                data=H3K36me3.TDH.other.chunk1$regions)+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                ymin=ymin, ymax=ymax,
                linetype=status),
            fill=NA,
            color="black",
            data=error.regions)+
  scale_color_manual(values=c(PeakSeg="deepskyblue",
                       PeakSegJoint="black"))+
  geom_segment(aes(chromStart/1e3, y,
                   color=model,
                   xend=chromEnd/1e3, yend=y),
               size=3,
               data=data.frame(peaks.df, model="PeakSeg"))+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=H3K36me3.TDH.other.chunk1$counts,
            color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })
