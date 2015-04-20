## -*- compile-command: "make HOCKING-PeakSegJoint-slides.pdf" -*-
works_with_R("3.2.0",
             xtable="1.7.4",
             data.table="1.9.5",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakSegDP@fe06a5b91d68c5d1ec471cb15c3ec3935dc2624d",
             "tdhock/PeakSegJoint@e6bc386c555e203cc80343814939d51785c03af1",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732")

load("~/projects/chip-seq-paper/chunks/H3K4me3_TDH_other/18/regions.RData")
load("~/projects/chip-seq-paper/chunks/H3K4me3_TDH_other/18/counts.RData")
counts$count <- counts$coverage

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

profile.list <- split(counts, counts$sample.id)
region.list <- split(regions, regions$sample.id)

chrom.range <- with(counts, {
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
  sample.regions <- region.list[[sample.id]]
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
    error <- PeakErrorChrom(peaks.df, sample.regions)
    sample.region.list[[peaks.str]] <-
      data.frame(sample.id, peaks, error)
    sample.error.list[[peaks.str]] <-
      data.frame(sample.id, peaks,
                 errors=with(error, sum(fp+fn)))
  }
  errors <- do.call(rbind, sample.error.list)
  best.error <- errors[which.min(errors$errors),]
  peaks.str <- paste(best.error$peaks)
  max.count <- max(profile$count)
  h <- max.count/10
  y <- -h
  error.region.list[[sample.id]] <-
    data.frame(ymin=y-h, ymax=y+h, sample.region.list[[peaks.str]])
  peak.list[[sample.id]] <-
    data.frame(sample.id, fit$fit$peaks[[peaks.str]], y)
  timing.list[[sample.id]] <- data.frame(seconds, sample.id)
}
error.regions <- do.call(rbind, error.region.list)
peaks.df <- do.call(rbind, peak.list)

positive.regions <- subset(regions, annotation %in% c("peakStart", "peakEnd"))
positive.bases <- with(positive.regions, chromEnd-chromStart)
positive.quartile <- quantile(positive.bases)

bases.per.problem.vec <-
  as.integer(sort(c((4.5) * (2^(9:14))
                    )))

count.dt <- data.table(counts)
setkey(count.dt, chromStart, chromEnd)

regions.dt <- data.table(regions)
regions.dt$region.i <- 1:nrow(regions.dt)
setkey(regions.dt, chromStart, chromEnd)

joint.time.list <- list()
best.peak.list <- list()
resolution.err.list <- list()
problem.list <- list()
for(bases.per.problem.i in seq_along(bases.per.problem.vec)){
  bases.per.problem <- bases.per.problem.vec[[bases.per.problem.i]]
  cat(sprintf("%4d / %4d %s\n", bases.per.problem.i,
              length(bases.per.problem.vec), bases.per.problem))
  problemSeq <- seq(0, chrom.range[2], by=bases.per.problem)
  problemStart <-
    as.integer(c(problemSeq,
                 problemSeq+bases.per.problem/3,
                 problemSeq+2*bases.per.problem/3))
  problemStart <-
    as.integer(c(problemSeq,
                 problemSeq+bases.per.problem/4,
                 problemSeq+2*bases.per.problem/4,
                 problemSeq+3*bases.per.problem/4))
  problemStart <-
    as.integer(c(problemSeq,
                 problemSeq+bases.per.problem/2))
  problemEnd <- problemStart+bases.per.problem
  peakStart <- as.integer(problemStart+bases.per.problem/4)
  peakEnd <- as.integer(problemEnd-bases.per.problem/4)
  is.overlap <- chrom.range[1] < problemEnd
  problem.name <- sprintf("chr21:%d-%d", problemStart, problemEnd)
  problems <- 
    data.frame(problem.name, peakStart, peakEnd,
               bases.per.problem, problemStart, problemEnd)[is.overlap,]
  all.problems.dt <- data.table(problems)
  setkey(all.problems.dt, problemStart, problemEnd)

  ggplot()+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  alpha=0.5,
                  color="grey",
                  data=regions)+
    scale_fill_manual(values=ann.colors)+
    geom_segment(aes(problemStart/1e3, problem.name,
                     xend=problemEnd/1e3, yend=problem.name),
                 data=problems)

  overlap.regions <- foverlaps(regions.dt, all.problems.dt, nomatch=0L)
  region.counts <- table(overlap.regions$region.i)
  region.weights <- 1/region.counts
  overlap.regions[, weight := region.weights[paste(region.i)]]
  table(paste(overlap.regions$problem.name))

  problems.with.regions <- paste(unique(overlap.regions$problem.name))
  setkey(overlap.regions, problem.name)
  setkey(all.problems.dt, problem.name)
  problems.dt <- all.problems.dt[problems.with.regions]
  problems.dt$problem.i <- 1:nrow(problems.dt)
  problem.list[[paste(bases.per.problem)]] <- problems.dt
  resolution.stat.list <- list()
  for(problem.i in seq_along(problems.with.regions)){
    problem.name <- problems.with.regions[[problem.i]]
    cat(sprintf("%4d / %4d problems %s\n", problem.i,
                length(problems.with.regions), problem.name))
    problem.dt <- problems.dt[problem.name]
    setkey(problem.dt, problemStart, problemEnd)
    problem.counts <- foverlaps(count.dt, problem.dt, nomatch=0L)
    problem.regions <- overlap.regions[problem.name]
    regions.by.sample <- split(problem.regions, problem.regions$sample.id)

    problem.count.list <- ProfileList(problem.counts)
    ## max(sapply(problem.count.list, with, min(chromStart)))
    ## min(sapply(problem.count.list, with, max(chromEnd)))
    seconds <- system.time({
      fit <- PeakSegJointHeuristic(problem.count.list)
      converted <- ConvertModelList(fit)
    })[["elapsed"]]

    joint.time.list[[problem.name]] <-
      data.frame(problem.dt, seconds,
                 data=nrow(problem.counts))

    segments.by.peaks <- if(!is.null(converted$segments$peaks)){
      with(converted, split(segments, segments$peaks))
    }
    peaks.by.peaks <- if(!is.null(converted$peaks$peaks)){
      with(converted, split(peaks, peaks$peaks))
    }else{
      list()
    }
    peaks.by.peaks[["0"]] <- Peaks()
    problem.err.list <- list()
    problem.stats.list <- list()
    for(peaks.str in names(peaks.by.peaks)){
      peaks.all.samples <- peaks.by.peaks[[peaks.str]]
      peaks.by.sample <-
        split(peaks.all.samples, peaks.all.samples$sample.id, drop=TRUE)
      peaks.stats.list <- list()
      for(sample.id in names(regions.by.sample)){
        sample.regions <- regions.by.sample[[sample.id]]
        sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
          peaks.by.sample[[sample.id]]
        }else{
          Peaks()
        }
        error <- PeakErrorChrom(sample.peaks, sample.regions)
        error$weight <- sample.regions$weight
        error$errors <- with(error, fp+fn)
        error$weighted.error <- with(error, errors * weight)
        problem.err.list[[peaks.str]][[sample.id]] <-
          data.frame(error, sample.id)
        peaks.stats.list[[sample.id]] <-
          data.frame(sample.id, weighted.error=sum(error$weighted.error),
                     total.weight=sum(error$weight))
      }
      peaks.stats <- do.call(rbind, peaks.stats.list)
      problem.stats.list[[peaks.str]] <-
        data.frame(peaks=as.numeric(peaks.str),
                   weighted.error=sum(peaks.stats$weighted.error),
                   total.weight=sum(peaks.stats$total.weight))
    }
    problem.stats <- do.call(rbind, problem.stats.list)
    problem.stats <- problem.stats[order(problem.stats$peaks), ]
    best.error <- problem.stats[which.min(problem.stats$weighted.error),]
    if(best.error$weighted.error > 0){
      ##stop("non-zero error")
    }
    resolution.stat.list[[problem.name]] <- best.error

    ggplot(,aes(peaks, weighted.error))+
      geom_line(data=problem.stats)+
      geom_point(data=best.error)
    
    peaks.str <- paste(best.error$peaks)
    best.regions <- do.call(rbind, problem.err.list[[peaks.str]])
    best.peaks <- peaks.by.peaks[[peaks.str]]
    best.segments <- segments.by.peaks[[peaks.str]]

    if(nrow(best.peaks)){
      best.peak.list[[paste(bases.per.problem)]][[problem.name]] <-
        data.frame(best.peaks, problem.i=problem.dt$problem.i)
    }

    count.na.list <- list()
    for(sample.id in names(problem.count.list)){
      df <- problem.count.list[[sample.id]]
      r <- tail(df, 1)
      r$chromStart <- r$chromEnd
      count.na.list[[sample.id]] <-
        data.frame(sample.id, rbind(df, r))
    }
    count.na <- do.call(rbind, count.na.list)

    resPlot <- 
    ggplot()+
      scale_linetype_manual("error type",
                            limits=c("correct", 
                              "false negative",
                              "false positive"
                                     ),
                            values=c(correct=0,
                              "false negative"=3,
                              "false positive"=1))+
      geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                        fill=annotation),
                    alpha=0.5,
                    color="grey",
                    data=best.regions)+
      geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                        linetype=status),
                    fill=NA,
                    size=2,
                    color="black",
                    data=best.regions)+
      scale_y_continuous("aligned read coverage",
                         breaks=function(limits){
                           floor(limits[2])
                         })+
      xlab("position on chromosome (kilo bases = kb)")+
      geom_step(aes(chromStart/1e3, count),
                data=count.na,
                size=1,
                color="grey50")+
      ## geom_segment(aes(chromStart/1e3, 0,
      ##                  xend=chromEnd/1e3, yend=0),
      ##              size=3,
      ##              color="deepskyblue",
      ##              data=best.peaks)+
      geom_segment(aes(chromStart/1e3, mean,
                       xend=chromEnd/1e3, yend=mean),
                   data=best.segments,
                   color="green")+
      theme_bw()+
      scale_fill_manual(values=ann.colors)+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
        sub("McGill0", "", val)
      })
  }#problem.name
  resolution.stats <- do.call(rbind, resolution.stat.list)
  stopifnot(all.equal(sum(resolution.stats$total.weight), nrow(regions.dt)))
  resolution.err.list[[paste(bases.per.problem)]] <-
    data.frame(bases.per.problem,
               weighted.error=sum(resolution.stats$weighted.error))
}
resolution.err <- do.call(rbind, resolution.err.list)
print(resolution.err)
best.res <- resolution.err[which.min(resolution.err$weighted.error), ]
best.res.str <- paste(best.res$bases.per.problem)
best.peaks <- do.call(rbind, best.peak.list[[best.res.str]])
best.problems <- problem.list[[best.res.str]]
prob.peaks <- best.peaks
prob.peaks$sample.id <- "problems"

all.problems <- do.call(rbind, problem.list)
ggplot()+
  geom_point(aes(problemStart/1e3, bases.per.problem),
               data=all.problems)+
  geom_segment(aes(problemStart/1e3, bases.per.problem,
                   xend=problemEnd/1e3, yend=bases.per.problem),
               data=all.problems)+
  scale_y_log10()

problemsPlot <- 
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
                data=regions)+
  scale_color_manual(values=c(PeakSeg="deepskyblue",
                       PeakSegJoint="black"),
                     limits=c("PeakSegJoint", "PeakSeg"))+
  geom_segment(aes(chromStart/1e3, y,
                   color=model,
                   xend=chromEnd/1e3, yend=y),
               size=2,
               data=data.frame(peaks.df, model="PeakSeg"))+
  geom_segment(aes(chromStart/1e3, 0,
                   color=model,
                   xend=chromEnd/1e3, yend=0),
               size=2,
               data=data.frame(best.peaks, model="PeakSegJoint"))+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=counts,
            color="grey50")+
  theme_bw()+
  geom_segment(aes(chromStart/1e3, problem.i,
                   color=model,
                   xend=chromEnd/1e3, yend=problem.i),
               size=2,
               data=data.frame(prob.peaks, model="PeakSegJoint"))+
  geom_segment(aes(problemStart/1e3, problem.i,
                   xend=problemEnd/1e3, yend=problem.i),
               data=data.frame(sample.id="problems", best.problems))+
  geom_text(aes(max(best.problems$problemStart)/1e3, 1,
                label=paste(bases.per.problem, "bases/problem")),
            vjust=0, 
            data=data.frame(sample.id="problems", best.res))+
  theme(panel.margin=grid::unit(0, "cm"))+
  coord_cartesian(xlim=(chrom.range+c(-1,1)*chrom.bases*0/8)/1e3)+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })

png("figure-H3K4me3-profiles.png", width=9, h=7, res=200, units="in")
print(problemsPlot)
dev.off()

joint.times <- do.call(rbind, joint.time.list)
joint.by.res <- split(joint.times, joint.times$bases.per.problem)
best.times <- joint.by.res[[best.res.str]]
best.times.df <-
  rbind(best.times[, c("data", "seconds")],
        data.frame(data=NA, seconds=sum(best.times$seconds), row.names="total"))

single.times <- do.call(rbind, timing.list)
ggplot()+
  geom_point(aes(data, seconds), data=joint.times)
resolution.times <- joint.times %>%
  group_by(bases.per.problem) %>%
  summarise(seconds=sum(seconds))

single.times.df <-
  rbind(single.times,
        data.frame(seconds=sum(single.times$seconds), sample.id="total"))
single.times.xt <- xtable(single.times.df)
print(single.times.xt, include.rownames=FALSE, floating=FALSE,
      file="table-H3K4me3-PeakSeg.tex")

xt <- xtable(best.times.df)
print(xt, file="table-H3K4me3.tex", floating=FALSE)
