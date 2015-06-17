## -*- compile-command: "make HOCKING-PeakSegJoint-slides.pdf" -*-
works_with_R("3.2.0",
             xtable="1.7.4",
             data.table="1.9.5",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakSegDP@6520faad59fcf8c94ed4345d2b0d826b7f61faf9",
             "tdhock/PeakSegJoint@0390d43851e5427d3b2b45cb5a2b80b8007db67a",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732")

data(H3K36me3.TDH.other.chunk1)

all.regions <- H3K36me3.TDH.other.chunk1$regions
to.rep <- all.regions$chromEnd==43223588
all.regions$chromStart[to.rep] <- 43205000
all.regions$chromEnd[to.rep] <- 43260000

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

profile.list <- with(H3K36me3.TDH.other.chunk1, split(counts, counts$sample.id))
region.list <- split(all.regions, all.regions$sample.id)
## Comment for() loop below to ignore extra peakEnd annotation.
## for(sample.id in names(region.list)){
##   df <- region.list[[sample.id]]
##   r <- df[1,]
##   ##r$chromStart <- 43320000
##   r$chromStart <- 43270000
##   r$chromEnd <- 43380000
##   r$annotation <- "peakEnd"
##   region.list[[sample.id]] <- rbind(df, r)
## }

chrom.range <- with(H3K36me3.TDH.other.chunk1$counts, {
  c(min(chromStart), max(chromEnd))
})
chrom.bases <- chrom.range[2]-chrom.range[1]
n.bins <- 2000L
bases.per.bin <- as.integer(chrom.bases/n.bins)
load("../PeakSeg-paper/dp.timings.RData")
chunk.timings <-
  subset(dp.timings, set.name=="H3K36me3_TDH_other" & chunk.id==1,
         select=c(sample.id, bases, data, seconds))
rownames(chunk.timings) <- chunk.timings$sample.id
chunk.timings$seconds2000 <- NA
fit.list <- list()
error.region.list <- list()
peak.list <- list()
for(sample.id in names(profile.list)){
  profile <- profile.list[[sample.id]]
  regions <- region.list[[sample.id]]
  chunk.timings[sample.id, "seconds2000"] <- system.time({
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
  max.count <- max(profile$count)
  h <- max.count/10
  y <- -h
  error.region.list[[sample.id]] <-
    data.frame(ymin=y-h, ymax=y+h, sample.region.list[[peaks.str]])
  peak.list[[sample.id]] <-
    data.frame(sample.id, fit$fit$peaks[[peaks.str]], y)
}
error.regions <- do.call(rbind, error.region.list)
peaks.df <- do.call(rbind, peak.list)

positive.regions <- subset(regions, annotation %in% c("peakStart", "peakEnd"))
positive.bases <- with(positive.regions, chromEnd-chromStart)
positive.quartile <- quantile(positive.bases)

bases.per.problem.vec <-
  as.integer(sort(c((4.5) * (2^(12:17))
                    ##,(4.9) * (2^(12:16))
                    )))
bases.per.problem.vec <- 200000

count.dt <- data.table(H3K36me3.TDH.other.chunk1$counts)
setkey(count.dt, chromStart, chromEnd)

regions.dt <- data.table(all.regions)
regions.dt$region.i <- 1:nrow(regions.dt)
setkey(regions.dt, chromStart, chromEnd)

joint.time.list <- list()
best.peak.list <- list()
resolution.err.list <- list()
problem.list <- list()
for(bases.per.problem in bases.per.problem.vec){
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
  problems.dt <- data.table(problems)
  setkey(problems.dt, problemStart, problemEnd)
  problems.dt$problem.i <- 1:nrow(problems.dt)

  ggplot()+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  alpha=0.5,
                  color="grey",
                  data=all.regions)+
    scale_fill_manual(values=ann.colors)+
    geom_segment(aes(problemStart/1e3, problem.name,
                     xend=problemEnd/1e3, yend=problem.name),
                 data=problems)

  overlap.regions <- foverlaps(regions.dt, problems.dt, nomatch=0L)
  region.counts <- table(overlap.regions$region.i)
  region.weights <- 1/region.counts
  overlap.regions[, weight := region.weights[paste(region.i)]]
  table(paste(overlap.regions$problem.name))

  problems.with.regions <- paste(unique(overlap.regions$problem.name))
  setkey(overlap.regions, problem.name)
  setkey(problems.dt, problem.name)
  problem.list[[paste(bases.per.problem)]] <-
    problems.dt[problems.with.regions]
  resolution.stat.list <- list()
  for(problem.name in problems.with.regions){
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
        regions <- regions.by.sample[[sample.id]]
        peaks <- if(sample.id %in% names(peaks.by.sample)){
          peaks.by.sample[[sample.id]]
        }else{
          Peaks()
        }
        error <- PeakErrorChrom(peaks, regions)
        error$weight <- regions$weight
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
    vline.df <-
      data.frame(base=c(fit$data_start_end[2], max(best.segments$chromStart)),
                 what=c("data end", "peak end"))
    ## this next plot shows the right edge of the data, and that in
    ## particular the peak does not end at the exact same location as
    ## the data.
    
    ## bases.per.problem <- 294912
    ##problem.name <- "chr21:43057152-43352064"
    resPlot+coord_cartesian(xlim=c(43351900, 43352100)/1e3)+
      geom_vline(aes(xintercept=base/1e3), data=vline.df)+
      geom_text(aes(base/1e3, 0, label=paste0(what, " ")),
                hjust=1,
                data=data.frame(sample.id="label", vline.df))
  }#problem.name
  resolution.stats <- do.call(rbind, resolution.stat.list)
  stopifnot(all.equal(sum(resolution.stats$total.weight), nrow(regions.dt)))
  resolution.err.list[[paste(bases.per.problem)]] <-
    data.frame(bases.per.problem,
               weighted.error=sum(resolution.stats$weighted.error))
}
resolution.err <- do.call(rbind, resolution.err.list)
print(resolution.err)
best.res <- resolution.err[1, ]
best.res.str <- paste(best.res$bases.per.problem)
best.peaks <- do.call(rbind, best.peak.list[[best.res.str]])
best.problems <- problem.list[[best.res.str]]
prob.peaks <-
  data.frame(unique(best.peaks[, c("problem.i", "chromStart", "chromEnd")]),
             sample.id="step 1")

all.problems <- do.call(rbind, problem.list)
ggplot()+
  geom_point(aes(problemStart/1e3, bases.per.problem),
               data=all.problems)+
  geom_segment(aes(problemStart/1e3, bases.per.problem,
                   xend=problemEnd/1e3, yend=bases.per.problem),
               data=all.problems)+
  scale_y_log10()

uniq.peaks <- unique(best.peaks[, c("chromStart", "chromEnd")])
clustered.peaks <- clusterPeaks(uniq.peaks)
cluster.list <- split(clustered.peaks, clustered.peaks$cluster)
newprob.list <- list()
refined.list <- list()
for(cluster.name in names(cluster.list)){
  cluster <- cluster.list[[cluster.name]]
  cluster.chromStart <- min(cluster$chromStart)
  cluster.chromEnd <- max(cluster$chromEnd)
  cluster.mid <- as.integer((cluster.chromEnd + cluster.chromStart)/2)
  half.bases <- as.integer(best.res$bases.per.problem/2)
  cluster.num <- as.numeric(cluster.name)
  before.name <- paste(cluster.num-1)
  chromEnd.before <- if(before.name %in% names(cluster.list)){
    max(cluster.list[[before.name]]$chromEnd)
  }else{
    0
  }
  after.name <- paste(cluster.num+1)
  chromStart.after <- if(after.name %in% names(cluster.list)){
    min(cluster.list[[after.name]]$chromStart)
  }else{
    Inf
  }
  problemStart <- cluster.mid - half.bases
  if(problemStart < chromEnd.before){
    problemStart <- as.integer((chromEnd.before+cluster.chromStart)/2)
  }
  problemEnd <- cluster.mid + half.bases
  if(chromStart.after < problemEnd){
    problemEnd <- as.integer((chromStart.after+cluster.chromEnd)/2)
  }
  stopifnot(problemStart < cluster.chromStart)
  stopifnot(cluster.chromEnd < problemEnd)
  
  problem.i <- cluster.num+1
  newprob.list[[cluster.name]] <- newprob <- 
    data.table(problem.i, cluster.name, cluster.num,
               problemStart,
               problemEnd)

  setkey(newprob, problemStart, problemEnd)
  problem.counts <- foverlaps(count.dt, newprob, nomatch=0L)
  problem.regions <- foverlaps(regions.dt, newprob, nomatch=0L)
  regions.by.sample <- split(problem.regions, problem.regions$sample.id)
  problem.count.list <- ProfileList(problem.counts)
  seconds <- system.time({
    fit <- PeakSegJointHeuristic(problem.count.list)
    converted <- ConvertModelList(fit)
  })[["elapsed"]]

  peaks.by.peaks <- if(!is.null(converted$peaks$peaks)){
    with(converted, split(peaks, peaks$peaks))
  }else{
    list()
  }
  peaks.by.peaks[["0"]] <- Peaks()
  problem.stats.list <- list()
  for(peaks.str in names(peaks.by.peaks)){
    peaks.all.samples <- peaks.by.peaks[[peaks.str]]
    peaks.by.sample <-
      split(peaks.all.samples, peaks.all.samples$sample.id, drop=TRUE)
    peaks.stats.list <- list()
    for(sample.id in names(regions.by.sample)){
      regions <- regions.by.sample[[sample.id]]
      peaks <- if(sample.id %in% names(peaks.by.sample)){
        peaks.by.sample[[sample.id]]
      }else{
        Peaks()
      }
      error <- PeakErrorChrom(peaks, regions)
      error$weight <- regions$weight
      error$errors <- with(error, fp+fn)
      peaks.stats.list[[sample.id]] <-
        data.frame(sample.id, errors=sum(error$errors),
                   regions=nrow(error))
    }
    peaks.stats <- do.call(rbind, peaks.stats.list)
    problem.stats.list[[peaks.str]] <-
      data.frame(peaks=as.numeric(peaks.str),
                 errors=sum(peaks.stats$errors),
                 regions=sum(peaks.stats$regions))
  }
  problem.stats <- do.call(rbind, problem.stats.list)
  sorted.stats <- problem.stats[order(problem.stats$peaks),]
  problem.param <- sorted.stats[which.min(sorted.stats$errors), "peaks"]
  refined.list[[cluster.name]] <-
    data.frame(problem.i, peaks.by.peaks[[problem.param]])
}
newprobs <- do.call(rbind, newprob.list)
refined <- do.call(rbind, refined.list)
refined.below <-
  data.frame(unique(refined[, c("problem.i", "chromStart", "chromEnd")]),
             sample.id="step 2")

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
                data=all.regions)+
  scale_color_manual(values=c(PeakSeg="black",
                       PeakSegJoint="deepskyblue"),
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
               data=data.frame(refined, model="PeakSegJoint"))+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=H3K36me3.TDH.other.chunk1$counts,
            color="grey50")+
  theme_bw()+
  geom_segment(aes(chromStart/1e3, problem.i,
                   color=model,
                   xend=chromEnd/1e3, yend=problem.i),
               size=2,
               data=data.frame(prob.peaks, model="PeakSegJoint"))+
  geom_segment(aes(chromStart/1e3, problem.i,
                   color=model,
                   xend=chromEnd/1e3, yend=problem.i),
               size=2,
               data=data.frame(refined.below, model="PeakSegJoint"))+
  geom_segment(aes(problemStart/1e3, problem.i,
                   xend=problemEnd/1e3, yend=problem.i),
               data=data.frame(sample.id="step 1", best.problems))+
  geom_segment(aes(problemStart/1e3, problem.i,
                   xend=problemEnd/1e3, yend=problem.i),
               data=data.frame(sample.id="step 2", newprobs))+
  geom_text(aes(max(best.problems$problemStart)/1e3, 1,
                label=sprintf("%d bases/problem", bases.per.problem)),
            vjust=0, 
            data=data.frame(sample.id="step 1", best.res))+
  theme(panel.margin=grid::unit(0, "cm"))+
  coord_cartesian(xlim=(chrom.range+c(-1,1)*chrom.bases/3)/1e3)+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })

problemsPlot+
  coord_cartesian(xlim=c(43214142, 43223588)/1e3)

png("figure-H3K36me3-profiles.png", width=9, h=6, res=200, units="in")
print(problemsPlot)
dev.off()

joint.times <- do.call(rbind, joint.time.list)
joint.by.res <- split(joint.times, joint.times$bases.per.problem)
best.times <- joint.by.res[[best.res.str]]
best.times.df <-
  rbind(best.times[, c("data", "seconds")],
        data.frame(data=NA, seconds=sum(best.times$seconds), row.names="total"))

single.times <- chunk.timings
ggplot()+
  geom_point(aes(data, seconds), data=joint.times)
resolution.times <- joint.times %>%
  group_by(bases.per.problem) %>%
  summarise(seconds=sum(seconds))

single.times.df <-
  rbind(single.times,
        data.frame(sample.id="total",
                   bases="",
                   data="",
                   seconds=sum(single.times$seconds),
                   seconds2000=sum(single.times$seconds2000)))
single.times.xt <- xtable(single.times.df)
print(single.times.xt, include.rownames=FALSE, floating=FALSE,
      file="table-H3K36me3-PeakSeg.tex")

xt <- xtable(best.times.df)
print(xt, file="table-H3K36me3.tex", floating=FALSE)
