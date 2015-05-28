works_with_R("3.2.0",
             ggplot2="1.0",
             "tdhock/animint@3ee3b962ba5737d8b291813fe6e73a75705f5cdf",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@d99652f043f738c25fad22737999eb794c264c54")

##load("TF.benchmark.corrected.RData")

set.seed(4)

doSim <- function(sample.id){
  data.vec <- do.call(c, lapply(c(8, 4, 15, 12, 5, 7), function(mean.value){
    rpois(4, mean.value)
  }))
  data.chromEnd <- seq_along(data.vec)
  data.frame(sample.id,
             chromStart=data.chromEnd-1L,
             chromEnd=data.chromEnd,
             count=data.vec)
}
sim <- do.call(rbind, lapply(c("sample1", "sample2"), doSim))
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ .)+
  geom_point(aes(chromEnd, count), data=sim)

## Begin R implementation of multiple sample constrained
## segmentation heuristic. Input: profiles data.frame.
profiles <- sim
unfilled.profile.list <- split(profiles, profiles$sample.id, drop=TRUE)
unfilled.chromStart <- max(sapply(unfilled.profile.list, with, chromStart[1]))
unfilled.chromEnd <-
  min(sapply(unfilled.profile.list, with, chromEnd[length(chromEnd)]))
unfilled.bases <- unfilled.chromEnd-unfilled.chromStart
bin.factor <- 2L
bases.per.bin <- 1L
while(unfilled.bases/bases.per.bin/bin.factor >= 4){
  bases.per.bin <- bases.per.bin * bin.factor
}
print(bases.per.bin)
n.bins <- as.integer(unfilled.bases %/% bases.per.bin + 1L)

extra.bases <- n.bins * bases.per.bin - unfilled.bases
extra.before <- as.integer(extra.bases/2)
extra.after <- extra.bases - extra.before
max.chromStart <- unfilled.chromStart-extra.before
min.chromEnd <- unfilled.chromEnd + extra.after
profile.list <- list()
for(sample.id in names(unfilled.profile.list)){
  one.sample <- 
    subset(unfilled.profile.list[[sample.id]],
           unfilled.chromStart < chromEnd &
             chromStart <= unfilled.chromEnd)
  one.sample$chromStart[1] <- unfilled.chromStart
  one.sample$chromEnd[nrow(one.sample)] <- unfilled.chromEnd
  stopifnot(with(one.sample, sum(chromEnd-chromStart)) == unfilled.bases)
  first.row <- last.row <- one.sample[1,]
  first.row$chromStart <- max.chromStart
  first.row$chromEnd <- unfilled.chromStart
  first.row$count <- 0L
  last.row$chromStart <- unfilled.chromEnd
  last.row$chromEnd <- min.chromEnd
  last.row$count <- 0L
  profile.list[[sample.id]] <-
    rbind(one.sample)
    rbind(first.row, one.sample, last.row)
}
bases <- min.chromEnd-max.chromStart
## End pre-processing to add zeros.

## Small bins are just for testing the computation of the loss
## function in the R implementation, and should not be ported to C
## code.
n.samples <- length(profile.list)

na.mat <- 
  matrix(NA, n.bins, n.samples,
         dimnames=list(bin=NULL, sample.id=names(profile.list)))
first.cumsums <- list(count=na.mat)
bin.list <- list()
norm.list <- list()
for(sample.i in seq_along(profile.list)){
  sample.id <- names(profile.list)[sample.i]
  one <- profile.list[[sample.i]]
  max.count <- max(one$count)
  bins <- binSum(one, max.chromStart, bases.per.bin, n.bins)
  
  stopifnot(n.bins == nrow(bins))
  bins[nrow(bins), "chromEnd"] <- min.chromEnd
  bins$mean <- with(bins, count/(chromEnd-chromStart))
  bins$mean.norm <- bins$mean/max.count
  bin.list[[sample.id]] <- data.frame(sample.id, rbind(bins, NA))
  bases.vec <- with(bins, chromEnd-chromStart)
  stopifnot(bins$count >= 0)
  first.cumsums$count[, sample.i] <- cumsum(bins$count)
  one$count.norm <- one$count/max.count
  norm.list[[sample.i]] <- one
}
bin.df <- do.call(rbind, bin.list)
norm.df <- do.call(rbind, norm.list)

ggplot()+
  scale_color_manual(values=c(data="grey50",
                       bins="black", segments="green"))+
  geom_point(aes(chromEnd, count, color=what),
            data=data.frame(norm.df, what="data"))+
  geom_segment(aes(chromStart+0.5, mean,
                   xend=chromEnd+0.5, yend=mean,
                   color=what),
               data=data.frame(bin.df, what="bins"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })

## The formula for the optimal Poisson loss 
## for 1 segment with d integer data points x_j is
## \sum_{j=1}^d m - x_j \log m_j =
##   ( \sum_{j=1}^d x_j ) (1-\log m)
## where the segment mean m = (\sum x_j)/d,
OptimalPoissonLoss <- function(mean.value, cumsum.value){
  ifelse(mean.value == 0, 0, cumsum.value * (1-log(mean.value)))
}
loss.list <- list()
loss.show.list <- list()
peak.list <- list()
seg.list <- list()
best.loss.list <- list()
flat.cumsums <- first.cumsums$count[n.bins, ]
flat.means <- flat.cumsums/unfilled.bases

flat.loss.vec <- OptimalPoissonLoss(flat.means, flat.cumsums)
best.loss.list[["0"]] <- sum(flat.loss.vec)
for(seg1.last in 1:(n.bins-2)){
  seg1.cumsums <- first.cumsums$count[seg1.last, ]
  seg1.bases <- seg1.last*bases.per.bin
  seg1.chromEnd <- seg1.bases + max.chromStart
  seg1.corrected <- seg1.bases - extra.before
  seg1.means <- seg1.cumsums/seg1.corrected
  seg1.loss.vec <- OptimalPoissonLoss(seg1.means, seg1.cumsums)
  seg1.loss <- sum(seg1.loss.vec)
  for(seg2.last in (seg1.last+1):(n.bins-1)){
    cumsum.seg2.end <- first.cumsums$count[seg2.last, ]
    seg2.cumsums <- cumsum.seg2.end-seg1.cumsums
    seg12.bases <- seg2.last*bases.per.bin
    seg2.bases <- seg12.bases - seg1.bases
    seg2.chromEnd <- seg1.chromEnd + seg2.bases
    seg2.means <- seg2.cumsums/seg2.bases
    seg2.loss.vec <- OptimalPoissonLoss(seg2.means, seg2.cumsums)
    seg2.loss <- sum(seg2.loss.vec)
    
    seg3.cumsums <- first.cumsums$count[n.bins, ]-cumsum.seg2.end
    seg3.bases <- bases - seg12.bases
    seg3.corrected <- seg3.bases - extra.after
    seg3.means <- seg3.cumsums/seg3.corrected

    mean.mat <- rbind(seg1.means, seg2.means, seg3.means)

    this.seg.list <- list()
    for(sample.id in colnames(mean.mat)){
      this.seg.list[[sample.id]] <-
        data.frame(sample.id,
                   chromStart=c(max.chromStart, seg1.chromEnd, seg2.chromEnd,
                     max.chromStart),
                   chromEnd=c(seg1.chromEnd, seg2.chromEnd, min.chromEnd,
                     min.chromEnd),
                   mean=c(mean.mat[, sample.id], flat.means[[sample.id]]),
                   segments=c(3, 3, 3, 1))
    }
    these.segs <- do.call(rbind, this.seg.list)
    ggplot()+
      scale_color_manual(values=c(data="grey50",
                           bins="black", segments="green"))+
      geom_step(aes(chromStart, count, color=what),
                data=data.frame(norm.df, what="data"))+
      geom_segment(aes(chromStart, mean,
                       xend=chromEnd, yend=mean,
                       color=what),
                   size=1,
                   data=data.frame(these.segs, what="segments"))+
      geom_segment(aes(chromStart, mean,
                       xend=chromEnd, yend=mean,
                       color=what),
                   data=data.frame(bin.df, what="bins"))+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
        sub("McGill0", "", val)
      })
    
    seg3.loss.vec <- OptimalPoissonLoss(seg3.means, seg3.cumsums)
    seg3.loss <- sum(seg3.loss.vec)
    seg123.loss.vec <- seg1.loss.vec + seg2.loss.vec + seg3.loss.vec

    peak.feasible <- seg1.means < seg2.means & seg2.means > seg3.means
    diff.loss.vec <- flat.loss.vec - seg123.loss.vec
    possible.df <- 
      data.frame(flat.loss.vec, seg123.loss.vec,
                 diff.loss.vec, peak.feasible)
    ordered.df <- 
      possible.df[order(possible.df$diff.loss.vec, decreasing = TRUE), ]
    for(peaks in 1:nrow(ordered.df)){
      with.peaks <- ordered.df[1:peaks, ]
      with.segs <-
        subset(these.segs,
               sample.id %in% rownames(with.peaks) & segments==3)
      without.segs <-
        subset(these.segs,
               !sample.id %in% rownames(with.peaks) & segments==1)
      without.peaks <-
        possible.df[!rownames(possible.df) %in% rownames(with.peaks),]
      with.loss <- with.peaks$seg123.loss.vec
      without.loss <- without.peaks$flat.loss.vec
      total.loss <- sum(with.loss, without.loss)
      loss.show.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <-
        data.frame(seg1.last, seg2.last, peaks, total.loss,
                   sample.id=bases.per.bin,
                   feasible=all(peak.feasible),
                   chromStart=seg1.last*bases.per.bin+max.chromStart,
                   chromEnd=seg2.last*bases.per.bin+max.chromStart)
    }
    if(any(peak.feasible)){
      feasible.df <- subset(possible.df, peak.feasible)
      ordered.df <- 
        feasible.df[order(feasible.df$diff.loss.vec, decreasing = TRUE), ]
      for(peaks in 1:nrow(ordered.df)){
        with.peaks <- ordered.df[1:peaks, ]
        with.segs <-
          subset(these.segs,
                 sample.id %in% rownames(with.peaks) & segments==3)
        without.segs <-
          subset(these.segs,
                 !sample.id %in% rownames(with.peaks) & segments==1)
        without.peaks <-
          possible.df[!rownames(possible.df) %in% rownames(with.peaks),]
        with.loss <- with.peaks$seg123.loss.vec
        without.loss <- without.peaks$flat.loss.vec
        total.loss <- sum(with.loss, without.loss)
        loss.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <-
          data.frame(seg1.last, seg2.last, peaks, total.loss,
                     sample.id=bases.per.bin,
                     chromStart=seg1.last*bases.per.bin+max.chromStart,
                     chromEnd=seg2.last*bases.per.bin+max.chromStart)
        peak.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <- 
          data.frame(sample.id=rownames(with.peaks),
                     chromStart=seg1.last*bases.per.bin+max.chromStart,
                     chromEnd=seg2.last*bases.per.bin+max.chromStart)
        seg.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <-
          rbind(without.segs, with.segs)
      }#peaks
    }#any(peak.feasible)
  }#seg2.last
}#seg1.last
best.seg.list <- list()
best.peak.list <- list()
best.indices.list <- list()
for(peaks.str in names(loss.list)){
  loss.df <- do.call(rbind, loss.list[[peaks.str]])
  loss.best <- loss.df[which.min(loss.df$total.loss), ]
  best.indices.list[[peaks.str]] <- loss.best
  last.str <- with(loss.best, paste(seg1.last, seg2.last))
  peaks <- as.numeric(peaks.str)
  model.i <- peaks + 1
  peak.df <- peak.list[[peaks.str]][[last.str]]
  seg.df <- seg.list[[peaks.str]][[last.str]]

  ggplot()+
    ggtitle(paste0("best model with ", peaks,
                   " peak", ifelse(peaks==1, "", "s")))+
    scale_color_manual(values=c(data="grey50",
                         bins="black", segments="green"))+
    geom_step(aes(chromStart/1e3, count, color=what),
              data=data.frame(norm.df, what="data"))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 size=1,
                 data=data.frame(seg.df, what="segments"))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 data=data.frame(bin.df, what="bins"))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })

  best.seg.list[[peaks.str]] <- data.frame(peaks, seg.df)
  best.peak.list[[peaks.str]] <-
    data.frame(peaks, y=peaks*-0.1, peak.df)
  best.loss.list[[peaks.str]] <- loss.best$total.loss
}

best.peaks <- do.call(rbind, best.peak.list)
by.sample.loc <-
  split(best.peaks, with(best.peaks, paste(sample.id, chromStart, chromEnd)))
short.label.list <- list()
for(sample.loc.name in names(by.sample.loc)){
  sample.loc <- by.sample.loc[[sample.loc.name]]
  peaks.txt <- paste(sample.loc$peaks, collapse=",")
  short.label.list[[sample.loc.name]] <-
    data.frame(sample.loc[1,], peaks.txt)
}
short.labels <- do.call(rbind, short.label.list)
best.peaks$sample.id <- factor(best.peaks$sample.id, names(profile.list))
sample.counts <- table(best.peaks$sample.id)
dftype <- function(what, df){
  df$sample.id <- factor(df$sample.id, names(sample.counts))
  data.frame(what, df)
}
ggplot()+
  scale_color_manual(values=c(data="grey50",
                       bins="black", peaks="deepskyblue"))+
  geom_step(aes(chromStart/1e3, count.norm, color=what),
            data=dftype("data", norm.df))+
  geom_segment(aes(chromStart/1e3, y,
                   xend=chromEnd/1e3, yend=y,
                   color=what),
               size=1,
               data=dftype("peaks", best.peaks))+
  geom_text(aes(chromStart/1e3, y,
                label=paste0(peaks, " peak",
                  ifelse(peaks==1, "", "s"), " "),
                color=what),
            hjust=1,
            size=3,
            vjust=0.5,
            data=dftype("peaks", best.peaks))+
  geom_segment(aes(chromStart/1e3, mean.norm,
                   xend=chromEnd/1e3, yend=mean.norm,
                   color=what),
               data=dftype("bins", bin.df))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })
## for 0, 1, ..., maxPeaks, run the bin pyramid grid search,
## around the peaks found in this first step.
zoom.peak.list <- list("0"=Peaks())
zoom.loss.list <-
  list("0"=data.frame(peaks=0, loss=sum(flat.loss.vec)))
for(peaks.str in names(best.indices.list)){
  loss.best <- best.indices.list[[peaks.str]]
  best.peak.df <- best.peak.list[[peaks.str]]
  samples.with.peaks <- paste(best.peak.df$sample.id)
  flat.without.peaks <- sum(flat.loss.vec)-
    sum(flat.loss.vec[samples.with.peaks])
  sub.norm.df <- subset(norm.df, sample.id %in% samples.with.peaks)
  sub.bin.df <- subset(bin.df, sample.id %in% samples.with.peaks)
  n.samples <- length(samples.with.peaks)
  last.cumsums <- list()
  before.cumsums <- list(left=list(), right=list())
  for(data.type in names(first.cumsums)){
    data.mat <- first.cumsums[[data.type]]
    last.cumsums[[data.type]] <-
      data.mat[nrow(data.mat),][samples.with.peaks]
    before.cumsums$left[[data.type]] <- if(loss.best$seg1.last == 1){
      structure(rep(0, n.samples), names=samples.with.peaks)
    }else{
      data.mat[loss.best$seg1.last-1,][samples.with.peaks]
    }
    before.cumsums$right[[data.type]] <-
      data.mat[loss.best$seg2.last-1,][samples.with.peaks]
  }
  last.chromEnd <- min.chromEnd
  n.bins.zoom <- bin.factor * 2L
  n.cumsum <- n.bins.zoom + 1L

  ## These will change at the end of each iteration.
  peakStart <- best.peak.df$chromStart[1]
  peakEnd <- best.peak.df$chromEnd[1]
  search.list <- list()
  best.list <- list()
  bases.per.bin.zoom <- bases.per.bin
  while(bases.per.bin.zoom > 1){
    left.chromStart <- peakStart - bases.per.bin.zoom
    right.chromStart <- peakEnd-bases.per.bin.zoom
    bases.per.bin.zoom <- as.integer(bases.per.bin.zoom / bin.factor)

    cumsum.list <- function(chromStart){
      limits <- bases.per.bin.zoom*(0:(n.cumsum-1))+chromStart
      chromStart.vec <- limits[-length(limits)]
      chromEnd.vec <- limits[-1]
      intervals <-
        paste0(chromStart.vec, "-", chromEnd.vec)
      dn <- list(bin=c("before", intervals), sample.id=samples.with.peaks)
      m <- matrix(NA, n.cumsum, n.samples, dimnames=dn)
      list(count=m, 
           chromStart=chromStart.vec, chromEnd=chromEnd.vec)
    }
    cumsum.mats <-
      list(left=cumsum.list(left.chromStart),
           right=cumsum.list(right.chromStart))

    for(sample.i in seq_along(samples.with.peaks)){
      sample.id <- samples.with.peaks[[sample.i]]
      one <- profile.list[[sample.id]]
      lr.list <-
        list(left=binSum(one, left.chromStart,
               bases.per.bin.zoom, n.bins.zoom),
             right=binSum(one, right.chromStart,
               bases.per.bin.zoom, n.bins.zoom,
               empty.as.zero=TRUE))
      for(lr in names(lr.list)){
        lr.bins <- lr.list[[lr]]
        stopifnot(nrow(lr.bins) == n.bins.zoom)
        lr.bases <- with(lr.bins, chromEnd-chromStart)
        lr.before <- before.cumsums[[lr]]
        lr.counts <- list(count=lr.bins$count)
        for(data.type in names(lr.counts)){
          lr.count.vec <-
            c(lr.before[[data.type]][[sample.id]], lr.counts[[data.type]])
          cumsum.mats[[lr]][[data.type]][, sample.id] <-
            cumsum(lr.count.vec)
        }
      }
    }
    ##print(t(cumsum.mats$left$count[-1,]))
    possible.grid <- 
      expand.grid(left.cumsum.row=2:n.cumsum, right.cumsum.row=2:n.cumsum)
    possible.grid$left.chromStart <-
      cumsum.mats$left$chromStart[possible.grid$left.cumsum.row-1]
    possible.grid$left.chromEnd <-
      cumsum.mats$left$chromEnd[possible.grid$left.cumsum.row-1]
    possible.grid$right.chromStart <-
      cumsum.mats$right$chromStart[possible.grid$right.cumsum.row-1]
    possible.grid$right.chromEnd <-
      cumsum.mats$right$chromEnd[possible.grid$right.cumsum.row-1]
    feasible.grid <-
      subset(possible.grid,
             left.chromStart <= right.chromStart)
    feasible.grid$model.i <- 1:nrow(feasible.grid)
    model.list <- list()
    seg.list <- list()
    sample.loss.list <- list()
    for(model.i in feasible.grid$model.i){
      model.row <- feasible.grid[model.i, ]

      seg1.i <- model.row$left.cumsum.row-1
      seg1.cumsums <- cumsum.mats$left$count[seg1.i, ]
      seg1.chromEnd <- cumsum.mats$left$chromStart[seg1.i]
      seg1.bases <- seg1.chromEnd-max.chromStart
      seg1.corrected <- seg1.bases - extra.before
      seg1.means <- seg1.cumsums/seg1.corrected
      seg1.loss <- OptimalPoissonLoss(seg1.means, seg1.cumsums)
      seg.list[[paste(model.i, 1)]] <-
        data.frame(chromStart=unfilled.chromStart, chromEnd=seg1.chromEnd,
                   mean=seg1.means, sample.id=samples.with.peaks,
                   model.i,
                   row.names=NULL)

      cumsum.seg2.end <-
        cumsum.mats$right$count[model.row$right.cumsum.row, ]
      seg2.cumsums <- cumsum.seg2.end-seg1.cumsums
      seg2.chromEnd <-
        cumsum.mats$right$chromEnd[model.row$right.cumsum.row-1]
      seg2.bases <- seg2.chromEnd-seg1.chromEnd
      seg2.means <- seg2.cumsums/seg2.bases
      seg2.loss <- OptimalPoissonLoss(seg2.means, seg2.cumsums)
      seg.list[[paste(model.i, 2)]] <-
        data.frame(chromStart=seg1.chromEnd, chromEnd=seg2.chromEnd,
                   mean=seg2.means, sample.id=samples.with.peaks,
                   model.i,
                   row.names=NULL)

      seg3.cumsums <- last.cumsums$count-cumsum.seg2.end
      seg3.bases <- last.chromEnd-seg2.chromEnd
      seg3.corrected <- seg3.bases - extra.after
      seg3.means <- seg3.cumsums/seg3.corrected
      seg3.loss <- OptimalPoissonLoss(seg3.means, seg3.cumsums)
      seg.list[[paste(model.i, 3)]] <-
        data.frame(chromStart=seg2.chromEnd, chromEnd=unfilled.chromEnd,
                   mean=seg3.means, sample.id=samples.with.peaks,
                   model.i,
                   row.names=NULL)

      total.bases <- sum(seg1.bases + seg2.bases + seg3.bases)
      
      total.loss <- sum(seg1.loss + seg2.loss + seg3.loss)
      with.without.loss <- flat.without.peaks + total.loss
      total.loss.vec <- seg1.loss+seg2.loss+seg3.loss
      feasible.vec <- seg1.means < seg2.means & seg2.means > seg3.means
      model.list[[model.i]] <-
        data.frame(model.row, total.loss, with.without.loss,
                   feasible=all(feasible.vec))
      sample.loss.list[[model.i]] <-
        data.frame(sample.id=samples.with.peaks, loss=total.loss.vec)
    }
    ## Plot the segment means as a reality check.
    seg.df <- do.call(rbind, seg.list)
    ggplot()+
      geom_step(aes(chromStart/1e3, count),
                data=data.frame(sub.norm.df, what="data"),
                color="grey")+
      geom_segment(aes(chromStart/1e3, mean,
                       xend=chromEnd/1e3, yend=mean),
                   data=data.frame(seg.df, what="models"),
                   size=1, color="green")+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ model.i, scales="free")

    ## Then plot the peaks only, colored by total cost of the model.
    model.df <- do.call(rbind, model.list)
    model.df$y <- with(model.df, model.i/max(model.i))
    feasible.models <- subset(model.df, feasible)
    best.model <- feasible.models[which.min(feasible.models$total.loss), ]

    ggplot()+
      coord_equal()+
      geom_point(aes(with.without.loss, total.loss),
                 data=model.df)

    ggplot()+
      xlab("position on chromosome (kilobases = kb)")+
      scale_y_continuous("", breaks=NULL)+
      geom_step(aes(chromStart/1e3, count.norm),
                data=data.frame(sub.norm.df, what="data"),
                color="grey")+
      geom_segment(aes(chromStart/1e3, mean.norm,
                       xend=chromEnd/1e3, yend=mean.norm),
                   data=data.frame(sub.bin.df, what="bins"),
                   color="black")+
      geom_segment(aes(peakStart/1e3, 0,
                       xend=peakEnd/1e3, yend=0),
                   data=data.frame(loss.best, what="peak"),
                   color="green")+
      scale_linetype_manual(values=c("TRUE"=1, "FALSE"=2))+
      geom_segment(aes(left.chromStart/1e3, y,
                       color=total.loss,
                       linetype=feasible,
                       xend=right.chromEnd/1e3, yend=y),
                   data=data.frame(model.df, sample.id="peaks"),
                   size=1)+
      geom_text(aes(left.chromStart/1e3, y,
                    label="optimal "),
                data=data.frame(best.model, sample.id="peaks"),
                hjust=1)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ .)

    search.list[[paste(bases.per.bin.zoom)]] <-
      data.frame(bases.per.bin.zoom, model.df)
    best.list[[paste(bases.per.bin.zoom)]] <-
      data.frame(bases.per.bin.zoom, best.model)

    peakStart <- best.model$left.chromStart
    peakEnd <- best.model$right.chromEnd
    before.i.list <-
      list(left=best.model$left.cumsum.row-2,
           right=best.model$right.cumsum.row-1)
    for(lr in names(before.i.list)){
      before.i <- before.i.list[[lr]]
      mats <- cumsum.mats[[lr]]
      before.cumsums[[lr]]$count <-
        structure(mats$count[before.i,],
                  names=colnames(mats$count))
    }
  }##while(bases.per.bin.zoom)
  samplefac <- function(L){
    df <- do.call(rbind, L)
    sample.ids <- unique(norm.df$sample.id)
    bases.num <- sort(unique(df$bases.per.bin.zoom), decreasing=TRUE)
    bases.num <- c(bases.num[1]*2, bases.num)
    levs <- c(paste(sample.ids), paste(bases.num))
    df$sample.id <- factor(df$bases.per.bin.zoom, levs)
    df
  }
  search.df <- samplefac(search.list)
  best.df <- samplefac(best.list)
  
  loss.with.peaks <- sample.loss.list[[best.model$model.i]]
  samples.without.peaks <-
    names(profile.list)[!names(profile.list) %in% samples.with.peaks]
  loss.without.peaks <-
    data.frame(sample.id=samples.without.peaks,
               loss=flat.loss.vec[samples.without.peaks])
  sample.loss.df <- rbind(loss.with.peaks, loss.without.peaks)
  peaks <- as.numeric(peaks.str)
  zoom.loss.list[[peaks.str]] <- 
    data.frame(peaks, loss=sum(sample.loss.df$loss))
  zoom.peak.list[[peaks.str]] <-
    data.frame(peaks, sample.id=samples.with.peaks,
               chromStart=peakStart, chromEnd=peakEnd)
}#peaks.str
zoom.peaks <- do.call(rbind, zoom.peak.list)
zoom.loss <- do.call(rbind, zoom.loss.list)

first.loss <- do.call(rbind, loss.list[["2"]])
first.loss$model.i <- 1:nrow(first.loss)
first.loss$sample.id <-
  factor(first.loss$sample.id, c("sample1", "sample2", "4", "2", "1"))

loss.show <- do.call(rbind, loss.show.list[["2"]])
loss.show$sample.id <- 
  factor(4, c("sample1", "sample2", "4", "2", "1"))
loss.show$model.i <- 1:nrow(loss.show)
loss.show.feas <- subset(loss.show, feasible)
loss.show.ord <- loss.show.feas[order(loss.show.feas$total.loss), ][1,]

best.segs <- subset(seg.df, model.i == best.model$model.i)
best.breaks <-
  subset(best.segs,
         min(chromStart) < chromStart &
          sample.id == sample.id[1],
         -sample.id)
makePeak <- function(start, end){
  sprintf("(%d,%d]", start, end)
}
best.df$peak <- with(best.df, makePeak(left.chromStart, right.chromEnd))
loss.show.ord$peak <- with(loss.show.ord, makePeak(chromStart, chromEnd))
loss.show$peak <- with(loss.show, makePeak(chromStart, chromEnd))
search.df$peak <- with(search.df, makePeak(left.chromStart, right.chromEnd))

viz.peaks <-
  data.frame(chromStart=c(loss.show$chromStart, search.df$left.chromStart),
             chromEnd=c(loss.show$chromEnd, search.df$right.chromEnd))
viz.segs.list <- list()
for(peak.i in 1:nrow(viz.peaks)){
  p <- viz.peaks[peak.i, ]
  p.vec <- unlist(p)
  peak <- with(p, makePeak(chromStart, chromEnd))
  seg.end.vec <- c(p.vec, unfilled.chromEnd)
  seg.start.vec <- c(0, p.vec)
  for(sample.id in names(profile.list)){
    pro <- profile.list[[sample.id]]
    for(seg.i in seq_along(seg.end.vec)){
      seg.end <- seg.end.vec[[seg.i]]
      seg.start <- seg.start.vec[[seg.i]]+1
      seg.data <- pro$count[seg.start:seg.end]
      seg.mean <- mean(seg.data)
      viz.segs.list[[paste(peak.i, sample.id, seg.i)]] <- 
        data.frame(chromStart=seg.start-1,
                   chromEnd=seg.end,
                   sample.id,
                   peak,
                   mean=seg.mean)
    }
  }
}
viz.segs <- unique(do.call(rbind, viz.segs.list))
viz.breaks <-
  subset(viz.segs,
         min(chromStart) < chromStart &
         sample.id == sample.id[1], select=-sample.id)

viz <-
  list(panels=ggplot()+
         scale_x_continuous("data to segment (base position along chromosome)",
                            breaks=1:24)+
         scale_color_continuous("Poisson loss")+
         scale_y_continuous("", breaks=NULL)+
         geom_point(aes(chromEnd, count),
                    data=data.frame(sub.norm.df, what="data"),
                    color="black")+
         geom_segment(aes(chromStart+0.5, model.i,
                          clickSelects=peak,
                          xend=chromEnd+0.5, yend=model.i,
                          linetype=ifelse(feasible, "feasible", "infeasible"),
                          color=total.loss),
                      size=4,
                      data=loss.show)+
         geom_segment(aes(left.chromStart+0.5, y,
                          clickSelects=peak,
                          color=total.loss,
                          linetype=ifelse(feasible, "feasible", "infeasible"),
                          xend=right.chromEnd+0.5, yend=y),
                      data=search.df,
                      size=4)+
         scale_linetype_discrete("feasible?")+
         ##guides(linetype=guide_legend(order=2))+
         geom_text(aes(right.chromEnd+0.5, y-0.05,
                       clickSelects=peak,
                       label=" selected"),
                   data=best.df,
                   ##size=3,
                   hjust=0)+
         geom_text(aes(chromStart+0.5, model.i-0.5,
                       clickSelects=peak,
                       label="selected "),
                   data=loss.show.ord,
                   ##size=3,
                   hjust=1)+
         theme_bw()+
         theme(panel.margin=grid::unit(0, "cm"))+
         facet_grid(sample.id ~ ., labeller=function(var, val){
           ifelse(grepl("sample", val), paste(val), paste("bin size", val))
         }, scales="free")+
         geom_segment(aes(chromStart+0.5, mean,
                          showSelected=peak,
                          xend=chromEnd+0.5, yend=mean),
                      data=viz.segs,
                      color="green")+
         geom_vline(aes(xintercept=chromStart+0.5,
                        showSelected=peak),
                    data=viz.breaks,
                    color="green",
                    linetype="dotted")+
         ggtitle("click blue models to select peak start/end")+
         theme_animint(width=1000, height=600),

       duration=list(peak=500),

       first=list(peak=subset(best.df, sample.id==1)$peak),

       title="PeakSegJoint fast heuristic segmentation algorithm")

animint2dir(viz, "figure-heuristic-algo")

fig <-
  ggplot()+
    scale_x_continuous("data to segment (base position along chromosome)",
                       breaks=c(1,seq(2, 24, by=2)))+
    scale_color_continuous("Poisson loss")+
    geom_segment(aes(chromStart+0.5, model.i,
                     xend=chromEnd+0.5, yend=model.i,
                     linetype=ifelse(feasible, "feasible", "infeasible"),
                     color=total.loss),
                 size=1,
                 data=loss.show)+
    scale_y_continuous("", breaks=NULL)+
    geom_point(aes(chromEnd, count),
               data=data.frame(sub.norm.df, what="data"),
               color="black")+
    geom_segment(aes(left.chromStart+0.5, y,
                     color=total.loss,
                     linetype=ifelse(feasible, "feasible", "infeasible"),
                     xend=right.chromEnd+0.5, yend=y),
                 data=search.df,
                 size=1)+
    scale_linetype_discrete("feasible?")+
    ##guides(linetype=guide_legend(order=2))+
    geom_text(aes(right.chromEnd+0.5, y,
                  label=" selected"),
              data=best.df,
              size=3,
              hjust=0)+
    geom_text(aes(chromStart+0.5, model.i,
                  label="selected "),
              data=loss.show.ord,
              size=3,
              hjust=1)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      ifelse(grepl("sample", val), paste(val), paste("bin size", val))
    }, scales="free")+
    geom_segment(aes(chromStart+0.5, mean,
                     xend=chromEnd+0.5, yend=mean),
                 data=best.segs,
                 color="green")+
    geom_vline(aes(xintercept=chromStart+0.5),
               data=best.breaks,
               color="green",
               linetype="dotted")

pdf("figure-heuristic-algo.pdf", h=5)
print(fig)
dev.off()
