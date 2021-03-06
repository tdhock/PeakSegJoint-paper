works_with_R("3.2.0",
             ggplot2="1.0",
             microbenchmark="1.3.0",
             "tdhock/PeakSegDP@fe06a5b91d68c5d1ec471cb15c3ec3935dc2624d",
             "tdhock/PeakSegJoint@ff5a7c58e297b54b328047f4e02285f0cb5d2838")

data(H3K36me3.TDH.other.chunk1)

some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         43379893 < chromEnd)
profile.list <- split(some.counts, some.counts$sample.id)
sample.peak.list <- list()

param.values <-
  c(2, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100, 200, 400, 800)
cDPA.seconds.list <- list()
heuristic.seconds.list <- list()
loss.list <- list()
for(sample.i in seq_along(profile.list)){
  sample.id <- names(profile.list)[[sample.i]]
  cat(sprintf("%4d / %4d %s\n", sample.i, length(profile.list), sample.id))
  sample.counts <- profile.list[[sample.id]]
  m.args <- list()
  result.list <- list()
  for(param in param.values){
    param.name <- paste(param)
    param.int <- as.integer(param)
    m.args[[param.name]] <- substitute({
      fit <- PeakSegJointHeuristic(sample.counts, PINT)
      model <- fit$models[[2]]
      peak <- model$peak_start_end
      peak.df <- data.frame(chromStart=peak[1], chromEnd=peak[2])
      result.list[[paste(PINT)]] <- data.frame(param=PINT, peak.df)
    }, list(PINT=param.int))
  }
  times <- microbenchmark(list=m.args, times=3)
  results <- do.call(rbind, result.list)
  times$seconds <- times$time/1e9
  heuristic.seconds.list[[sample.id]] <-
    data.frame(algorithm="heuristic", sample.id, times)

  heuristic <- result.list[["2"]]
  dp.seconds <- system.time({
    dp <- PeakSegDP(sample.counts, maxPeaks=1L)
  })[["elapsed"]]
  cDPA.seconds.list[[sample.id]] <-
    data.frame(algorithm="cDPA", sample.id, seconds=dp.seconds)

  dp.peak <- dp$peaks[["1"]]

  all.peaks <-
    rbind(data.frame(algorithm="heuristic", results),
          data.frame(algorithm="cDPA",
                     param=NA,
                     chromStart=dp.peak$chromStart,
                     chromEnd=dp.peak$chromEnd))

  count.vec <- with(sample.counts, rep(count, chromEnd-chromStart))
  sample.bases <- sum(with(sample.counts, chromEnd-chromStart))
  stopifnot(all.equal(sample.bases, length(count.vec)))
  seg1.chromStart <- sample.counts$chromStart[1]
  seg3.chromEnd <- sample.counts[nrow(sample.counts), "chromEnd"]
  seg1.chromEnd <- all.peaks$chromStart
  seg2.chromEnd <- all.peaks$chromEnd
  seg1.first <- seg1.chromStart - seg1.chromStart + 1
  seg1.last <- seg1.chromEnd - seg1.chromStart
  seg2.first <- seg1.last + 1
  seg2.last <- seg2.chromEnd - seg1.chromStart
  seg3.first <- seg2.last + 1
  seg3.last <- length(count.vec)
  data.sizes <- rowSums(cbind(seg3.last-seg2.last,
                              seg2.last-seg1.last,
                              seg1.last))
  stopifnot(data.sizes == length(count.vec))
  all.peaks$loss <- NA
  for(model.i in 1:nrow(all.peaks)){
    model <- all.peaks[model.i, ]
    seg.indices <-
      data.frame(first=c(seg1.first, seg2.first[model.i], seg3.first[model.i]),
                 last=c(seg1.last[model.i], seg2.last[model.i], seg3.last))
    seg.loss.list <- list()
    for(seg.i in 1:nrow(seg.indices)){
      s <- seg.indices[seg.i, ]
      seg.data <- count.vec[s$first:s$last]
      seg.mean <- mean(seg.data)
      seg.loss.list[[seg.i]] <- PoissonLoss(seg.data, seg.mean)
    }
    all.peaks$loss[model.i] <- sum(unlist(seg.loss.list))
  }
  all.peaks$status <- with(all.peaks, ifelse(loss == min(loss), "best", "not"))
  best.peak <- subset(all.peaks, loss == min(loss))[1, ]
  all.peaks$chromStart.diff <- all.peaks$chromStart - best.peak$chromStart
  all.peaks$chromEnd.diff <- all.peaks$chromEnd - best.peak$chromEnd

  loss.list[[sample.id]] <- data.frame(sample.id, all.peaks)
  
  sample.peak.list[[sample.id]] <-
    data.frame(sample.id,
               y=-0.1 * (1:2) * max(sample.counts$count),
               algorithm=c("heuristic.2", "cDPA"),
               chromStart=c(heuristic$chromStart, dp.peak$chromStart),
               chromEnd=c(heuristic$chromEnd, dp.peak$chromEnd))
}
sample.peaks <- do.call(rbind, sample.peak.list)
h.prof <- 
ggplot()+
  scale_color_discrete(limits=c("heuristic.2", "cDPA"))+
  geom_step(aes(chromStart/1e3, count),
            data=some.counts,
            color="grey50")+
  geom_segment(aes(chromStart/1e3, y,
                   color=algorithm,
                   xend=chromEnd/1e3, yend=y),
               data=sample.peaks,
               size=2)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })
png("figure-heuristic-profiles.png", width=9, h=7, res=200, units="in")
print(h.prof)
dev.off()

cDPA.seconds <- do.call(rbind, cDPA.seconds.list)
heuristic.seconds <- do.call(rbind, heuristic.seconds.list)
h.times <- 
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., labeller=function(var, val){
    sub("McGill0", "", val)
  })+
  geom_point(aes(as.numeric(as.character(expr)), seconds, color=algorithm),
             data=heuristic.seconds,
             pch=1)+
  scale_x_log10("bin factor suboptimality parameter",
                breaks=c(2, 10, 100, max(param.values)))+
  scale_y_log10()+
  geom_hline(aes(yintercept=seconds, color=algorithm),
             data=cDPA.seconds)
pdf("figure-heuristic-times.pdf", h=6)
print(h.times)
dev.off()

loss <- do.call(rbind, loss.list)
param.num <- sort(unique(as.numeric(loss$param)))
levs <- c("cDPA", param.num)
loss$param.fac <-
  factor(ifelse(is.na(loss$param), "cDPA", loss$param), levs)
h.loss <- 
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id + what ~ ., labeller=function(var, val){
    sub("McGill0", "", val)
  }, scales="free")+
  xlab("bin factor suboptimality parameter")+
  scale_size_manual(values=c(best=3, not=1))+
  geom_point(aes(param.fac, loss, color=algorithm, size=status),
             data=data.frame(loss, what="loss"))+
  geom_segment(aes(param.fac, chromStart,
                   xend=param.fac, yend=chromEnd,
                   color=algorithm, size=status),
               data=data.frame(loss, what="peak"))+
  geom_text(aes(param.fac, chromStart,
                label=chromStart.diff),
            vjust=0,
            size=2,
            data=data.frame(loss, what="peak"))+
  geom_text(aes(param.fac, chromEnd,
                label=chromEnd.diff),
            vjust=1,
            size=2,
            data=data.frame(loss, what="peak"))+
  ggtitle("distance to best peak in bases")

pdf("figure-heuristic-loss.pdf", h=6)
print(h.loss)
dev.off()
