works_with_R("3.1.3",
             ggplot2="1.0",
             "tdhock/PeakSegJoint@e6bc386c555e203cc80343814939d51785c03af1",
             microbenchmark="1.3.0")

data(H3K36me3.TDH.other.chunk1)
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         43100000 < chromEnd & chromStart < 43205000)
count.list <- ProfileList(some.counts)
##gctorture(TRUE)
m.args <- list()
results <- list()
param.values <-
  c(2, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100, 125, 150,
    200, 300, 400, 500, 600)
for(param in param.values){
  param.name <- paste(param)
  param.int <- as.integer(param)
  m.args[[param.name]] <- substitute({
    print(PINT)
    fit <- PeakSegJointHeuristic(count.list, PINT)
    model <- fit$models[[3]]
    peak <- model$peak_start_end
    peak.df <- data.frame(chromStart=peak[1], chromEnd=peak[2])
    results[[paste(PINT)]] <- data.frame(param=PINT, peak.df)
  }, list(PINT=param.int))
}
set.seed(1)
times <- microbenchmark(list=m.args, times=3)

for(param.name in names(results)){
  result.row <- results[[param.name]]
  start.end.vec <-
    c(seg1.chromStart=fit$data_start_end[1],
      peakStart=result.row$chromStart,
      peakEnd=result.row$chromEnd,
      seg3.chromEnd=fit$data_start_end[2])
  seg.bases <- diff(start.end.vec)
  seg.starts <- start.end.vec[1:3]
  seg.list <- list()
  for(sample.id in names(count.list)){
    sample.counts <- count.list[[sample.id]]
    for(seg.i in seq_along(seg.bases)){
      bases.int <- as.integer(seg.bases[[seg.i]])
      start.int <- as.integer(seg.starts[[seg.i]])
      bin <- binSum(sample.counts, start.int, bases.int, 1L)
      loss <- bin$count*(1-log(bin$mean))
      seg.list[[paste(sample.id, seg.i)]] <-
        data.frame(sample.id, seg.i, loss)
    }
  }
  all.segs <- do.call(rbind, seg.list)
  results[[param.name]]$loss <- sum(all.segs$loss)
}

results.df <- do.call(rbind, results)

ggplot()+
  ##coord_equal()+
  geom_text(aes(chromStart, chromEnd, label=param),
             data=results.df)

min.loss <- subset(results.df, loss==min(loss))

refs <- data.frame(seconds=1, label="1 second", what="log10(seconds)")

param.plot <- 
ggplot()+
  ylab("")+
  geom_segment(aes(param, chromStart/1e3,
                   xend=param, yend=chromEnd/1e3),
               data=data.frame(results.df, what="peak position (kb)"))+
  geom_line(aes(param, loss),
            data=data.frame(results.df, what="Poisson Loss"))+
  geom_point(aes(param, loss),
            data=data.frame(min.loss, what="Poisson Loss"))+
  scale_x_log10("bin factor suboptimality parameter",
                breaks=c(2, 3, 5, 10, 100, max(param.values)))+
  ##scale_y_log10()+
  facet_grid(what ~ ., scales = "free")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_hline(aes(yintercept=log10(seconds)),
             color="grey",
             data=refs)+
  geom_text(aes(10, log10(seconds), label=label),
            vjust=-0.5,
            color="grey",
             data=refs)+
  geom_point(aes(as.numeric(as.character(expr)), log10(time/1e9)),
             data=data.frame(times, what="log10(seconds)"), pch=1)

pdf("figure-bin-factor.pdf", h=4)
print(param.plot)
dev.off()
