works_with_R("3.1.3",
             ggplot2="1.0",
             "tdhock/PeakSegJoint@361a3f1a9037947a1c79a3e754fded04f712008d",
             microbenchmark="1.3.0")

data(H3K36me3.TDH.other.chunk1)
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         43100000 < chromEnd & chromStart < 43205000)
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
    fit <- PeakSegJointHeuristic(some.counts, PINT)
    model <- fit$models[[3]]
    peak <- model$peak_start_end
    peak.df <- data.frame(chromStart=peak[1], chromEnd=peak[2])
    loss <- model$loss
    results[[paste(PINT)]] <- data.frame(param=PINT, peak.df, loss)
  }, list(PINT=param.int))
}
set.seed(1)
times <- microbenchmark(list=m.args, times=3)
results.df <- do.call(rbind, results)

ggplot()+
  ##coord_equal()+
  geom_text(aes(chromStart, chromEnd, label=param),
             data=results.df)

param.plot <- 
ggplot()+
  ylab("")+
  geom_segment(aes(param, chromStart/1e3,
                   xend=param, yend=chromEnd/1e3),
               data=data.frame(results.df, what="peak position (kb)"))+
  scale_x_log10("bin factor suboptimality parameter")+
  ##scale_y_log10()+
  facet_grid(what ~ ., scales = "free")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_point(aes(as.numeric(as.character(expr)), log10(time/1e9)),
             data=data.frame(times, what="log10(seconds)"), pch=1)

pdf("figure-bin-factor.pdf", h=4)
print(param.plot)
dev.off()
