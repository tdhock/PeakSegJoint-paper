works_with_R("3.1.3",
             ##"tdhock/PeakSegJoint@547ca81ce92c38b7c7f78cd071efc5afa96cb289",
             microbenchmark="1.3.0")

library(microbenchmark)
library(PeakSegJoint)

data(H3K36me3.TDH.other.chunk1)
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         43100000 < chromEnd & chromStart < 43205000)

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
optimal.timing <-
  data.frame(seconds=optimal.seconds, what="seconds")

library(ggplot2)
ggplot()+
  ggtitle(paste("heuristic segmentation accurate within 10 bases",
                "and much faster than optimal segmentation",
                sep="\n"))+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(as.numeric(as.character(param)), diff.bases),
             data=data.frame(results.df, what="diff.bases"),
             pch=1)+
  facet_grid(what ~ ., scales="free")+
  geom_hline(aes(yintercept=seconds), data=optimal.timing)+
  geom_text(aes(10, seconds,
                label=sprintf("optimal segmentation = %.1f seconds",
                  seconds)),
            hjust=0,
            vjust=1.5,
            data=optimal.timing)+
  geom_point(aes(as.numeric(as.character(expr)), time/1e9),
             data=data.frame(times, what="seconds"), pch=1)
