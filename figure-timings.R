works_with_R("3.2.0",
             tikzDevice="0.7.0",
             ggplot2="1.0")

load("timings.RData")

by.expr <- split(timings$seconds, timings$seconds$expr)
last <- function(algo){
  df <- by.expr[[algo]]
  df <- subset(df, n.data == max(n.data))
  mean(df$time)/1e9
}

label.df <-
  data.frame(data="seconds",
             algorithm=c("cDPA", "pDPA", "PeakSegJoint"),
             complexity=c("$O(B^2)$", "$O(B\\log B)$", ""),
             n.data=4000,
             seconds=c(last("cDPA"), last("pDPA"), 0.01))
gg <- 
ggplot()+
  scale_x_log10("data size to segment $B$")+
  scale_y_log10("")+
  geom_point(aes(n.data, time/1e9, color=expr),
             pch=1,
             data=data.frame(timings$seconds, data="seconds"))+
  scale_size_manual(values=c(PeakSegJoint=2, cDPA=2, pDPA=1))+
  geom_line(aes(n.data, diff, color=algorithm, size=algorithm),
            data=data.frame(timings$results, data="distance to true peak"))+
  theme_bw()+
  guides(size="none", color="none")+
  geom_text(aes(n.data, seconds, label=paste(algorithm, complexity),
                color=algorithm),
            show_guide=FALSE,
            size=3,
            hjust=0,
            data=label.df)+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(data ~ ., scales="free")

options(tikzMetricsDictionary="tikzMetrics",
        tikzDocumentDeclaration="\\documentclass{article}")
tikz("figure-timings.tex", h=3, w=6)
print(gg)
dev.off()
