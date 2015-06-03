works_with_R("3.2.0",
             "tdhock/PeakSegJoint@86bee0a4620160e2d4b904e7819b5792280d51de",
             data.table="1.9.4",
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
             complexity=c("$O(B^2)$", "$O(B\\log B)$", "$O(B\\log B)$"),
             n.data=c(8000, 8000, 2000),
             seconds=c(last("cDPA"), last("pDPA"), 0.05))

timings$problems[, bases := problemEnd-problemStart]
problem.range <- timings$problems[, .(min=min(bases), max=max(bases))]

gg <- 
ggplot()+
  geom_tallrect(aes(xmin=min, xmax=max),
                data=problem.range,
                fill="black",
                alpha=0.1)+
  scale_x_log10("data size to segment $B$")+
  scale_y_log10("seconds")+
  geom_point(aes(n.data, time/1e9, color=expr),
             pch=1,
             data=data.frame(timings$seconds, data="seconds"))+
  scale_size_manual(values=c(PeakSegJoint=2, cDPA=2, pDPA=1))+
  theme_bw()+
  guides(size="none", color="none")+
  geom_text(aes(n.data, seconds, label=paste(algorithm, complexity, sep="\n"),
                color=algorithm),
            show_guide=FALSE,
            size=3,
            vjust=1,
            hjust=0,
            data=label.df)+
  ##facet_grid(data ~ ., scales="free")+
  theme(panel.margin=grid::unit(0, "cm"))

options(tikzMetricsDictionary="tikzMetrics",
        tikzDocumentDeclaration="\\documentclass{article}\\usepackage{nips15submit_e,times}")
tikz("figure-timings-small.tex", h=2, w=4)
print(gg)
dev.off()

gg <- 
ggplot()+
  geom_tallrect(aes(xmin=min, xmax=max),
                data=problem.range,
                fill="black",
                alpha=0.1)+
  ## geom_vline(aes(xintercept=problemEnd-problemStart),
  ##            data=timings$problems)+
  scale_x_log10("data size to segment $B$")+
  scale_y_log10("seconds")+
  geom_point(aes(n.data, time/1e9, color=expr),
             pch=1,
             data=data.frame(timings$seconds, data="seconds"))+
  scale_size_manual(values=c(PeakSegJoint=2, cDPA=2, pDPA=1))+
  ## geom_line(aes(n.data, diff, color=algorithm, size=algorithm),
  ##           data=data.frame(timings$results, data="distance to\ntrue peak"))+
  ## geom_line(aes(n.data, 10^mean.loss, color=algorithm, size=algorithm),
  ##           data=data.frame(timings$results, data="mean\nPoisson loss"))+
  theme_bw()+
  guides(size="none", color="none")+
  geom_text(aes(n.data, seconds, label=paste(algorithm, complexity),
                color=algorithm),
            show_guide=FALSE,
            size=3,
            vjust=1,
            hjust=0,
            data=label.df)+
  ##facet_grid(data ~ ., scales="free")+
  theme(panel.margin=grid::unit(0, "cm"))

tikz("figure-timings.tex", h=2, w=5)
print(gg)
dev.off()
