works_with_R("3.2.0",
             "tdhock/PeakSegJoint@191489bb2c754da8abe1396f9812490ef7b11eed",
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
expr2algo <-
  c(cDPA="cDPA PeakSeg",
    pDPA="pDPA",
    PeakSegJoint="JointZoom PeakSegJoint")
timings$seconds$algorithm <- expr2algo[paste(timings$seconds$expr)]

label.df <-
  data.frame(data="seconds",
             algorithm=c("cDPA PeakSeg", "pDPA", "JointZoom PeakSegJoint"),
             complexity=c("$O(B^2)$", "$O(B\\log B)$", "$O(B\\log B)$"),
             n.data=c(8000, 8000, 500),
             seconds=c(last("cDPA"), last("pDPA"), 1e-4))

timings$problems[, bases := problemEnd-problemStart]
problem.range <- timings$problems[, .(min=min(bases), max=max(bases))]

labfun <- function(val){
  paste(val)
}
labfun <- scales::scientific

gg <- 
ggplot()+
  geom_tallrect(aes(xmin=min, xmax=max),
                data=problem.range,
                fill="black",
                alpha=0.1)+
  scale_x_log10("data size to segment $B$", labels=labfun)+
  scale_y_log10("seconds", limits=c(5e-5, 2e1),labels=labfun)+
  geom_point(aes(n.data, time/1e9, color=algorithm),
             pch=1,
             data=data.frame(timings$seconds, data="seconds"))+
  scale_size_manual(values=c("JointZoom PeakSegJoint"=2,
                      "cDPA PeakSeg"=2, pDPA=1))+
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

options(tikzMetricsDictionary="tikzMetrics",
        tikzDocumentDeclaration="\\documentclass{article}\\usepackage{nips15submit_e,times}")
tikz("figure-timings-small.tex", h=2, w=4)
print(gg)
dev.off()


label.df <-
  data.frame(data="seconds",
             algorithm=c("cDPA PeakSeg", "pDPA", "JointZoom PeakSegJoint"),
             complexity=c("$O(B^2)$", "$O(B\\log B)$", "$O(B\\log B)$"),
             n.data=c(8000, 8000, 5e3),
             seconds=c(last("cDPA"), last("pDPA"), 5e-4))

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
  geom_point(aes(n.data, time/1e9, color=algorithm),
             pch=1,
             data=data.frame(timings$seconds, data="seconds"))+
  scale_size_manual(values=c("JointZoom PeakSegJoint"=2,
                      "cDPA PeakSeg"=2, pDPA=1))+
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

tikz("figure-timings.tex", h=1.7, w=5)
print(gg)
dev.off()
