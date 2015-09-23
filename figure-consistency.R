works_with_R("3.2.2",
             data.table="1.9.6",
             ggplot2="1.0.1",
             "tdhock/PeakSegJoint@55d3e66c77e1aee43f008e0c4e0c3ac35b21d4e4",
             "tdhock/PeakSegDP@6520faad59fcf8c94ed4345d2b0d826b7f61faf9")

color.code <-
  c(truth="#1B9E77", #teal
    PeakSeg="#D95F02", #orange
    PeakSegJoint="#7570B3", #violet
    "#E7298A", #pink
    "#66A61E", #green
    "#E6AB02", #tan
    "#A6761D", #brown
    "#666666") #grey

data.per.segment <- 100
mu.low <- 10
mu.high <- 11
n.samples <- 10
set.seed(1)
signal.dt.list <- list()
model.dt.list <- list()
break.dt.list <- list()
for(sample.id in 1:n.samples){
  count <-
    c(rpois(data.per.segment, mu.low),
      rpois(data.per.segment, mu.high),
      rpois(data.per.segment, mu.low))
  chromEnd <- seq_along(count)
  chromStart <- chromEnd - 1L
  one.sample <-
    data.table(sample.id,
               chromStart,
               chromEnd,
               count)
  fit <- PeakSegDP(one.sample, maxPeaks=1L)
  seg.df <- subset(fit$segments, peaks == 1)
  model.dt.list[[sample.id]] <- 
    data.table(sample.id, model="PeakSeg",
               seg.df[, c("chromStart", "chromEnd", "mean")])
  break.dt.list[[sample.id]] <-
    data.table(sample.id, model="PeakSeg",
               subset(seg.df, chromStart != 0)[, "chromStart", drop=FALSE])
  signal.dt.list[[sample.id]] <- one.sample
}
signal.dt <- do.call(rbind, signal.dt.list)
truth.dt <- data.table(chromStart=c(100, 200), model="truth")

error.list <- list()
for(sample.size in 1:n.samples){
  these.signals <- signal.dt[sample.id <= sample.size, ]
  fit <- PeakSegJointSeveral(these.signals)
  converted <- ConvertModelList(fit)
  joint.segs <- subset(converted$segments, peaks == sample.size)
  joint.segs$model <- "PeakSegJoint"
  sub.model.list <- model.dt.list[1:sample.size]
  sub.model.list$joint <- joint.segs[, names(model.dt.list[[1]])]
  model.dt <- do.call(rbind, sub.model.list)
  joint.breaks <- subset(joint.segs, chromStart != 0)
  sub.break.list <- break.dt.list[1:sample.size]
  sub.break.list$joint <- joint.breaks[, names(break.dt.list[[1]])]
  break.dt <- do.call(rbind, sub.break.list)

  middle.segs <-
    model.dt[model=="PeakSeg" & chromStart != 0 & chromEnd != 300, ]
  mean.break <- function(x){
    floor(mean(x))
  }
  PeakSeg.breaks <- 
    middle.segs[, list(chromStart=mean.break(chromStart),
                       chromEnd=mean.break(chromEnd))]
  guess.mat <- 
    cbind(PeakSegJoint=unlist(joint.breaks[1, c("chromStart", "chromEnd")]),
          PeakSeg=unlist(PeakSeg.breaks))
  errors <- colSums(abs(guess.mat - truth.dt$chromStart))
  error.list[[sample.size]] <-
    data.table(sample.size, model=names(errors), errors)
  
  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ .)+
    geom_point(aes(chromEnd, count),
               color="grey50",
               data=signal.dt)+
    geom_vline(aes(xintercept=chromStart+0.5, color=model),
               data=truth.dt)+
    geom_segment(aes(chromStart+0.5, mean,
                     xend=chromEnd+0.5, yend=mean,
                     color=model, size=model),
                 data=model.dt)+
    geom_vline(aes(xintercept=chromStart+0.5, color=model, size=model),
               data=break.dt,
               linetype="dashed")+
    scale_size_manual(values=c(PeakSegJoint=0.5, PeakSeg=1))+
    scale_color_manual(values=color.code)
  
}

error <- do.call(rbind, error.list)

ggseed <- 
  ggplot()+
    scale_color_manual(values=color.code)+
    geom_line(aes(sample.size, errors, group=model, color=model),
              data=error)

pdf("figure-consistency.pdf")
print(ggseed)
dev.off()
