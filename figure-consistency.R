works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@f7d566dc6afe1e9fdbd636cc78199258c541b7f9",
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
increase <- 1
n.samples <- 10

viz.error.list <- list()
viz.model.list <- list()
viz.guess.list <- list()
viz.truth.list <- list()
viz.signal.list <- list()
for(increase in 1:10){
  print(increase)
  mu.high <- mu.low + increase
  for(seed in 1:4){
    set.seed(seed)
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
    viz.signal.list[[paste(increase, seed)]] <-
      data.table(increase, seed, signal.dt)
    truth.dt <- data.table(chromStart=c(100, 200), model="truth")
    viz.truth.list[[paste(increase, seed)]] <-
      data.table(increase, seed, truth.dt)

    error.list <- list()
    for(sample.size in 1:n.samples){
      these.signals <- signal.dt[sample.id <= sample.size, ]
      fit <- PeakSegJointSeveral(these.signals)
      converted <- ConvertModelList(fit)
      joint.segs <- subset(converted$segments, peaks == max(peaks))
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
      guess.by.model <- list()
      for(model in colnames(guess.mat)){
        guess.by.model[[model]] <-
          data.table(model, chromStart=guess.mat[, model],
                     sample.id=rep(1:sample.size, each=2))
      }
      guess.dt <- do.call(rbind, guess.by.model)
      errors <- colSums(abs(guess.mat - truth.dt$chromStart))
      error.list[[sample.size]] <-
        data.table(increase, sample.size, model=names(errors), errors)
      
      ggplot()+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.id ~ ., labeller=function(var, val){
          paste("sample", val)
        })+
        geom_point(aes(chromEnd, count),
                   color="grey50",
                   data=signal.dt)+
        geom_vline(aes(xintercept=chromStart+0.5, color=model),
                   show_guide=TRUE,
                   linetype="dashed",
                   data=truth.dt)+
        guides(size="none")+
        geom_segment(aes(chromStart+0.5, mean,
                         xend=chromEnd+0.5, yend=mean,
                         color=model, size=model),
                     data=model.dt)+
        geom_vline(aes(xintercept=chromStart+0.5, color=model, size=model),
                   show_guide=TRUE,
                   linetype="dashed",
                   data=guess.dt)+
        scale_size_manual(values=c(PeakSegJoint=0.5, PeakSeg=1))+
        scale_color_manual(values=color.code)

      viz.model.list[[paste(increase, sample.size, seed)]] <-
        data.table(increase, sample.size, seed, model.dt)
      viz.guess.list[[paste(increase, sample.size, seed)]] <-
        data.table(increase, sample.size, seed, guess.dt)
    }

    seed.error <- do.call(rbind, error.list)
    ggseed <- 
      ggplot()+
        ylab("distance from true peaks to estimated peaks")+
        scale_color_manual(values=color.code)+
        geom_line(aes(sample.size, errors, group=model, color=model),
                  data=seed.error)

    viz.error.list[[paste(increase, seed)]] <-
      data.table(increase, seed, seed.error)
  }
}

viz.model <- do.call(rbind, viz.model.list)
viz.truth <- do.call(rbind, viz.truth.list)
viz.signal <- do.call(rbind, viz.signal.list)
viz.guess <- do.call(rbind, viz.guess.list)
viz.error <- do.call(rbind, viz.error.list)
setkey(viz.error, model)
(gg <- ggplot()+
  ylab("distance from true peaks to estimated peaks")+
  scale_color_manual(values=color.code)+
  geom_line(aes(sample.size, errors,
                group=interaction(model, seed), color=model),
            data=viz.error)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ increase))

viz.diff <- viz.error["PeakSegJoint"]
viz.diff$PeakSeg <- viz.error["PeakSeg"]$errors
viz.diff[, diff := PeakSeg - errors ]
ggplot()+
  ylab("<- PeakSeg better ------- PeakSegJoint better ->")+
  geom_line(aes(sample.size, diff, group=seed),
            data=viz.diff)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ increase)

increase.diffs <- viz.diff[, list(mean.diff=mean(diff)), by=increase]
ggplot()+
  make_tallrect(increase.diffs, "increase")+
  geom_line(aes(increase, mean.diff), data=increase.diffs)

## PeakConsistency <-
##   list(model=data.frame(viz.model),
##        truth=data.frame(viz.truth),
##        signal=data.frame(viz.signal),
##        guess=data.frame(viz.guess),
##        error=data.frame(all.seeds.error))
## save(PeakConsistency, file="~/R/animint/data/PeakConsistency.RData")
## prompt(PeakConsistency, file="~/R/animint/man/PeakConsistency.Rd")

viz <-
  list(increase=ggplot()+
         make_tallrect(increase.diffs, "increase")+
         geom_line(aes(increase, mean.diff), data=increase.diffs),
       
         errors=ggplot()+
         ylab("distance from true peaks to estimated peaks")+
         scale_color_manual(values=color.code)+
         make_tallrect(viz.error, "sample.size")+
         geom_line(aes(sample.size, errors,
                       clickSelects=seed,
                       showSelected=increase,
                       group=interaction(model, seed),
                       color=model),
                   size=5,
                   alpha=0.7,
                   data=viz.error),

       signals=ggplot()+
         theme_bw()+
         theme_animint(width=1000, height=800)+
         theme(panel.margin=grid::unit(0, "cm"))+
         facet_grid(sample.id ~ ., labeller=function(var, val){
           paste("sample", val)
         })+
         geom_point(aes(chromEnd, count,
                        showSelected3=increase,
                        showSelected=seed),
                    color="grey50",
                    data=viz.signal)+
         geom_vline(aes(xintercept=chromStart+0.5, color=model,
                        showSelected3=increase,
                        showSelected=seed),
                    show_guide=TRUE,
                    linetype="dashed",
                    data=viz.truth)+
         guides(size="none")+
         geom_segment(aes(chromStart+0.5, mean,
                          xend=chromEnd+0.5, yend=mean,
                          showSelected3=increase,
                          showSelected=seed, showSelected2=sample.size,
                          color=model, size=model),
                      data=viz.model)+
         geom_vline(aes(xintercept=chromStart+0.5,
                        showSelected=seed, showSelected2=sample.size,
                        showSelected3=increase,
                        color=model, size=model),
                    show_guide=TRUE,
                    linetype="dashed",
                    data=viz.guess)+
         scale_size_manual(values=c(PeakSegJoint=0.5, PeakSeg=1))+
         scale_color_manual(values=color.code))

unlink("figure-consistency", recursive=TRUE)
animint2dir(viz, "figure-consistency")

pdf("figure-consistency.pdf")
print(gg)
dev.off()
