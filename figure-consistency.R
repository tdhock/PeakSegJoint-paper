works_with_R("3.2.2",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@96cefe9e2208779266ed73e55841c46d7778a478")

color.code <-
  c(truth="#1B9E77", #teal
    PeakSeg="#D95F02", #orange
    PeakSegJoint="#7570B3", #violet
    "#E7298A", #pink
    "#66A61E", #green
    "#E6AB02", #tan
    "#A6761D", #brown
    "#666666") #grey

load("PeakConsistency.RData")

viz <- list(
  increase=ggplot()+
    ggtitle("Select true increase in peak mean")+
    ylab("<- PeakSeg better ------- PeakSegJoint better ->")+
    make_tallrect(PeakConsistency$increase, "increase")+
    geom_line(aes(increase, mean.diff), data=PeakConsistency$increase),
  
  errors=ggplot()+
    ggtitle("Select random seed and sample size")+
    scale_x_continuous("sample size", breaks=c(1, 5, 10, 15))+
    ylab("distance from true peaks to estimated peaks")+
    scale_color_manual(values=color.code)+
    make_tallrect(PeakConsistency$error, "sample.size")+
    geom_line(aes(sample.size, errors,
                  clickSelects=seed,
                  showSelected=increase,
                  group=interaction(model, seed),
                  color=model),
              size=5,
              alpha=0.7,
              data=PeakConsistency$error),
  title="Consistency of peak recovery in simulated Poisson data",
  ##duration=list(seed=1000, increase=1000, sample.size=1000),
  signals=ggplot()+
    ggtitle("True and estimated peak positions")+
    theme_bw()+
    theme_animint(width=1000, height=1000)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      paste("sample", val)
    })+
    geom_point(aes(chromEnd, count,
                   key=chromEnd,
                   showSelected3=increase,
                   showSelected=seed),
               color="grey50",
               data=PeakConsistency$signal)+
    geom_vline(aes(xintercept=chromStart+0.5, color=model,
                   showSelected3=increase,
                   showSelected=seed),
               show_guide=TRUE,
               linetype="dashed",
               data=PeakConsistency$truth)+
    guides(size="none")+
    geom_segment(aes(chromStart+0.5, mean,
                     xend=chromEnd+0.5, yend=mean,
                     key=paste(model, seg.i),
                     showSelected3=increase,
                     showSelected=seed, showSelected2=sample.size,
                     color=model, size=model),
                 data=PeakConsistency$model)+
    geom_vline(aes(xintercept=chromStart+0.5,
                   key=paste(model, break.i),
                   showSelected=seed, showSelected2=sample.size,
                   showSelected3=increase,
                   color=model, size=model),
               show_guide=TRUE,
               linetype="dashed",
               data=PeakConsistency$guess)+
    scale_size_manual(values=c(PeakSegJoint=1, PeakSeg=2))+
    scale_color_manual(values=color.code))

viz$errors+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ increase)

unlink("figure-consistency", recursive=TRUE)
animint2dir(viz, "figure-consistency")

pdf("figure-consistency.pdf")
print(gg)
dev.off()
