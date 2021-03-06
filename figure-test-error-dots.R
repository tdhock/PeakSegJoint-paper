works_with_R("3.2.0",
             ggplot2="1.0",
             dplyr="0.4.0")

load("cheating.error.RData")
load("step1.error.RData")
load("step2.error.RData")
load("weighted.error.RData")
PeakSeg.results <- read.csv("PeakSeg-results.csv") %>%
  filter(learning != "unsupervised") %>%
  filter(!grepl("oracle.[13]", algorithm)) %>%
  mutate(algorithm=ifelse(algorithm=="oracle.41", "PeakSeg.41",
           paste(algorithm))) %>%
  mutate(algorithm=sub("[.][0-9]*$", "", algorithm))

step2.stats <- step2.error %>%
  mutate(algorithm="step2") %>%
  group_by(set.name, split.i, algorithm) %>%
  summarise(errors=sum(errors))
step2.stats$regions <- cheating.error$regions
step2.stats$percent <- with(step2.stats, errors/regions*100)
step2.stats$algo.type <- "PeakSegJoint"
step2.stats$learning <- "interval\nregression"

step2.all.stats <- step2.error.all %>%
  mutate(algorithm="PeakSegJoint") %>%
  group_by(set.name, split.i, algorithm) %>%
  summarise(errors=sum(fn+fp),
            regions=n()) %>%
  mutate(percent=errors/regions*100,
         algo.type="PeakSegJoint",
         learning="interval\nregression")

weighted.all.stats <- weighted.error.all %>%
  mutate(algorithm="weighted") %>%
  group_by(set.name, split.i, algorithm) %>%
  summarise(errors=sum(fn+fp),
            regions=n()) %>%
  mutate(percent=errors/regions*100,
         algo.type="PeakSegJoint",
         learning="interval\nregression")

step1.stats <- step1.error %>%
  mutate(algorithm="step1") %>%
  group_by(set.name, split.i, algorithm) %>%
  summarise(errors=sum(fp+fn),
            regions=n()) %>%
  mutate(percent=errors/regions*100,
         algo.type="PeakSegJoint",
         learning="interval\nregression")

step1.best.stats <- data.frame(best.for.train.res) %>%
  mutate(algorithm="best.for\nselected\nresolution",
         percent=errors/regions*100,
         algo.type="PeakSegJoint",
         learning="cheating")

cheating.stats <- data.frame(cheating.error) %>%
  mutate(algorithm="test.res",
         percent=errors/regions*100,
         algo.type="PeakSegJoint",
         learning="cheating")

common.names <- names(PeakSeg.results)
all.stats <-
  rbind(
    PeakSeg.results,
    ##cheating.stats[, common.names], #comment to hide cheaters.
    ##step1.best.stats[, common.names], #comment to hide cheaters.
    ##step1.stats,
    ##step2.stats,
    ##weighted.all.stats,
    step2.all.stats)

show.stats <- all.stats %>%
  filter(!grepl("AIC/BIC", algorithm),
         !grepl("NTNU", set.name))

region.range <- show.stats %>%
  group_by(set.name, split.i) %>%
  summarise(min=min(regions),
            max=max(regions)) %>%
  mutate(diff=max-min)
##stopifnot(with(region.range, min == max))
data.frame(region.range)

show.means <- show.stats %>%
  group_by(set.name, algorithm, learning, algo.type) %>%
  summarise(percent=mean(percent))

## Some quantitative results for the epigenomics abstract.
med.stats <- show.stats %>%
  group_by(set.name, algorithm, learning, algo.type) %>%
  summarise(median=median(percent),
            mean=mean(percent),
            sd=sd(percent),
            quartile25=quantile(percent, 0.25),
            quartile75=quantile(percent, 0.75)) %>%
  filter(grepl("TDH_immune", set.name),
         algorithm != "PeakSeg") %>%
  group_by()

med.stats %>%
  select(set.name, algorithm, quartile25, median, quartile75) 
med.stats %>%
  select(set.name, algorithm, mean, sd) 

show.vlines <- show.means %>%
  group_by(set.name) %>%
  filter(seq_along(percent) == which.min(percent)) %>%
  select(-algo.type)
show.vlines <- show.means %>%
  group_by(set.name) %>%
  filter(algorithm=="PeakSegJoint") %>%
  select(-algo.type)

algo.colors <-
  c("cheating"="grey",
    "interval\nregression"="#D95F02",
    "grid\nsearch"="#1B9E77",
    "unsupervised"="#7570B3")

dots <-  #with 1 set of facets.
ggplot()+
  geom_vline(aes(xintercept=percent),
             data=show.vlines)+
  geom_point(aes(percent, algorithm, color=learning),
             data=show.means,
             alpha=0.2,
             size=3)+
  geom_point(aes(percent, algorithm, color=learning),
             data=show.stats, pch=1)+
  facet_grid(. ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  }, scales="free_y", space="free_y")+
  scale_y_discrete("model")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_color_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_fill_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_x_continuous("percent incorrect peak region labels (test error)",
                     breaks=seq(0, 100, by=25))

pdf("figure-test-error-dots.pdf", h=2.5, w=7)
print(dots)
dev.off()
