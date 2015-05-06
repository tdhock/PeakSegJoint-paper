works_with_R("3.2.0",
             ggplot2="1.0",
             dplyr="0.4.0")

load("step1.error.RData")
PeakSeg.results <- read.csv("PeakSeg-results.csv")

step1.stats <- step1.error %>%
  mutate(algorithm="PeakSegJoint") %>%
  group_by(set.name, split.i, algorithm) %>%
  summarise(errors=sum(fp+fn),
            regions=n()) %>%
  mutate(percent=errors/regions*100,
         algo.type="PeakSegJoint",
         learning="interval\nregression")

all.stats <- rbind(PeakSeg.results, step1.stats)

region.range <- all.stats %>%
  group_by(set.name, split.i) %>%
  summarise(min=min(regions),
            max=max(regions))
stopifnot(with(region.range, min == max))

show.stats <- all.stats %>%
  filter(!grepl("AIC/BIC", algorithm),
         !grepl("NTNU", set.name))

show.means <- show.stats %>%
  group_by(set.name, algorithm, learning, algo.type) %>%
  summarise(percent=mean(percent))

show.vlines <- show.means %>%
  group_by(set.name) %>%
  filter(seq_along(percent) == which.min(percent)) %>%
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
  facet_grid(algo.type ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  }, scales="free_y", space="free_y")+
  scale_y_discrete("model . parameters learned")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_color_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_fill_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_x_continuous("percent incorrect peak region labels (test error)",
                     breaks=seq(0, 100, by=20))

pdf("figure-test-error-dots.pdf", h=4.2, w=8)
print(dots)
dev.off()
