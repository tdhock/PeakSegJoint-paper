works_with_R("3.2.0",
             ggplot2="1.0",
             dplyr="0.4.0")

load("cheating.error.RData")
load("step1.RData")

resolutions <- do.call(rbind, best.res.list)
resolutions$cheating <- resolutions$bases.per.problem
resolutions$step1 <- NA
for(split.name in names(step1)){
  split.data <- step1[[split.name]]
  resolutions[split.name, "step1"] <- as.numeric(split.data$res.str)
}
resolutions$experiment <- sub("_.*", "", resolutions$set.name)
resolutions$experiment.type <-
  ifelse(grepl("^H", resolutions$experiment),
         resolutions$experiment,
         "TF")

abline.df <- data.frame(slope=1, intercept=0)

resPlot <- 
ggplot()+
  geom_text(aes(step1, cheating,
                color=experiment.type,
                label=split.i),
             data=resolutions)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  coord_equal()+
  geom_abline(aes(slope=slope, intercept=intercept),
              color="grey",
              data=abline.df)+
  scale_x_log10("resolution selected in training set")+
  scale_y_log10("best resolution for test set")+
  ggtitle("resolution in bases per problem (6 train/test splits plotted)")+
  facet_wrap("set.name")

pdf("figure-scatter-cheating-step1.pdf", h=5)
print(resPlot)
dev.off()
