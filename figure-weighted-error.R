works_with_R("3.2.0",
             dplyr="0.4.0",
             ggplot2="1.0")

error.list <- list()
for(set.name in dir("chunk.problems")){
  set.path <- file.path("chunk.problems", set.name)
  for(chunk.id in dir(set.path)){
    chunk.path <- file.path(set.path, chunk.id)
    RData.files <- dir(chunk.path)
    for(RData.file in RData.files){
      RData.path <- file.path(chunk.path, RData.file)
      load(RData.path)
      error.list[[RData.path]] <- error.row
    }
  }
}

all.error <- do.call(rbind, error.list)
set.res.stats <- all.error %>%
  group_by(set.name, bases.per.problem) %>%
  summarise(weighted.error=sum(weighted.error),
            total.weight=sum(total.weight))

ggplot()+
  geom_line(aes(bases.per.problem, weighted.error),
             data=set.res.stats)+
  theme_bw()+
  scale_x_log10()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(set.name ~ ., scales="free")
