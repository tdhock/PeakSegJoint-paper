works_with_R("3.2.0",
             dplyr="0.4.0",
             ggplot2="1.0")

load("chunk.problems.RData")

all.error.list <- list()
for(set.name in names(weighted.error.list)){
  one.set <- weighted.error.list[[set.name]]
  for(chunk.name in names(one.set)){
    one.chunk <- one.set[[chunk.name]]
    for(res.str in names(one.chunk)){
      all.error.list[[paste(set.name, chunk.name, res.str)]] <-
        one.chunk[[res.str]]
    }
  }
}
all.error <- do.call(rbind, all.error.list)

with(all.error, table(set.name, total.weight))

chunk.err <- all.error %>%
  group_by(set.name, chunk.name, bases.per.problem) %>%
  summarise(weighted.error=sum(weighted.error),
            total.weight=sum(total.weight))
with(chunk.err, table(set.name, total.weight))

chunk.lims <- chunk.err %>%
  group_by(set.name, chunk.name) %>%
  summarise(min=min(total.weight),
            max=max(total.weight)) %>%
  group_by(set.name) %>%
  mutate(set.min=min(min))
stopifnot(with(chunk.lims, abs(min-max) < 1e-6))

small.chunks <- chunk.lims %>%
  filter(abs(min-set.min) < 1e-6)

set.res.stats <- chunk.err %>%
  ##filter(!chunk.name %in% small.chunks$chunk.name) %>%
  group_by(set.name, bases.per.problem) %>%
  summarise(weighted.error=sum(weighted.error),
            total.weight=sum(total.weight))
set.res.min <- set.res.stats %>%
  group_by(set.name) %>%
  filter(weighted.error==min(weighted.error))

set.res.lims <- set.res.stats %>%
  group_by(set.name) %>%
  summarise(min=min(total.weight),
            max=max(total.weight))
stopifnot(with(set.res.lims, min == max))

err.plot <- 
ggplot()+
  geom_point(aes(bases.per.problem, weighted.error),
             data=set.res.min)+
  geom_line(aes(bases.per.problem, weighted.error),
             data=set.res.stats)+
  theme_bw()+
  scale_x_log10()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(set.name ~ ., scales="free", labeller=function(var, val){
    gsub("_", "\n", val)
  })

pdf("figure-weighted-error.pdf", h=10, w=14)
print(err.plot)
dev.off()
