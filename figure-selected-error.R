works_with_R("3.2.0",
             dplyr="0.4.0",
             ggplot2="1.0")

load("selected.by.set.RData")

selected <- do.call(rbind, selected.by.set)

regions.range <- selected.error %>%
  group_by(set.name, split.i) %>%
  summarise(min=min(regions),
            max=max(regions))
stopifnot(with(regions.range, min==max))

ggplot()+
  geom_point(aes(bases.per.problem, errors),
             data=selected,
             pch=1)+
  geom_line(aes(bases.per.problem, errors, group=split.i),
            data=selected.error)+
  scale_x_log10()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_wrap("set.name")

pdf("figure-selected-error.pdf")
