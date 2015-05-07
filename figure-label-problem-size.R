works_with_R("3.2.0",
             ggplot2="1.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@ff5a7c58e297b54b328047f4e02285f0cb5d2838")

load("step1.RData")
load("cheating.error.RData")
load("chunk.problems.RData")

regions.file.list <- list()
regions.file.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/regions.RData")
for(regions.file in regions.file.vec){
  load(regions.file)
  regions.file.list[[regions.file]] <- regions
}

chosen.list <- list()
peak.size.list <- list()
guess.list <- list()
for(split.name in names(step1)){
  split.i <- sub(".* ", "", split.name)
  best.rows <- 
  data.frame(split.i,
             what=c("cheating", "step1"),
             bases.per.problem=as.numeric(c(step1[[split.name]]$res.str,
               best.res.list[[split.name]]$bases.per.problem)))
  set.name <- sub(" .*", "", split.name)
  region.files <- Sys.glob(sprintf("../chip-seq-paper/chunks/%s/*/regions.RData", set.name))
  bases.list <- list()
  for(region.file in region.files){
    regions <- regions.file.list[[region.file]]
    positive <- subset(regions, annotation %in% c("peakStart", "peakEnd"))
    bases.list[[region.file]] <- with(positive, chromEnd-chromStart)
  }
  bases.vec <- do.call(c, bases.list)
  med <- quantile(bases.vec, 0.25)
  guess.list[[set.name]] <- data.frame(set.name, min=med, max=med*100)
  peak.size.list[[set.name]] <- data.frame(set.name, bases.per.problem=bases.vec)
  data.by.chunk <- chunk.problems[[set.name]]
  set.name <- sub(" .*", "", split.name)
  chosen.list[[split.name]] <- data.frame(set.name, split.name, best.rows)
}
chosen.tall <- do.call(rbind, chosen.list)
peak.size <- do.call(rbind, peak.size.list)
guess <- do.call(rbind, guess.list)

psize <- 
ggplot()+
  ggtitle("25% quantile * 100 bases")+
  ylab("")+
  geom_tallrect(aes(xmin=min, xmax=max),
                alpha=0.5,
                data=guess)+
  scale_x_log10()+
  facet_grid(set.name ~ ., labeller=function(var, val){
    gsub("_", "\n", val)
  })+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_point(aes(bases.per.problem, what),
             data=data.frame(what="regions", peak.size))+
  geom_text(aes(bases.per.problem, what,
                label=split.i),
             data=chosen.tall)

pdf("figure-label-problem-size.pdf", h=10, w=14)
print(psize)
dev.off()
