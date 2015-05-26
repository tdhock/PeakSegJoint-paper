works_with_R("3.2.0",
             data.table="1.9.4",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@a1ea491f49e9bdb347f1caadebe7b750de807ac4")

error.list <- list()
seg.list <- list()
for(set.name in dir("PeakSegJoint-chunks")){
  set.dir <- file.path("PeakSegJoint-chunks", set.name)
  problems.RData.vec <- Sys.glob(file.path(set.dir, "*", "problems.RData"))
  chunk.i <- 0
  for(problems.RData in problems.RData.vec){
    chunk.i <- chunk.i + 1
    chunk.dir <- dirname(problems.RData)
    chunk.id <- basename(chunk.dir)
    objs <- load(problems.RData)
    error.list[[problems.RData]] <- 
      data.table(set.name, chunk.i, chunk.id, res.error)
    seg.list[[problems.RData]] <-
      data.table(set.name, chunk.i, chunk.id,
                 res.error[errors==min(errors), ])
  }
}

segs <- do.call(rbind, seg.list)

errors.all <- do.call(rbind, error.list)
error.totals <-
  errors.all[, list(errors=sum(errors)), by=.(set.name, bases.per.problem)]

ggplot()+
  scale_x_log10()+
  geom_line(aes(bases.per.problem, errors),
            data=data.frame(error.totals, what="errors"))+
  geom_point(aes(bases.per.problem, chunk.i),
             data=data.frame(segs, what="chunk best"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ set.name, scales="free_y")
