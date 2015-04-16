works_with_R("3.1.3",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@547ca81ce92c38b7c7f78cd071efc5afa96cb289")

data(H3K36me3.TDH.other.chunk1)
sample.ids <- sprintf("McGill%04d", c(22, 23, 19))
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         sample.id %in% sample.ids &
         43100000 < chromEnd & chromStart < 43205000)
some.counts$what <- "coverage"
some.regions <-
  subset(H3K36me3.TDH.other.chunk1$regions,
         sample.id %in% sample.ids &
         chromStart < 43205000)

library(ggplot2)
ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

sort.samples <- function(df){
  df$sample.id <- factor(df$sample.id, sample.ids)
  df
}
regions.plot <- 
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         c(0, floor(limits[2]))
                       })+
  xlab("position on chromosome (kilo bases = kb)")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                alpha=0.5,
                color="grey",
                data=sort.samples(some.regions))+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=sort.samples(some.counts),
            color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    paste("sample", seq_along(val))
  })

fit <- PeakSegJointHeuristic(some.counts)
converted <- ConvertModelList(fit)
zoom.peak.list <- with(converted, split(peaks, peaks$peaks))
zoom.peak.list[["0"]] <- Peaks()

## Compute PeakError on this sequence of models.
regions.by.sample <- split(some.regions, some.regions$sample.id, drop=TRUE)
counts.by.sample <- split(some.counts, some.counts$sample.id, drop=TRUE)
error.region.list <- list()
error.peaks.list <- list()
peak.label.list <- list()
for(peaks.str in names(zoom.peak.list)){
  several.samples <- zoom.peak.list[[peaks.str]]
  peaks.by.sample <- split(several.samples, several.samples$sample.id)
  peaks <- as.numeric(peaks.str)
  for(sample.id in names(regions.by.sample)){
    one.sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
      peaks.by.sample[[sample.id]]
    }else{
      Peaks()
    }
    max.count <- max(counts.by.sample[[sample.id]]$count)
    y <- (peaks+1) * -0.15 * max.count
    h <- -0.05 * max.count
    peak.label.list[[paste(sample.id, peaks)]] <-
      data.frame(sample.id, peaks, y)
    one.sample.regions <- regions.by.sample[[sample.id]]
    error <- PeakErrorChrom(one.sample.peaks, one.sample.regions)
    error.region.list[[paste(sample.id, peaks)]] <-
      data.frame(sample.id, peaks, error, y, h)
    if(nrow(one.sample.peaks)){
      error.peaks.list[[paste(sample.id, peaks)]] <-
        data.frame(one.sample.peaks, y, h)
    }
  }
}

error.peaks <- do.call(rbind, error.peaks.list)
error.regions <- do.call(rbind, error.region.list)
peak.labels <- do.call(rbind, peak.label.list)

err.plot <- 
regions.plot+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                linetype=status,
                ymin=y-h, ymax=y+h),
            fill=NA,
            color="black",
            data=error.regions)+
  scale_linetype_manual("error type",
                        limits=c("correct", 
                          "false negative",
                          "false positive"
                                 ),
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  geom_text(aes(43140,
                y, label=paste0(peaks, " peak",
                     ifelse(peaks==1, "", "s"),
                                " ")),
            size=2.5,
            hjust=0,
            data=peak.labels)+
  geom_segment(aes(chromStart/1e3, y,
                   xend=chromEnd/1e3, yend=y),
               size=1,
               color="deepskyblue",
               data=data.frame(what="peaks", error.peaks))

png("figure-PeakSegJoint.png", width=7, height=4, units="in", res=200)
print(err.plot)
dev.off()
