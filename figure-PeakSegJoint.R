works_with_R("3.2.0",
             "tdhock/PeakSegDP@fe06a5b91d68c5d1ec471cb15c3ec3935dc2624d",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@86bee0a4620160e2d4b904e7819b5792280d51de")

data(H3K36me3.TDH.other.chunk1)
sample.ids <- sprintf("McGill%04d", c(23, 22, 37))
sample.y <- -seq_along(sample.ids)*100
sample.h <- 40
names(sample.y) <- sample.ids
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         sample.id %in% sample.ids &
         43100000 < chromEnd & chromStart < 43205000)
some.counts$what <- "coverage"
some.regions <-
  subset(H3K36me3.TDH.other.chunk1$regions,
         sample.id %in% sample.ids &
         chromStart < 43205000)
is.first <- some.regions$chromEnd == min(some.regions$chromEnd)
some.regions[is.first, "chromEnd"] <- 43168000

library(ggplot2)
ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

plevs <- c("sample 1", "sample 2", "sample 3",
           "0 peaks", "1 peak", "2 peaks", "3 peaks")

id2sample <- function(x){
  sample.id <- factor(paste(x), sample.ids)
  panel.str <- paste("sample", as.integer(sample.id))
  factor(panel.str, plevs)
}
  
sort.samples <- function(df){
  df$sample.id <- factor(df$sample.id, sample.ids)
  panel.str <- paste("sample", as.integer(df$sample.id))
  df$panel <- factor(panel.str, plevs)
  df
}
regions.plot <- 
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         hi <- floor(limits[2])
                         if(hi==10){
                           c(0, hi)
                         }else if(0 < hi){
                           hi
                         }else{
                           -Inf
                         }
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
  facet_grid(panel ~ ., scales="free")

profile.list <- ProfileList(some.counts)
step1 <- PeakSegJointHeuristicStep1(profile.list)
unconstrained.list <- list()
for(sample.id in names(profile.list)){
  pro <- profile.list[sample.id]
  fit <- PeakSegJointSeveral(pro)
  converted <- ConvertModelList(fit)
  unconstrained.list[[sample.id]] <-
    data.frame(converted$peaks, diff=diff(converted$loss$loss))
}
unconstrained.peaks.list <- list()
for(peaks in seq_along(unconstrained.list)){
  id.vec <- c("McGill0022", "McGill0023", "McGill0037")[1:peaks]
  peak.df <- do.call(rbind, unconstrained.list[id.vec])
  peak.df$peaks <- peaks
  unconstrained.peaks.list[[paste(peaks)]] <- peak.df
}
fit <- PeakSegJointHeuristic(profile.list)
converted <- ConvertModelList(fit)
converted$peaks$diff <- NA
peaks.by.model <- 
  list(
    ##Unconstrained=unconstrained.peaks.list,
    PeakSegJoint=with(converted, split(peaks, peaks$peaks)))
error.region.list <- list()
error.peaks.list <- list()
peak.label.list <- list()
for(model in names(peaks.by.model)){
  zoom.peak.list <- peaks.by.model[[model]]
  zoom.peak.list[["0"]] <- Peaks()

  ## Compute PeakError on this sequence of models.
  regions.by.sample <- split(some.regions, some.regions$sample.id, drop=TRUE)
  counts.by.sample <- split(some.counts, some.counts$sample.id, drop=TRUE)
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
      peak.label.list[[paste(sample.id, peaks, model)]] <-
        data.frame(sample.id, peaks, model, y)
      one.sample.regions <- regions.by.sample[[sample.id]]
      error <- PeakErrorChrom(one.sample.peaks, one.sample.regions)
      panel.str <- sprintf("%d peak%s", peaks, ifelse(peaks==1, "", "s"))
      panel <- factor(panel.str, plevs)
      error.region.list[[paste(sample.id, peaks, model)]] <-
        data.frame(sample.id, peaks, error, y, h,
                   sample=id2sample(sample.id),
                   panel,
                   model,
                   sample.y=sample.y[[sample.id]])
      if(nrow(one.sample.peaks)){
        error.peaks.list[[paste(sample.id, peaks, model)]] <-
          data.frame(one.sample.peaks, y, h,
                     sample=id2sample(sample.id),
                     panel,
                     model,
                     sample.y=sample.y[[sample.id]])
      }
    }
  }
}

error.peaks <- do.call(rbind, error.peaks.list)
error.regions <- do.call(rbind, error.region.list)
peak.labels <- do.call(rbind, peak.label.list)

first.regions <- subset(error.regions, chromStart==min(chromStart))

err.plot <- 
regions.plot+
  geom_text(aes(chromStart/1e3, sample.y,
                label=paste0(sample, " ")),
            hjust=1,
            size=3.5,
            data=first.regions)+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                fill=annotation,
                ymin=sample.y-sample.h, ymax=sample.y+sample.h),
            color="grey",
            alpha=0.5,
            data=error.regions)+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                linetype=status,
                ymin=sample.y-sample.h, ymax=sample.y+sample.h),
            fill=NA,
            size=1,
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
  ## geom_text(aes(43140,
  ##               y, label=paste0(peaks, " peak",
  ##                    ifelse(peaks==1, "", "s"),
  ##                               " ")),
  ##           size=2.5,
  ##           hjust=0,
  ##           data=peak.labels)+
  scale_size_manual(values=c(Unconstrained=3, PeakSegJoint=1.5))+
  scale_color_manual(values=c(Unconstrained="royalblue",
                       PeakSegJoint="deepskyblue"))+
  geom_segment(aes(chromStart/1e3, sample.y,
                   color=model, size=model,
                   xend=chromEnd/1e3, yend=sample.y),
               data=data.frame(what="peaks", error.peaks))
if(length(peaks.by.model) < 2){
  err.plot <- err.plot + guides(color="none", size="none")
}

png("figure-PeakSegJoint.png", width=8, height=5, units="in", res=200)
print(err.plot)
dev.off()
