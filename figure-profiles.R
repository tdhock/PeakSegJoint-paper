## -*- compile-command: "make HOCKING-PeakSegJoint-slides.pdf" -*-
works_with_R("3.2.0",
             ggplot2="1.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegDP@6520faad59fcf8c94ed4345d2b0d826b7f61faf9",
             "tdhock/PeakSegJoint@0390d43851e5427d3b2b45cb5a2b80b8007db67a")

data(H3K36me3.TDH.other.chunk1)

all.regions <- H3K36me3.TDH.other.chunk1$regions
to.rep <- all.regions$chromEnd==43223588
all.regions$chromStart[to.rep] <- 43205000
all.regions$chromEnd[to.rep] <- 43260000

make.tex <- function(prefix){

  tex.file <- paste0(prefix, ".tex")

  tex.list <- list()
  for(figure.i in seq_along(figure.list)){
    tit <- names(figure.list)[[figure.i]]
    png.file <- sprintf("%s-%d.png", prefix, figure.i)
    cat(sprintf("%d / %d figures\n", figure.i, length(figure.list)))
    tex.list[[figure.i]] <- sprintf("
\\begin{frame}
  \\frametitle{%s}
  \\includegraphics[width=\\textwidth]{%s}
\\end{frame}
", tit, png.file)
    gg <- figure.list[[figure.i]]
    png(png.file, width=7, height=5, units="in", res=200)
    print(gg)
    dev.off()
  }

  tex <- paste(tex.list, collapse="\n")
  cat(tex, file=tex.file)

}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

lims <- c(43100000, 43205000)
lims.df <- data.frame(chromStart=lims[1], chromEnd=lims[2])
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         lims[1] < chromEnd & chromStart < lims[2])
some.regions <- subset(all.regions,
                       chromStart < lims[2])

fit <- PeakSegJointSeveral(some.counts)
converted <- ConvertModelList(fit)

zoom.peak.list <-
  c(list("0"=Peaks()),
    with(converted, split(peaks, peaks$peaks)))
zoom.seg.list <- with(converted, split(segments, segments$peaks))

sample.order.list <- sapply(zoom.peak.list, function(df)paste(df$sample.id))
sample.order <- unique(unlist(sample.order.list))
sort.samples <- function(df){
  df$sample.id <- factor(df$sample.id, sample.order)
  df
}

figure.list <- list()

figure.list[["Peaks visually obvious in H3K36me3 data"]] <- 
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_blank(aes(linetype=status),
               data=data.frame(status=c("correct", "false positive",
                                 "false negative")))+
  xlab("position on chromosome (kilo bases = kb)")+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=sort.samples(H3K36me3.TDH.other.chunk1$counts),
            color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })

figure.list[["H3K36me3 data and visually determined labels"]] <- 
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_blank(aes(linetype=status),
               data=data.frame(status=c("correct", "false positive",
                                 "false negative")))+
  xlab("position on chromosome (kilo bases = kb)")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                alpha=0.5,
                color="grey",
                data=sort.samples(all.regions))+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=sort.samples(H3K36me3.TDH.other.chunk1$counts),
            color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })

one.id <- "McGill0022"
one.sample.counts <- 
  subset(H3K36me3.TDH.other.chunk1$counts,
         sample.id==one.id)
bases.per.bin <- with(one.sample.counts, {
  as.integer((max(chromEnd)-chromStart[1])/2000)
})
fit <-
  segmentBins(one.sample.counts,
              bin.chromStart=one.sample.counts$chromStart[1],
              bin.size=bases.per.bin)
segs.by.peaks <- with(fit$fit, split(segments, segments$peaks))
one.sample.regions <- 
  subset(all.regions,
         sample.id==one.id)
figure.list[["H3K36me3 data and labels (zoom to one sample)"]] <- one.gg <- 
  ggplot()+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          guide=guide_legend(direction="horizontal"),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_blank(aes(linetype=status),
               data=data.frame(status=c("correct", "false positive",
                                 "false negative")))+
    scale_y_continuous("aligned read coverage")+
    xlab("position on chromosome (kilo bases = kb)")+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  alpha=0.5,
                  color="grey",
                  data=one.sample.regions)+
    scale_fill_manual(values=ann.colors,
                    guide=guide_legend(direction="horizontal"))+
    geom_step(aes(chromStart/1e3, count),
              data=one.sample.counts,
              color="grey50")+
    theme_bw()+
    theme(legend.position="top")+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })

load("../PeakSeg-paper/dp.peaks.RData")
dp.model <- dp.peaks[["H3K36me3_TDH_other/1"]][[one.id]]

## Compute PeakError on this sequence of models.
regions.by.sample <- split(some.regions, some.regions$sample.id)
error.by.peaks <- list()
for(peaks.str in names(dp.model)){
  one.sample.peaks <- dp.model[[peaks.str]]
  error <- PeakErrorChrom(one.sample.peaks, one.sample.regions)
  peaks <- as.numeric(peaks.str)
  tit <-
    paste0("PeakSeg model with ",
           peaks, " peak",
           ifelse(peaks==1, "", "s"))
  zoom.segs <- segs.by.peaks[[peaks.str]]
  breaks <- subset(zoom.segs, min(chromStart) < chromStart)
  gg <- 
  one.gg +
    guides(linetype=guide_legend(order=2,
             override.aes=list(fill="white")))+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      linetype=status),
                  fill=NA,
                  color="black",
                  size=2,
                  data=error)+
    geom_vline(aes(xintercept=chromStart/1e3),
               data=breaks,
               color="green",
               linetype="dashed")+
    geom_segment(aes(chromStart/1e3, mean/bases.per.bin,
                     xend=chromEnd/1e3, yend=mean/bases.per.bin),
                 data=zoom.segs,
                 color="green",
                 size=1)
  if(0 < peaks){
    gg <- gg+
    geom_segment(aes(chromStart/1e3, -1,
                     xend=chromEnd/1e3, yend=-1),
                 data=one.sample.peaks,
                 size=2,
                 color="deepskyblue")
  }
  figure.list[[tit]] <- gg
}
print(error)

make.tex("figure-profiles-PeakSeg")

figure.list <- list()

figure.list[["H3K36me3 data and visually determined labels"]] <- 
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_blank(aes(linetype=status),
               data=data.frame(status=c("correct", "false positive",
                                 "false negative")))+
  xlab("position on chromosome (kilo bases = kb)")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                alpha=0.5,
                color="grey",
                data=sort.samples(all.regions))+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=sort.samples(H3K36me3.TDH.other.chunk1$counts),
            color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
  })

figure.list[["H3K36me3 data and labels (zoom to one peak)"]] <- regions.plot <- 
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_blank(aes(linetype=status),
               data=data.frame(status=c("correct", "false positive",
                                 "false negative")))+
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
    sub("McGill0", "", val)
  })

## Compute PeakError on this sequence of models.
regions.by.sample <- split(some.regions, some.regions$sample.id)
error.list <- list()
for(peaks.str in names(zoom.peak.list)){
  several.samples <- zoom.peak.list[[peaks.str]]
  peaks.by.sample <- split(several.samples, several.samples$sample.id)
  peaks <- as.numeric(peaks.str)
  error.by.sample <- list()
  for(sample.id in names(regions.by.sample)){
    one.sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
      peaks.by.sample[[sample.id]]
    }else{
      Peaks()
    }
    one.sample.regions <- regions.by.sample[[sample.id]]
    error <- PeakErrorChrom(one.sample.peaks, one.sample.regions)
    error.by.sample[[sample.id]] <- data.frame(sample.id, error)
  }
  peaks.error <- do.call(rbind, error.by.sample)
  fp <- sum(peaks.error$fp)
  fn <- sum(peaks.error$fn)
  error.list[[peaks.str]] <-
    data.frame(peaks, errors=fp+fn, regions=nrow(peaks.error))
  tit <-
    paste0("PeakSegJoint model with ",
           peaks, " peak",
           ifelse(peaks==1, "", "s"))
  zoom.segs <- zoom.seg.list[[peaks.str]]
  breaks <- subset(zoom.segs, min(chromStart) < chromStart)
  figure.list[[tit]] <- 
  regions.plot +
    guides(linetype=guide_legend(order=2,
             override.aes=list(fill="white")))+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      linetype=status),
                  fill=NA,
                  color="black",
                  size=2,
                  data=sort.samples(peaks.error))+
    geom_vline(aes(xintercept=chromStart/1e3),
               data=breaks,
               color="green",
               linetype="dashed")+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 data=sort.samples(zoom.segs),
                 color="green",
                 size=1)
}
error.sum <- do.call(rbind, error.list)
show.error <- data.frame(what="incorrect regions", error.sum)
show.loss <- data.frame(what="Poisson loss", converted$loss)
##figure.list[["Select model with minimal number of incorrect regions"]] <- 
ggplot()+
  scale_y_continuous("", breaks=function(limits){
    if(limits[1] < -50000){
      c(0, -50000, -100000)
    }else{
      c(0, 4, 8, 12)
    }
  })+
  geom_point(aes(peaks, loss),
             data=show.loss)+
  geom_line(aes(peaks, loss),
            data=show.loss)+
  geom_point(aes(peaks, errors),
             data=show.error)+
  geom_line(aes(peaks, errors),
            data=show.error)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")

make.tex("figure-profiles")
