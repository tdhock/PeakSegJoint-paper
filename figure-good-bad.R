works_with_R("3.2.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@86bee0a4620160e2d4b904e7819b5792280d51de",
             ggplot2="1.0")
             

db <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db"

getURL <- function(data.type){
  u <- sprintf("%s/H3K4me3_TDH_immune/5/%s.RData", db, data.type)
  dest.file <- basename(u)
  if(!file.exists(dest.file)){
    download.file(u, dest.file)
  }
  load(dest.file, envir=globalenv())
}
getURL("peaks/macs.default")
getURL("regions")
getURL("counts")
counts$count <- counts$coverage
sample.id.vec <- sprintf("McGill%04d", c(2, 4, 91, 322))
some.regions <- subset(regions, sample.id %in% sample.id.vec)
some.counts <- subset(counts, sample.id %in% sample.id.vec)
some.peaks <- subset(peaks[[1]], sample.id %in% sample.id.vec)
levs <- c(paste0("sample", 1:4), "bad", "good", "better")

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

some.counts$peaks <- ifelse(some.counts$chromStart < 118110000, 2, 4)
counts.by.peaks <- split(some.counts, some.counts$peaks)
joint.list <- list()
for(peaks.str in names(counts.by.peaks)){
  peak.counts <- counts.by.peaks[[peaks.str]]
  fit <- PeakSegJointSeveral(peak.counts)
  converted <- ConvertModelList(fit)
  all.peaks <- converted$peaks
  show.peaks <- subset(all.peaks, peaks == peaks.str)
  joint.list[[peaks.str]] <- show.peaks
}
joint <- do.call(rbind, joint.list)
joint$right <- factor("better", levs)
joint$left <- "0 errors"

counts.by.type <- split(some.counts, some.counts$cell.type, drop=TRUE)
each.type.list <- list()
for(cell.type in names(counts.by.type)){
  type.counts <- counts.by.type[[cell.type]]
  type.counts$peaks <- ifelse(type.counts$chromStart < 118110000, 2, 4)
  counts.by.peaks <- split(type.counts, type.counts$peaks)
  for(peaks.str in names(counts.by.peaks)){
    peak.counts <- counts.by.peaks[[peaks.str]]
    fit <- PeakSegJointSeveral(peak.counts)
    converted <- ConvertModelList(fit)
    all.peaks <- converted$peaks
    peaks.num <- if(cell.type=="bcell" && peaks.str==2){
      0
    }else{
      max(all.peaks$peaks)
    }
    if(peaks.num > 0){
      show.peaks <- subset(all.peaks, peaks == peaks.num)
      each.type.list[[paste(cell.type, peaks.str)]] <- show.peaks
    }
  }
}

some.regions$left <- some.regions$cell.type
McGill2sample <- function(x){
  i <- as.integer(factor(x, sample.id.vec))
  factor(paste0("sample", i), levs)
}
some.regions$right <- McGill2sample(some.regions$sample.id)
some.counts$left <- some.counts$cell.type
some.counts$right <- McGill2sample(some.counts$sample.id)
each.type <- do.call(rbind, each.type.list)
each.type$right <- factor("good", levs)
each.type$left <- "0 errors"
sample.id.y <- -seq_along(sample.id.vec)
names(sample.id.y) <- sample.id.vec
half.y <- 0.4
sampley <- function(df){
  df$y <- sample.id.y[paste(df$sample.id)]
  df$ymax <- df$y + half.y
  df$ymin <- df$y - half.y
  df
}
each.type$chromStart[1] <- 118094873
combine <- function(...){
  in.peaks <- list(...)
  peak.cols <- c("sample.id", "chromStart", "chromEnd", "right", "left")
  out.peaks <- lapply(in.peaks, "[", peak.cols)
  do.call(rbind, out.peaks)
}

some.peaks$left <- "5 errors"
some.peaks$right <- factor("bad", levs)
fake.peaks <- subset(some.peaks, sample.id != "McGill0002")
all.peaks <-
  sampley(combine(each.type, joint, fake.peaks))
peaks.by.right <- split(all.peaks, all.peaks$right, drop=TRUE)
regions.by.right <- list()
for(right in names(peaks.by.right)){
  p <- peaks.by.right[[right]]
  error.reg <- PeakErrorSamples(p, some.regions)
  regions.by.right[[right]] <-
    data.frame(right=factor(right, levs), left=p$left[1], error.reg)
}
error.regions <- sampley(do.call(rbind, regions.by.right))
sample.labels <- subset(error.regions, chromStart == min(chromStart))
sample.labels$sample <- McGill2sample(sample.labels$sample.id)

good.bad <- 
ggplot()+
  scale_y_continuous("aligned read coverage",
                     breaks=function(limits){
                       maybe <- floor(limits[2])
                       if(maybe==74)return(c(0, maybe))
                       if(maybe != -1)maybe else -100
                     })+
  geom_text(aes(chromStart/1e3, y,
                label=paste0(sample, " ")),
            data=sample.labels,
            size=3,
            hjust=1)+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                ymin=ymin, ymax=ymax,
                fill=annotation),
            data=error.regions,
            alpha=0.5,
            color="grey")+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                ymin=ymin, ymax=ymax,
                linetype=status),
            data=error.regions,
            size=0.75,
            fill=NA,
            color="black")+
  geom_segment(aes(chromStart/1e3, y,
                   xend=chromEnd/1e3, yend=y),
               data=all.peaks,
               size=3,
               color="deepskyblue")+
  scale_linetype_manual("error type",
                        limits=c("correct", 
                          "false negative",
                          "false positive"
                                 ),
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  xlab("position on chromosome (kilo bases = kb)")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                alpha=0.5,
                color="grey",
                data=some.regions)+
  scale_fill_manual("label", values=ann.colors)+
  geom_step(aes(chromStart/1e3, count),
            data=some.counts,
            color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(right + left ~ ., scales="free")

pdf("figure-good-bad.pdf", h=5, w=8)
print(good.bad)
dev.off()
