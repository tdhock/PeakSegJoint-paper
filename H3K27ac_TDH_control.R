works_with_R("3.1.3",
             ggplot2="1.0",
             "tdhock/PeakSegJoint@e6bc386c555e203cc80343814939d51785c03af1",
             data.table="1.9.4")

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

g.pos.pattern <-
  paste0("(?<chrom>chr.+?)",
         ":",
         "(?<chromStart>[0-9 ,]+)",
         "-",
         "(?<chromEnd>[0-9 ,]+)",
         " ",
         "(?<annotation>[a-zA-Z]+)",
         "(?<types>.*)")

## Parse the first occurance of pattern from each of several strings
## using (named) capturing regular expressions, returning a matrix
## (with column names).
str_match_perl <- function(string,pattern){
  stopifnot(is.character(string))
  stopifnot(is.character(pattern))
  stopifnot(length(pattern)==1)
  parsed <- regexpr(pattern,string,perl=TRUE)
  captured.text <- substr(string,parsed,parsed+attr(parsed,"match.length")-1)
  captured.text[captured.text==""] <- NA
  captured.groups <- do.call(rbind,lapply(seq_along(string),function(i){
    st <- attr(parsed,"capture.start")[i,]
    if(is.na(parsed[i]) || parsed[i]==-1)return(rep(NA,length(st)))
    substring(string[i],st,st+attr(parsed,"capture.length")[i,]-1)
  }))
  result <- cbind(captured.text,captured.groups)
  colnames(result) <- c("",attr(parsed,"capture.names"))
  result
}

labels.file <- "H3K27ac_TDH_control.txt"

cat("Reading ", labels.file, "\n", sep="")
labels.lines <- readLines(labels.file)
is.blank <- labels.lines == ""
chunk.id <- cumsum(is.blank)
label.df <- data.frame(chunk.id, line=labels.lines)[!is.blank, ]
cat(length(unique(label.df$chunk.id)), " chunks, ",
    nrow(label.df), " label lines\n", sep="")

## Error checking.
raw.vec <- paste(label.df$line)
line.vec <- gsub(",", "", raw.vec)
match.mat <- str_match_perl(line.vec, g.pos.pattern)
stopifnot(!is.na(match.mat[,1]))
not.recognized <-
  ! match.mat[, "annotation"] %in% c("peakStart", "peakEnd", "peaks", "noPeaks")
if(any(not.recognized)){
  print(raw.vec[not.recognized])
  print(match.mat[not.recognized, drop=FALSE])
  stop("unrecognized annotation")
}
match.df <-
  data.frame(chrom=match.mat[, "chrom"],
             chromStart=as.integer(match.mat[, "chromStart"]),
             chromEnd=as.integer(match.mat[, "chromEnd"]),
             annotation=match.mat[, "annotation"],
             types=match.mat[, "types"],
             chunk.id=label.df$chunk.id,
             stringsAsFactors=FALSE)
match.by.chrom <- split(match.df, match.df$chrom)
for(chrom in names(match.by.chrom)){
  chrom.df <- match.by.chrom[[chrom]]
  sorted <- chrom.df[with(chrom.df, order(chromStart, chromEnd)), ]
  same.as.next <- diff(sorted$chromStart) <= 0
  if(any(same.as.next)){
    bad.i <- which(same.as.next)
    print(sorted[c(bad.i, bad.i+1), ])
    stop("chromStart not increasing")
  }
  if(any(with(sorted, chromStart >= chromEnd))){
    print(sorted)
    stop("chromStart >= chromEnd")
  }
  overlaps.next <-
    with(sorted, chromStart[-1] < chromEnd[-length(chromEnd)])
  if(any(overlaps.next)){
    print(data.frame(sorted, overlaps.next=c(overlaps.next, FALSE)))
    stop("overlapping regions")
  }
}

## determine total set of cell types with positive=Peak annotations.
stripped <- gsub(" *$", "", gsub("^ *", "", match.df$types))
type.list <- strsplit(stripped, split=" ")
names(type.list) <- rownames(match.df)
sample.types <- unique(unlist(type.list))
cat("sample types with peak annotations: ",
    paste(sample.types, collapse=", "),
    "\n",
    sep="")

match.by.chunk <- split(match.df, match.df$chunk.id)

sample.list <-
  list(skeletalMuscle=c(McGill0012="http://epigenomesportal.ca/public_data/donor/McGill0012/MS001201/H3K27ac/MS001201.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0014="http://epigenomesportal.ca/public_data/donor/McGill0014/MS001401/H3K27ac/MS001401.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0013="http://epigenomesportal.ca/public_data/donor/McGill0013/MS001301/H3K27ac/MS001301.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0015="http://epigenomesportal.ca/public_data/donor/McGill0015/MS001501/H3K27ac/MS001501.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0015="http://epigenomesportal.ca/public_data/donor/McGill0016/MS001601/H3K27ac/MS001601.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0019="http://epigenomesportal.ca/public_data/donor/McGill0019/MS001901/H3K27ac/MS001901.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0036="http://epigenomesportal.ca/public_data/donor/McGill0036/MS003601/H3K27ac/MS003601.skeletalMuscle.H3K27ac.signal.bigWig",
McGill0037="http://epigenomesportal.ca/public_data/donor/McGill0037/MS003701/H3K27ac/MS003701.skeletalMuscle.H3K27ac.signal.bigWig"),
       skeletalMuscleInput=c(McGill0012="http://epigenomesportal.ca/public_data/donor/McGill0012/MS001201/Input/MS001201.skeletalMuscle.Input.signal.bigWig",
McGill0013="http://epigenomesportal.ca/public_data/donor/McGill0013/MS001301/Input/MS001301.skeletalMuscle.Input.signal.bigWig",
McGill0014="http://epigenomesportal.ca/public_data/donor/McGill0014/MS001401/Input/MS001401.skeletalMuscle.Input.signal.bigWig",
McGill0015="http://epigenomesportal.ca/public_data/donor/McGill0015/MS001501/Input/MS001501.skeletalMuscle.Input.signal.bigWig",
McGill0016="http://epigenomesportal.ca/public_data/donor/McGill0016/MS001601/Input/MS001601.skeletalMuscle.Input.signal.bigWig",
McGill0019="http://epigenomesportal.ca/public_data/donor/McGill0019/MS001901/Input/MS001901.skeletalMuscle.Input.signal.bigWig",
McGill0036="http://epigenomesportal.ca/public_data/donor/McGill0036/MS003601/Input/MS003601.skeletalMuscle.Input.signal.bigWig",
McGill0037="http://epigenomesportal.ca/public_data/donor/McGill0037/MS003701/Input/MS003701.skeletalMuscle.Input.signal.bigWig"),
       thyroid=c(CEMT_44="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34711.IX2581.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_5_TATAAT.q5.F1028.PET.ucsc.bw",
CEMT_43="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34392.IX2580.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_4_AAGACT.q5.F1028.PET.ucsc.bw",
CEMT_45="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34718.IX2582.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_6_ATACGG.q5.F1028.PET.ucsc.bw",
CEMT_40="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34086.IX2577.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_1_CAAAAG.q5.F1028.PET.ucsc.bw",
CEMT_42="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34385.IX2579.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_3_CTATAC.q5.F1028.PET.ucsc.bw",
CEMT_41="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34093.IX2578.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_2_CGATGT.q5.F1028.PET.ucsc.bw"),
       thyroidInput=c(CEMT_44="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34712.IX2581.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_5_TGCTGG.q5.F1028.PET.ucsc.bw",
CEMT_43="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34393.IX2580.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_4_ACGATA.q5.F1028.PET.ucsc.bw",
CEMT_45="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34719.IX2582.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_6_CACGAT.q5.F1028.PET.ucsc.bw",
CEMT_40="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34087.IX2577.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_1_CCAACA.q5.F1028.PET.ucsc.bw",
CEMT_42="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34386.IX2579.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_3_GCAAGG.q5.F1028.PET.ucsc.bw",
CEMT_41="http://www.bcgsc.ca/downloads/edcc/data/CEMT/ChIP/A34094.IX2578.75nt.hg19a.bwa-0.5.7.C3LU8ACXX_2_TAGCTT.q5.F1028.PET.ucsc.bw"))

chunk.list <- list()
for(chunk.id in names(match.by.chunk)){
  ## Check that all regions are on the same chrom.
  chunk.df <- match.by.chunk[[chunk.id]]
  chunkChrom <- paste(chunk.df$chrom[1])
  if(any(chunk.df$chrom != chunkChrom)){
    print(chunk.df)
    stop("each chunk must span only 1 chrom")
  }
  chunk.name <- paste0("H3K27ac_TDH_control/", chunk.id)
  region.chromEnd <- max(chunk.df$chromEnd)
  region.chromStart <- min(chunk.df$chromStart)
  region.bases <- region.chromEnd  - region.chromStart
  expand <- region.bases/10
  chunkStart <- as.integer(region.chromStart - expand)
  chunkEnd <- as.integer(region.chromEnd + expand)
  chunk.sample.list <- list()
  for(sample.type in names(sample.list)){
    sample.vec <- sample.list[[sample.type]]
    for(sample.i in seq_along(sample.vec)){
      sample.id <- names(sample.vec)[[sample.i]]
      cat(sprintf("%4d / %4d %s samples %s\n",
                  sample.i, length(sample.vec), sample.type, sample.id))
      u <- sample.vec[[sample.id]]
      bg.file <-
        sprintf("H3K27ac_TDH_control/%s/%s/%s.bedGraph",
                chunk.id, sample.type, sample.id)
      bg.dir <- dirname(bg.file)
      dir.create(bg.dir, showWarnings = FALSE, recursive = TRUE)
      cmd <-
        sprintf("./bigWigToBedGraph %s -chrom=%s -start=%d -end=%d %s",
                u, chunkChrom, chunkStart, chunkEnd, bg.file)
      system(cmd)
      bg.data <- fread(bg.file)
      setnames(bg.data, c("chrom", "chromStart", "chromEnd", "norm"))
      bg.data[, count := as.integer(norm/min(norm)) ]
      target.chromEnd <- with(bg.data, {
        unique(sort(c(chromStart[-1], chromEnd)))
      })
      target.chromStart <-
        c(bg.data$chromStart[1], target.chromEnd[-length(target.chromEnd)])
      bg.zero <- 
        data.frame(chromStart=target.chromStart,
                   chromEnd=target.chromEnd,
                   count=0L,
                   row.names=target.chromEnd)
      bg.zero[paste(bg.data$chromEnd), "count"] <- bg.data$count
      ggplot()+
        geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      ymin=0, ymax=norm),
                  color="grey",
                  data=data.frame(bg.data, what="norm"))+
        geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      ymin=0, ymax=count),
                  color="grey",
                  data=data.frame(bg.data, what="count"))+
        geom_step(aes(chromStart/1e3, count),
                  color="green",
                  data=data.frame(bg.data, what="count"))+
        geom_step(aes(chromStart/1e3, count),
                  data=data.frame(bg.zero, what="count"))+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(what ~ ., scales="free")
      chunk.sample.list[[paste(sample.type, sample.id)]] <-
        data.frame(sample.type, sample.id, bg.zero)
    }
  }
  region.list <- list()
  for(ann.i in 1:nrow(chunk.df)){
    chunk.row <- chunk.df[ann.i, ]
    type.vec <- type.list[[rownames(chunk.row)]]
    is.observed <- sample.types %in% type.vec
    observed <- sample.types[is.observed]
    not.observed <- sample.types[!is.observed]
    to.assign <- list()
    ann <- chunk.row$annotation
    to.assign[observed] <- ann
    to.assign[not.observed] <- "noPeaks"
    for(sample.type in names(to.assign)){
      sample.id <- names(sample.list[[sample.type]])
      annotation <- to.assign[[sample.type]]
      region.list[[paste(sample.type, ann.i)]] <- 
        data.frame(sample.type,
                   sample.id,
                   chromStart=chunk.row$chromStart,
                   chromEnd=chunk.row$chromEnd,
                   annotation)
    }
  }#ann.i
  chunk.data <- 
    list(counts=do.call(rbind, chunk.sample.list),
         regions=do.call(rbind, region.list))
  ggplot()+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.type + sample.id ~ ., scales="free")+
    scale_fill_manual(values=ann.colors)+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  data=chunk.data$regions,
                  alpha=0.5,
                  color="grey")+
    geom_step(aes(chromStart/1e3, count),
              color="grey50",
              data=chunk.data$counts)
  chunk.list[[chunk.name]] <- chunk.data
}

H3K27ac_TDH_control <- chunk.list

save(H3K27ac_TDH_control, file="H3K27ac_TDH_control.RData")
