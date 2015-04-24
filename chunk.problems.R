works_with_R("3.2.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@58a0014e9446996fb967c21826b196d5e9a3935e",
             data.table="1.9.4",
             ggplot2="1.0")

## The script computes a list of PeakSegJoint problem regions for each
## data set and bases.per.problem size.

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

tf.size <- 73728

target.sizes <-
  c(H3K36me3=294912,
    H3K4me3=18432,
    nrsf=tf.size,
    srf=tf.size,
    max=tf.size)

bases.per.problem.all <- 4.5 * 2^(5:20)

set.dir.i <- 10
chunk.id <- "9"
bases.per.problem <- 144
bases.per.problem <- 36864

weighted.error.list <- list()
chunk.problems <- list()
set.dirs <- Sys.glob("../chip-seq-paper/chunks/*_*_*") #include TF data!
for(set.dir.i in seq_along(set.dirs)){
  set.dir <- set.dirs[[set.dir.i]]
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  target.bases <- target.sizes[[experiment]]
  target.i <- which(bases.per.problem.all==target.bases)
  stopifnot(length(target.i) == 1)
  target.i.vec <- (target.i-3):(target.i+3)
  bases.per.problem.vec <- bases.per.problem.all[target.i.vec]
  chunk.ids <- dir(set.dir)
  for(chunk.id in chunk.ids){
    chunk.path <- file.path(set.dir, chunk.id)
    chunk.name <- paste0(set.name, "/", chunk.id)
    counts.RData <- file.path(chunk.path, "counts.RData")
    load(counts.RData)
    regions.RData <- file.path(chunk.path, "regions.RData")
    load(regions.RData)
    regions.dt <- data.table(regions)
    counts.dt <- data.table(counts)
    counts.dt[, count := coverage]
    setkey(regions.dt, chromStart, chromEnd)
    setkey(counts.dt, chromStart, chromEnd)
    regions.dt$region.i <- 1:nrow(regions.dt)
    chrom <- paste(regions$chrom[1])
    max.chromEnd <- max(counts$chromEnd)
    min.chromStart <- min(counts$chromStart)
    for(bases.per.problem in bases.per.problem.vec){
      res.str <- paste(bases.per.problem)
      RData.file <-
        sprintf("chunk.problems/%s/%s/%s.RData",
                set.name, chunk.id, bases.per.problem)
      if(file.exists(RData.file)){
        load(RData.file)
      }else{
        cat(sprintf("%s %s\n", chunk.name, bases.per.problem))
        problemSeq <- seq(0, max.chromEnd, by=bases.per.problem)
        problemStart <-
          as.integer(sort(c(problemSeq,
                            problemSeq+bases.per.problem/2)))
        problemEnd <- problemStart+bases.per.problem
        is.overlap <- min.chromStart < problemEnd &
          problemStart < max.chromEnd
        problem.name <- sprintf("%s:%d-%d", chrom, problemStart, problemEnd)
        problems <- 
          data.frame(problem.name, 
                     bases.per.problem, problemStart, problemEnd)[is.overlap,]
        all.problems.dt <- data.table(problems)
        setkey(all.problems.dt, problemStart, problemEnd)

        ggplot()+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3),
                        data=data.frame(chromStart=min.chromStart,
                          chromEnd=max.chromEnd),
                        fill="grey")+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                            fill=annotation),
                        data=regions)+
          scale_fill_manual(values=ann.colors)+
          geom_segment(aes(problemStart/1e3, problem.name,
                           xend=problemEnd/1e3, yend=problem.name),
                       data=problems)

        overlap.regions <- foverlaps(regions.dt, all.problems.dt, nomatch=0L)
        region.counts <- table(overlap.regions$region.i)
        region.weights <- 1/region.counts
        overlap.regions[, weight := region.weights[paste(region.i)]]
        table(paste(overlap.regions$problem.name))

        problems.with.regions <- paste(unique(overlap.regions$problem.name))
        setkey(overlap.regions, problem.name)

        problem.list <- list() #[[problem.name]]
        problem.stats.list <- list() #[[problem.name]]
        for(problem.i in 1:nrow(all.problems.dt)){
          problem <- all.problems.dt[problem.i, ]
          problem.name <- paste(problem$problem.name)
          problem.regions <- overlap.regions[problem.name]
          regions.by.sample <-
            split(problem.regions, problem.regions$sample.id, drop=TRUE)
          ## microbenchmark(foverlaps={
          ##   problem.counts <- foverlaps(counts.dt, problem, nomatch=0L)
          ##   print(nrow(problem.counts))
          ## },vector.ops={
          ##   problem.counts <-
          ##     counts.dt[! (chromEnd < problem$problemStart |
          ##                 problem$problemEnd < chromStart), ]
          ##   print(nrow(problem.counts))
          ## }, times=10)
          problem.counts <-
            counts.dt[! (chromEnd < problem$problemStart |
                           problem$problemEnd < chromStart), ]
          profile.list <- ProfileList(problem.counts)
          converted <- tryCatch({
            fit <- PeakSegJointHeuristic(profile.list)
            ConvertModelList(fit)
          }, error=function(e){
            list(loss=data.frame(peaks=0, loss=0))
          })
          peaks.by.peaks <- list("0"=Peaks())
          if(!is.null(converted$peaks)){
            peaks.list <- with(converted, split(peaks, peaks$peaks))
            peaks.by.peaks[names(peaks.list)] <- peaks.list
          }
          error.by.peaks <- list()
          if(length(regions.by.sample))for(peaks.str in names(peaks.by.peaks)){
            model.peaks <- peaks.by.peaks[[peaks.str]]
            peaks.by.sample <- split(model.peaks, model.peaks$sample.id)
            error.by.sample <- list()
            for(sample.id in names(regions.by.sample)){
              sample.regions <- regions.by.sample[[sample.id]]
              sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
                peaks.by.sample[[sample.id]]
              }else{
                Peaks()
              }
              sample.error <- PeakErrorChrom(sample.peaks, sample.regions)
              sample.error$weight <- sample.regions$weight
              sample.error$weighted.error <- with(sample.error, (fp+fn)*weight)
              error.by.sample[[sample.id]] <-
                data.frame(sample.id, peaks=as.numeric(peaks.str), sample.error)
            }#sample.id
            error.by.peaks[[peaks.str]] <- do.call(rbind, error.by.sample)
          }#peaks.str
          loss.df <- converted$loss
          if(length(regions.by.sample)){
            loss.df[names(error.by.peaks), "weighted.error"] <-
              sapply(error.by.peaks, with, sum(weighted.error))
            loss.df[names(error.by.peaks), "total.weight"] <-
              sapply(error.by.peaks, with, sum(weight))
          }else{
            loss.df$weighted.error <- 0
            loss.df$total.weight <- 0
          }
          loss.df$cummin <- cummin(loss.df$loss)
          some.loss <- subset(loss.df, loss==cummin & !is.na(total.weight))
          exact <-
            with(some.loss, exactModelSelection(loss, peaks, peaks))
          exact$weighted.error <- loss.df[paste(exact$peaks), "weighted.error"]
          indices <- with(exact, {
            largestContinuousMinimum(weighted.error,
                                     max.log.lambda-min.log.lambda)
          })
          total.weight <- some.loss$total.weight[1]
          stopifnot(total.weight == some.loss$total.weight)
          problem.stats.list[[problem.name]] <-
            data.frame(problem.name,
                       weighted.error=min(some.loss$weighted.error),
                       total.weight)
          problem.list[[problem.name]] <- 
            list(features=featureMatrix(profile.list),
                 modelSelection=exact,
                 peaks=peaks.by.peaks,
                 target=c(min.log.lambda=exact$min.log.lambda[indices$start],
                   max.log.lambda=exact$max.log.lambda[indices$end]))
        }#problem.i
        problem.stats <- do.call(rbind, problem.stats.list)
        total.weight <- sum(problem.stats$total.weight)
        stopifnot(all.equal(nrow(regions), total.weight))
        error.row <- 
          data.frame(set.name, chunk.name, bases.per.problem,
                     weighted.error=sum(problem.stats$weighted.error),
                     total.weight)
        RData.dir <- dirname(RData.file)
        dir.create(RData.dir, showWarnings=FALSE, recursive=TRUE)
        save(problem.list, error.row, file=RData.file)
      }
      
      chunk.problems[[set.name]][[chunk.name]][[res.str]] <- problem.list
      weighted.error.list[[set.name]][[chunk.name]][[res.str]] <- error.row
    }#bases.per.problem
  }#chunk.id
}#set.dir.i

save(chunk.problems, weighted.error.list, file="chunk.problems.RData")
