works_with_R("3.2.0",
             data.table="1.9.4",
             ggplot2="1.0")

## The script computes a list of PeakSegJoint problem regions for each
## data set and bases.per.problem size.

target.sizes <-
  c(H3K36me3=294912,
    H3K4me3=18432,
    nrsf=576,
    srf=576,
    max=576)

bases.per.problem.all <- 4.5 * 2^(5:20)

set.dirs <- Sys.glob("../chip-seq-paper/chunks/*_*_*") #include TF data!
for(set.dir.i in seq_along(set.dirs)){
  set.dir <- set.dirs[[set.dir.i]]
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  target.bases <- target.sizes[[experiment]]
  target.i <- which(bases.per.problem.all==target.bases)
  stopifnot(length(target.i) == 1)
  target.i.vec <- (target.i-2):(target.i+2)
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
      problemSeq <- seq(0, max.chromEnd, by=bases.per.problem)
      problemStart <-
        as.integer(sort(c(problemSeq,
                          problemSeq+bases.per.problem/2)))
      problemEnd <- problemStart+bases.per.problem
      is.overlap <- min.chromStart < problemEnd
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

      for(problem.i in 1:nrow(all.problems.dt)){
        problem <- all.problems.dt[problem.i, ]
        problem.name <- paste(problem$problem.name)
        problem.regions <- overlap.regions[problem.name]
        regions.by.sample <- split(problem.regions, problem.regions$sample.id)
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
        fit <- PeakSegJointHeuristic(profile.list)
        converted <- ConvertModelList(fit)
        peaks.by.peaks <- 
          c(list("0"=Peaks()),
            with(converted, split(peaks, peaks$peaks)))
        error.by.peaks <- list()
        for(peaks.str in names(peaks.by.peaks)){
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
        sapply(error.by.peaks, with, sum(weighted.error))
        sapply(error.by.peaks, with, sum(weight))
        stop("compute target interval.")
      }#problem.i
    }#bases.per.problem
  }#chunk.id
}#set.dir.i
