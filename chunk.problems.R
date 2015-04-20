works_with_R("3.2.0",
             ggplot2="1.0")

## The script computes a list of PeakSegJoint problem regions for each
## data set and bases.per.problem size.

chunk.dirs <- Sys.glob("../chip-seq-paper/chunks/*/") #include TF data!

positive.region.bases <- stop("compute median of positive regions")

bases.per.problem.all <- 4.5 * 2^(1:20)
bases.per.problem.hi <-
  bases.per.problem.all[bases.per.problem.all > positive.region.bases]
bases.per.problem.vec <- bases.per.problem.hi[1:6]

bases.per.problem <- bases.per.problem.vec[[bases.per.problem.i]]
cat(sprintf("%4d / %4d %s\n", bases.per.problem.i,
            length(bases.per.problem.vec), bases.per.problem))
problemSeq <- seq(0, chrom.range[2], by=bases.per.problem)
problemStart <-
  as.integer(c(problemSeq,
               problemSeq+bases.per.problem/3,
               problemSeq+2*bases.per.problem/3))
problemStart <-
  as.integer(c(problemSeq,
               problemSeq+bases.per.problem/4,
               problemSeq+2*bases.per.problem/4,
               problemSeq+3*bases.per.problem/4))
problemStart <-
  as.integer(c(problemSeq,
               problemSeq+bases.per.problem/2))
problemEnd <- problemStart+bases.per.problem
peakStart <- as.integer(problemStart+bases.per.problem/4)
peakEnd <- as.integer(problemEnd-bases.per.problem/4)
is.overlap <- chrom.range[1] < problemEnd
problem.name <- sprintf("chr21:%d-%d", problemStart, problemEnd)
problems <- 
  data.frame(problem.name, peakStart, peakEnd,
             bases.per.problem, problemStart, problemEnd)[is.overlap,]
all.problems.dt <- data.table(problems)
setkey(all.problems.dt, problemStart, problemEnd)

ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                alpha=0.5,
                color="grey",
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
setkey(all.problems.dt, problem.name)
