works_with_R("3.2.0",
             microbenchmark="1.3.0",
             Segmentor3IsBack="1.8",
             "tdhock/PeakSegDP@fe06a5b91d68c5d1ec471cb15c3ec3935dc2624d",
             "tdhock/PeakSegJoint@a1ea491f49e9bdb347f1caadebe7b750de807ac4")

simPeak <- function(n.data){
  stopifnot(is.numeric(n.data))
  n.data <- as.integer(n.data[1])
  c(rpois(n.data, 5),
    rpois(n.data, 10),
    rpois(n.data, 5))
}

simSamples <- function(n.samples, n.data){
  data.list <- list()
  for(sample.i in 1:n.samples){
    count <- simPeak(n.data)
    chromEnd <- seq_along(count) + 1000L
    data.list[[paste(sample.i)]] <- 
      data.frame(chromStart=chromEnd-1L,
                 chromEnd,
                 count)
  }
  data.list
}

seg.funs <-
  list(pDPA=function(data.list){
    result.list <- list()
    for(sample.i in seq_along(data.list)){
      one <- data.list[[sample.i]]
      fit <- Segmentor(one$count, Kmax=3, compress=FALSE)
      end.i <- fit@breaks[3, 1:2]
      result.list[[sample.i]] <- one$chromEnd[end.i]
    }
    as.integer(colMeans(do.call(rbind, result.list)))
  }, cDPA=function(data.list){
    result.list <- list()
    for(sample.i in seq_along(data.list)){
      one <- data.list[[sample.i]]
      fit <- cDPA(one$count, maxSegments=3)
      end.mat <- getPath(fit)
      end.i <- end.mat[3, 1:2]
      result.list[[sample.i]] <- one$chromEnd[end.i]
    }
    as.integer(colMeans(do.call(rbind, result.list)))
  }, PeakSegJoint=function(data.list){
    fit <- PeakSegJointHeuristic(data.list)
    models <- fit$models
    big <- models[[length(models)]]
    big$peak_start_end
  })

## First test: just one sample, varying data set size.
set.seed(1)
time.list <- list()
result.list <- list()
for(N in as.integer(10^seq(1, 6, by=0.25))){
  print(N)
  some <- simSamples(1, N)
  one <- some[[1]]
  target <- c(1000 + N, 1000 + N * 2)
  m.args <- list()
  algo.vec <- if(10^4 <= N)"PeakSegJoint" else names(seg.funs)
  N.result.list <- list()
  for(algorithm in algo.vec){
    m.args[[algorithm]] <- substitute({
      seg.fun <- seg.funs[[ALGO]]
      N.result.list[[ALGO]] <- seg.fun(some)
    }, list(ALGO=algorithm))
  }
  times <- microbenchmark(list=m.args, times=3)
  time.list[[paste(N)]] <- 
    data.frame(n.data=N, times)
  for(algorithm in names(N.result.list)){
    peak <- N.result.list[[algorithm]]
    seg.end.vec <- c(peak, max(one$chromEnd))
    seg.start.vec <- c(min(one$chromStart), peak)
    loss.list <- list()
    for(seg.i in seq_along(seg.end.vec)){
      seg.start <- seg.start.vec[[seg.i]]
      seg.end <- seg.end.vec[[seg.i]]
      i.end <- seg.end - 1000L
      i.start <- seg.start - 999L
      seg.count.vec <- one$count[i.start:i.end]
      seg.mean <- mean(seg.count.vec)
      loss.list[[seg.i]] <- PoissonLoss(seg.count.vec, seg.mean)
    }
    total.loss <- sum(unlist(loss.list))
    mean.loss <- total.loss/nrow(one)
    diff <- sum(abs(peak - target))
    result.list[[paste(algorithm, N)]] <-
      data.frame(algorithm=algorithm, n.data=N, diff, total.loss, mean.loss)
  }
}
one.sample.times <- do.call(rbind, time.list)
results <- do.call(rbind, result.list)

timings <-
  list(seconds=one.sample.times,
       results=results)

save(timings, file="timings.RData")

## Second test: several samples, one data set size.
## set.seed(1)
## time.list <- list()
## result.list <- list()
## N <- 100
## for(S in seq(1, 20, by=1)){
##   print(S)
##   some <- simSamples(S, N)
##   target <- c(1000 + N, 1000 + N * 2)
##   m.args <- list()
##   for(algorithm in names(seg.funs)){
##     m.args[[algorithm]] <- substitute({
##       seg.fun <- seg.funs[[ALGO]]
##       peak <- seg.fun(some)
##       diff <- sum(abs(peak - target))
##       result.list[[paste(ALGO, S)]] <-
##         data.frame(algorithm=ALGO, n.samples=S, diff)
##     }, list(ALGO=algorithm))
##   }
##   times <- microbenchmark(list=m.args, times=3)
##   time.list[[paste(S)]] <- 
##     data.frame(n.samples=S, times)
## }
## several.samples.times <- do.call(rbind, time.list)
## several.results <- do.call(rbind, result.list)

## ggplot()+
##   scale_x_log10()+
##   scale_y_log10("")+
##   geom_point(aes(n.samples, time/1e9, color=expr),
##              pch=1,
##              data=data.frame(several.samples.times, data="seconds"))+
##   scale_size_manual(values=c(PeakSegJoint=2, cDPA=2, pDPA=1))+
##   geom_line(aes(n.samples, diff, color=algorithm, size=algorithm),
##             data=data.frame(several.results, data="distance to true peak"))+
##   theme_bw()+
##   guides(size="none")+
##   theme(panel.margin=grid::unit(0, "cm"))+
##   facet_grid(data ~ ., scales="free")

