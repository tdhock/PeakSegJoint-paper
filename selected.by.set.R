works_with_R("3.2.0",
             data.table="1.9.4",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@9d1a4fc6751334513f2abb07689a187e04e4db76")

load("train.sets.RData")

RData.vec <- Sys.glob("PeakSegJoint-chunks/*/*/problems.RData")
chunk.list <- list()
for(RData in RData.vec){
  objs <- load(RData)
  chunk.list[[RData]] <- step2.error
}

pick.best.index <- structure(function
### Minimizer for local models, described in article section 2.3
### "Picking the optimal model"
(err
### Vector of errors to minimize.
 ){
  nparam <- length(err)
  candidates <- which(err==min(err))
  if(length(err)==1)return(candidates)
  st <- abs(median(candidates)-candidates)
  middle <- candidates[which.min(st)]
  if(all(diff(err)==0))return(middle)
  if(nparam %in% candidates && 1 %in% candidates){
    cat("Warning: strange error profile, picking something near the center\n")
    print(as.numeric(err))
    d <- diff(candidates)>1
    if(any(d)){
      which(d)[1]
    }else{
      middle
    }
  }else if(1 %in% candidates){
    max(candidates)
  }else if(nparam %in% candidates){
    min(candidates)
  }else {
    middle
  }
### Integer index of the minimal error.
},ex=function(){
  stopifnot(pick.best.index(rep(0,100))==50)

  err <- rep(1,100)
  err[5] <- 0
  stopifnot(pick.best.index(err)==5)

  ## should pick the middle
  err <- rep(1,100)
  err[40:60] <- 0
  stopifnot(pick.best.index(err)==50)

  ## should pick the biggest
  err <- rep(1,100)
  err[1:60] <- 0
  stopifnot(pick.best.index(err)==60)

  ## should pick the smallest
  err <- rep(1,100)
  err[50:100] <- 0
  stopifnot(pick.best.index(err)==50)
})

selected.by.set <- list()
error.by.set <- list()
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  best.list <- list()
  all.list <- list()
  for(split.i in seq_along(splits)){
    train.validation <- splits[[split.i]]
    error.list <- list()
    for(chunk.name in train.validation){
      RData.glob <-
        sprintf("PeakSegJoint-chunks/%s/problems.RData", chunk.name)
      RData.vec <- Sys.glob(RData.glob)
      for(RData in RData.vec){
        error.list[[RData]] <- chunk.list[[RData]]
      }
    }
    if(length(error.list) > 0){
      error.chunks <- do.call(rbind, error.list)
      error.split <- data.frame(split.i, error.chunks) %>%
        group_by(split.i, bases.per.problem) %>%
        summarise(fp=sum(fp),
                  fn=sum(fn),
                  errors=sum(errors),
                  regions=sum(regions))
      all.list[[paste(split.i)]] <-
        data.frame(set.name, error.split)
      regions.range <- error.split %>%
        summarise(min=min(regions),
                  max=max(regions))
      if(with(regions.range, min == max)){
        min.errors <- error.split %>%
          filter(errors == min(errors))
        error.split$status <- ""
        picked.i <- pick.best.index(error.split$errors)
        error.split$status[picked.i] <- "picked"
        if(nrow(min.errors) > 1){
          print(error.split)
        }
        best.list[[paste(split.i)]] <-
          error.split[picked.i, ]
      }
    }
  }
  if(length(best.list) > 0){
    set.best <- do.call(rbind, best.list)
    selected.by.set[[set.name]] <- data.frame(set.name, set.best)
    error.by.set[[set.name]] <- do.call(rbind, all.list)
  }
}#set.name
selected.error <- do.call(rbind, error.by.set)

save(selected.by.set,
     selected.error,
     file="selected.by.set.RData")
