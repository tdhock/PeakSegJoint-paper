works_with_R("3.2.0",
             data.table="1.9.4",
             ggplot2="1.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@f1682b7e9a3d11c410c208325610ea1ede13edfa")

load("train.sets.RData")

RData.vec <- Sys.glob("chunk.problems/*/*/*.RData")
chunk.list <- list()
for(RData in RData.vec){
  objs <- load(RData)
  chunk.list[[RData]] <- best.step2.error
}

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
      RData.vec <- Sys.glob(sprintf("chunk.problems/%s/*.RData", chunk.name))
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
        if(nrow(min.errors) > 1){
          print(error.split)
        }
        best.list[[paste(split.i)]] <- min.errors[1,]
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
