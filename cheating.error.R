works_with_R("3.2.0",
             dplyr="0.4.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegJoint@f1682b7e9a3d11c410c208325610ea1ede13edfa")

load("train.sets.RData")
load("selected.by.set.RData")

regions.file.list <- list()
regions.file.vec <- Sys.glob("../chip-seq-paper/chunks/*/*/regions.RData")
for(regions.file in regions.file.vec){
  load(regions.file)
  regions.file.list[[regions.file]] <- regions
}

RData.vec <- Sys.glob("chunk.problems/*/*/*.RData")
chunk.list <- list()
for(RData in RData.vec){
  objs <- load(RData)
  chunk.list[[RData]] <- best.step2.error
}

cheating.error.list <- list()
train.res.list <- list()
for(set.name in names(train.sets)){
  splits <- train.sets[[set.name]]
  selected.df <- selected.by.set[[set.name]]
  chunk.ids <- dir(file.path("chunk.problems", set.name))
  set.chunks <- paste0(set.name, "/", chunk.ids)
  stats.by.split <- list()
  for(split.i in 1:6){
    train.chunks <- splits[[split.i]]
    bases.per.problem <- selected.df[split.i, "bases.per.problem"]
    selected.res <- paste(bases.per.problem)
    is.train <- set.chunks %in% train.chunks
    test.chunks <- set.chunks[!is.train]
    split.name <- paste(set.name, "split", split.i)
    print(split.name)
    error.list <- list()
    for(chunk.name in test.chunks){
      RData.vec <- Sys.glob(sprintf("chunk.problems/%s/*.RData", chunk.name))
      for(RData in RData.vec){
        error.list[[RData]] <- chunk.list[[RData]]
      }
    }
    if(length(error.list) > 0){
      error.chunks <- do.call(rbind, error.list)
      error.split <- data.frame(set.name, split.i, error.chunks) %>%
        group_by(set.name, split.i, bases.per.problem) %>%
        summarise(errors=sum(errors),
                  regions=sum(regions)) %>%
        filter(regions == max(regions))
      train.res.list[[split.name]] <-
        subset(error.split, bases.per.problem==selected.res)
      cheating.error.list[[split.name]] <-
        subset(error.split, seq_along(errors) == which.min(errors))
    }
  }#split.i
}#set.name

cheating.error <- do.call(rbind, cheating.error.list)
best.for.train.res <- do.call(rbind, train.res.list)

save(cheating.error, best.for.train.res,
     file="cheating.error.RData")
