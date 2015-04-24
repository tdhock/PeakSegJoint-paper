load("histone.sets.RData")
load("chunk.problems.RData")

is.histone <- names(chunk.problems) %in% names(histone.sets)
tf.set.names <- names(chunk.problems)[!is.histone]

train.sets <- histone.sets
set.seed(1)
for(set.name in tf.set.names){
  one.set <- chunk.problems[[set.name]]
  chunk.names <- names(one.set)
  for(split.i in 1:6){
    train.sets[[set.name]][[split.i]] <-
      sample(chunk.names, length(chunk.names)/2)
  }
}

save(train.sets, file="train.sets.RData")

