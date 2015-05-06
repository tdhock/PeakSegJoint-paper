works_with_R("3.2.0",
             ggplot2="1.0",
             directlabels="2014.6.13",
             "tdhock/PeakSegJoint@ff5a7c58e297b54b328047f4e02285f0cb5d2838")

data(H3K4me3.PGP.immune.4608)
chrom.vec <- sub(":.*", "", names(H3K4me3.PGP.immune.4608))
table(chrom.vec)
train.chroms <- c("chr1", "chr9")
sets <-
  list(train=chrom.vec %in% train.chroms,
       validation=! chrom.vec %in% train.chroms)
train.problems <- H3K4me3.PGP.immune.4608[sets$train]
fit <-
  IntervalRegressionProblems(train.problems,
                             initial.regularization=0.01,
                             factor.regularization=1.1,
                             threshold=1e-3)
set.error.list <- list()
for(set.name in names(sets)){
  in.set <- sets[[set.name]]
  problem.list <- H3K4me3.PGP.immune.4608[in.set]
  error.mat.list <- list()
  for(problem.name in names(problem.list)){
    problem <- problem.list[[problem.name]]
    pred.log.lambda <- fit$predict(problem$features)
    too.hi <- problem$target[2] < pred.log.lambda
    too.lo <- pred.log.lambda < problem$target[1]
    is.error <- too.hi | too.lo
    error.mat.list[[problem.name]] <- is.error
  }
  error.mat <- do.call(rbind, error.mat.list)
  percent.error <- colMeans(error.mat) * 100
  set.error.list[[set.name]] <-
    data.frame(set.name,
               regularization=fit$regularization,
               percent.error)
}
min.validation <- 
  subset(set.error.list$validation,
         percent.error==min(percent.error))
best.models <- fit$weight.mat[, rownames(min.validation)]
best.nonzero <- best.models[apply(best.models!=0, 1, any), ]
print(best.nonzero)

regularization.range <- with(min.validation, {
  data.frame(min=min(regularization), max=max(regularization))
})

weight.list <- list()
for(regularization.i in seq_along(fit$regularization.vec)){
  regularization <- fit$regularization.vec[[regularization.i]]
  weight <- fit$weight.mat[, regularization.i]
  weight.list[[regularization.i]] <-
    data.frame(regularization, weight, variable=names(weight))
}
weight.df <- do.call(rbind, weight.list)

set.error <- do.call(rbind, set.error.list)
error.lines <- data.frame(set.error, what="percent incorrect targets")
learned.weights <- data.frame(weight.df, what="learned weights")
nonzero.weights <- subset(learned.weights, weight != 0)
selected.vars <- subset(learned.weights, variable %in% rownames(best.nonzero))

pathPlot <- 
ggplot()+
  guides(color="none")+
  ggtitle("Non-zero weights in models with min validation error")+
  geom_line(aes(-log10(regularization), weight,
                group=variable, color=variable),
            data=learned.weights)+
  geom_point(aes(-log10(regularization), weight,
                color=variable),
             pch=1,
             data=nonzero.weights)+
  geom_dl(aes(-log10(regularization), weight,
              color=variable, label=variable),
          data=selected.vars, method="last.polygons")+
  ylab("")+
  scale_x_continuous("model complexity -log10(regularization)",
                     limits=c(min(-log10(learned.weights$regularization)),
                       2.5))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")+
  ## geom_point(aes(-log10(regularization), percent.error),
  ##            data=data.frame(min.validation, what="percent incorrect targets"))+
  geom_tallrect(aes(xmin=-log10(min), xmax=-log10(max)),
                alpha=0.2,
                data=regularization.range)+
  geom_dl(aes(-log10(regularization), percent.error, label=set.name),
          data=error.lines, method="last.qp")+
  geom_line(aes(-log10(regularization), percent.error,
                group=set.name, linetype=set.name),
            data=error.lines,
            show_guide=FALSE)

pdf("figure-lasso-path.pdf", height=5)
print(pathPlot)
dev.off()
