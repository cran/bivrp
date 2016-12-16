bivrp <-
function(obj, sim=99, conf=.95, diagfun, simfun, fitfun, verb=F,
                  add.dplots=T, theta.sort=T, add.polygon=F, reduce.polygon=T,
                  kernel=F, chp=F, superpose.points=F, one.dim=F,
                  xlab, ylab, main, clear.device=F, ...) {

  res.original <- diagfun(obj)
  res.original1 <- res.original[[1]]
  res.original2 <- res.original[[2]]
  res1 <- res2 <- NULL
  
  for(i in 1:sim) {
    new.obj <- simfun(obj)
    new.fit <- fitfun(new.obj)
    res1 <- c(res1, diagfun(new.fit)[[1]])
    res2 <- c(res2, diagfun(new.fit)[[2]])
    if(verb) cat("Simulation", i, "out of", sim, "\n")
  }
  
  res.original.2 <- data.frame(res.original1, res.original2)
  res.original.ord <- sorttheta(res.original.2)
  nres <- length(res.original1)
  res1mat <- matrix(res1, byrow=F, nrow=nres)
  res2mat <- matrix(res2, byrow=F, nrow=nres)
  reslist <- list()
  for(i in 1:sim) {
    reslist[[i]] <- data.frame("res1"=res1mat[,i], "res2"=res2mat[,i]) 
  }
  reslist.ord <- lapply(reslist, sorttheta)

  ret <- list("reslist.ord"=reslist.ord, "res.original.ord"=res.original.ord,
              "res1"=res1, "res2"=res2, "add.polygon"=add.polygon,
              "res.original1"=res.original1, "res.original2"=res.original2, "theta.sort"=theta.sort,
              "conf"=conf, "superpose.points"=superpose.points, "kernel"=kernel, "one.dim"=one.dim, 
              "chp"=chp, "add.dplots"=add.dplots, "reduce.polygon"=reduce.polygon)
  class(ret) <- "bivrp"
  plot.bivrp(ret, kernel=kernel, one.dim=one.dim, chp=chp, add.dplots=add.dplots, theta.sort=theta.sort, reduce.polygon=reduce.polygon,
             superpose.points=superpose.points, xlab=xlab, ylab=ylab, main=main, ...)
  if(clear.device) dev.off()
  return(invisible(ret))
}
