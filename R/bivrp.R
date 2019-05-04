bivrp <-
function(obj, sim = 99, conf = .95, diagfun, simfun, fitfun, verb = FALSE,
         sort.res = TRUE, closest.angle = TRUE, angle.ref = - pi,
         counter.clockwise = TRUE, xlab, ylab, main,
         clear.device = FALSE, point.col, point.pch, ...) {
  
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
  
  nres <- length(res.original1)
  res1mat <- matrix(res1, byrow=F, nrow=nres)
  res2mat <- matrix(res2, byrow=F, nrow=nres)
  reslist <- list()
  for(i in 1:sim) {
    reslist[[i]] <- data.frame("res1"=res1mat[,i], "res2"=res2mat[,i]) 
  }
  
  if(sort.res) {
    if(closest.angle) {
      res.original.ord <- sort_theta(res.original.2, reference = angle.ref,
                                     counter.clockwise = counter.clockwise)
      ref <- atan2(res.original.ord[1,2], res.original.ord[1,1])
    } else {
      res.original.ord <- sort_theta(res.original.2, reference = angle.ref,
                                     counter.clockwise = counter.clockwise)
      ref <- angle.ref
    }
    reslist.ord <- lapply(reslist, sort_theta, reference = ref,
                          counter.clockwise = counter.clockwise)
  } else {
    reslist.ord <- lapply(reslist, function(obj) data.frame(x = obj[,1], y = obj[,2]))
    res.original.ord <- res.original.2
  }
  
  ret <- list("reslist.ord"=reslist.ord, "res.original.ord"=res.original.ord,
              "res1"=res1, "res2"=res2,
              "res.original1"=res.original1, "res.original2"=res.original2,
              "conf"=conf)
  class(ret) <- "bivrp"
  plot.bivrp(ret, xlab = xlab, ylab = ylab, main = main,
             point.col = point.col, point.pch = point.pch, ...)
  if(clear.device) dev.off()
  return(invisible(ret))
}

print.bivrp <- function(x, ...) {
  do.call(data.frame, x$reslist.ord)
}