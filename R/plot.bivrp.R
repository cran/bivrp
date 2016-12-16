plot.bivrp <-
function(x, kernel, superpose.points, chp, add.dplots, 
         theta.sort, add.polygon, reduce.polygon, one.dim, 
         pch=16, cex=.8, conf, xlab, ylab, main, ...) {
  
  obj <- x
  rm(x)
  if(missing(kernel)) kernel <- obj$kernel
  if(missing(chp)) chp <- obj$chp
  if(missing(superpose.points)) superpose.points <- obj$superpose.points
  if(missing(conf)) conf <- obj$conf
  if(missing(one.dim)) one.dim <- obj$one.dim
  if(missing(add.dplots)) add.dplots <- obj$add.dplots
  if(missing(theta.sort)) theta.sort <- obj$theta.sort
  if(missing(reduce.polygon)) reduce.polygon <- obj$reduce.polygon
  if(one.dim) {
    if(missing(xlab)) xlab <- "Residuals"
    if(missing(ylab)) ylab <- "Density"
    if(missing(main)) main <- ""
  } else {
    if(missing(xlab)) xlab <- "Residuals 1"
    if(missing(ylab)) ylab <- "Residuals 2"
    if(missing(main)) main <- ""
  }

  if(missing(add.polygon)) add.polygon <- obj$add.polygon
  reslist.ord <- obj$reslist.ord
  res.original.ord <- obj$res.original.ord
  res.simlist <- list()
  for(i in 1:nrow(res.original.ord)) {
    res.simlist[[i]] <- matrix(unlist(
      lapply(reslist.ord, function(x) x[i,])), 
      byrow=T, ncol=2)
  }

  rxs <- range(res1 <- unlist(lapply(reslist.ord, function(a) a$x)))
  rys <- range(res2 <- unlist(lapply(reslist.ord, function(a) a$y)))
  rxo <- range(res.original1 <- res.original.ord[,1])
  ryo <- range(res.original2 <- res.original.ord[,2])
  
  if(add.dplots) {
    preps <- add.dplots.prep(res1=res1, res2=res2, 
                             res.original1=res.original1, res.original2=res.original2)
    d1 <- preps$d1
    d2 <- preps$d2
    cc1 <- preps$cc1
    c2 <- preps$c2
    startd1 <- preps$startd1
    startd2 <- preps$startd2
    d1.or <- preps$d1.or
    d2.or <- preps$d2.or
    range.x <- c(min(rxs, rxo), max(rxs, rxo)*1.7)
    range.y <- c(min(rys, ryo), max(rys, ryo)*1.7)
    axes. <- F
    bty. <- "n"
    xlab. <- xlab
    ylab. <- ylab
    xlab <- ""
    ylab <- ""
  } else {
    range.x <- c(min(rxs, rxo), max(rxs, rxo))
    range.y <- c(min(rys, ryo), max(rys, ryo))
    axes. <- T
    bty. <- "o"
  }
  
  if(theta.sort) {
    plot(0, 0, xlim=range.x, ylim=range.y, type="n",
         xlab=xlab, ylab=ylab, main=main, axes=axes., bty=bty., ...)
    pol.centre <- matrix(0, ncol=2, nrow=nrow(res.original.ord))
    pol.inside <- NULL
    if(add.polygon) ppoly <- T else ppoly <- F
    for(i in 1:nrow(res.original.ord)) {
      chpoly <- chp.perpoint(df.point=res.simlist[[i]], 
                       res.original=res.original.ord[i,], reduce.polygon=reduce.polygon,
                       pch=pch, cex=cex, conf=conf, ppoly=ppoly, ...)
      pol.centre[i,] <- chpoly$pol.centre
      pol.inside[i] <- chpoly$pol.inside
    }
    point.colours <- ifelse(pol.inside, 1, 2)
    points(res.original.ord, pch=pch, cex=cex, col=point.colours, ...)
    points(pol.centre, pch=1, cex=cex, col=point.colours, ...)
    arrows(pol.centre[,1], pol.centre[,2], 
           as.numeric(res.original.ord[,1]), as.numeric(res.original.ord[,2]),
           code=3, length=0, lty=2, col=point.colours)
    tot.points <- nrow(res.original.ord)
    n.out <- sum(!pol.inside)
    cat(n.out, " out of ", tot.points, " points out of polygons (", round(n.out/tot.points*100, 2), "%).", "\n", sep="")
    
    if(add.dplots) {
      add.dplots.plot(range.x, range.y, xlab., ylab.,
                      d1, d2, cc1, c2, startd1, startd2,
                      res1, res2, res.original1, res.original2, d1.or, d2.or)
    }
    
  } else {
    
    res1 <- obj$res1
    res2 <- obj$res2
    res.original1 <- obj$res.original1
    res.original2 <- obj$res.original2
  
    if(one.dim) {
      par(mfrow=c(1,2))
      dres1 <- density(res1)
      dres1.or <- density(res.original1)
      plot(res1, rep(0, length(res1)), xlim=range(res1, res.original1), pch=16, col="lightgray", ylim=c(0, max(dres1$y, dres1.or$y)), type="n", xlab=xlab, ylab=ylab, main=main)
      polygon(dres1, col="lightgray", border="lightgray")
      polygon(dres1.or, col="#00000050", border="#00000050")
      points(res1, rep(0, length(res1)), pch=21, col="white", bg="lightgray")
      points(res.original1, rep(0, length(res.original1)), pch=16, col=1)

      dres2 <- density(res2)
      dres2.or <- density(res.original2)
      plot(res2, rep(0, length(res2)), xlim=range(res2, res.original2), pch=16, col="lightgray", ylim=c(0, max(dres2$y, dres2.or$y)), type="n", xlab=xlab, ylab=ylab, main=main)
      polygon(dres2, col="lightgray", border="lightgray")
      polygon(dres2.or, col="#00000050", border="#00000050")
      points(res2, rep(0, length(res2)), pch=21, col="white", bg="lightgray")
      points(res.original2, rep(0, length(res.original2)), pch=16, col=1)
      return(invisible()) 
    }

    
    if(chp) {
      chp.xy <- chull(res1, res2)
      dxy <- data.frame(res1, res2)
      pol <- dxy[chp.xy,]
      pol.area <- polygon.area(pol)$area
      if(reduce.polygon) {
        k <- get.k(pol, conf)
        pol <- get.newpolygon(k, pol)
      } else {
        dif <- 1
        while(dif > conf) {
          dxy <- dxy[-chp.xy,]
          chp.xy <- with(dxy, chull(res1, res2))
          pol <- dxy[chp.xy,]
          pol.area2 <- polygon.area(pol)$area
          dif <- pol.area2/pol.area
        }
      }
      plot(res1, res2, xlim=range.x, ylim=range.y, type="n",
           xlab=xlab, ylab=ylab, main=main, axes=axes., bty=bty., ...)
      col.polygon <- "lightgray"
      if(superpose.points) {
        points(res1, res2, col="lightgray", pch=16)
        col.polygon <- NA
      }
      polygon(pol, col=col.polygon, lty=2)
      points(res.original1, res.original2, pch=pch, cex=cex, ...)
      if(add.dplots) {
        add.dplots.plot(range.x, range.y, xlab., ylab.,
                        d1, d2, cc1, c2, startd1, startd2,
                        res1, res2, res.original1, res.original2, d1.or, d2.or)
      }
      return(invisible()) 
    }
    
    if(kernel) {
      #library(MASS)
      plot(res1, res2, xlim=range.x, ylim=range.y, type="n",
           xlab=xlab, ylab=ylab, main=main, axes=axes., bty=bty., ...)
      post1 <- kde2d(res1, res2, n=100, lims=c(range(res1)*1.2, range(res2)*1.2))
      dx <- diff(post1$x[1:2])
      dy <- diff(post1$y[1:2])
      sz <- sort(post1$z)
      c1 <- cumsum(sz)*dx*dy
      level <- approx(c1, sz, xout=1-conf)$y
      ckern <- contourLines(post1$x, post1$y, post1$z, levels=level)
      col.polygon <- "lightgray"
      if(superpose.points) {
        points(res1, res2, col="lightgray", pch=16)
        col.polygon <- NA
      }
      for(i in 1:length(ckern)) polygon(ckern[[i]]$x, ckern[[i]]$y, col=col.polygon, lty=2)
      points(res.original1, res.original2, pch=pch, cex=cex, ...)
      if(add.dplots) {
        add.dplots.plot(range.x, range.y, xlab., ylab.,
                        d1, d2, cc1, c2, startd1, startd2,
                        res1, res2, res.original1, res.original2, d1.or, d2.or)
      }
    }
  }
}
