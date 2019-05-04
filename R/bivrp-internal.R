add_dplots_plot <-
  function(range.x, range.y, xlab., ylab.,
           d1, d2, cc1, c2, startd1, startd2,
           res1, res2, res.original1, res.original2, d1.or, d2.or,
           transparent.colors) {
    #axis(1, round(seq(range.x[1], range.x[2]/1.7, length=5), 1))
    #axis(2, round(seq(range.y[1], range.y[2]/1.7, length=5), 1))
    
    mtext(side=1, text=xlab., line=3, adj=.375)
    mtext(side=2, text=ylab., line=3, adj=.375)
    
    if(transparent.colors) {
      fill_col <- c("lightgray", "#00000050")
      border_col <- c("lightgray", "#00000050")
      lty_sim <- 1
    } else {
      fill_col <- c(NULL,NULL)
      border_col <- c(1,1)
      lty_sim <- 2
    }
    
    polygon(d1$x, (d1$y * c2) + startd1, col=fill_col[1], border=border_col[1], lty = lty_sim)
    polygon(d1.or$x, (d1.or$y * c2) + startd1, col=fill_col[2], border=border_col[2])
    points(res1, rep(startd1, length(res1)), pch=21, col="white", bg="lightgray", cex=.8)
    points(res.original1, rep(startd1, length(res.original1)), pch=16, col=1, cex=.7)
    
    polygon((d2$y * cc1) + startd2, d2$x, col=fill_col[1], border=border_col[1], lty = lty_sim)
    polygon((d2.or$y * cc1) + startd2, d2.or$x, col=fill_col[2], border=border_col[2])
    points(rep(startd2, length(res2)), res2, pch=21, col="white", bg="lightgray", cex=.8)
    points(rep(startd2, length(res.original2)), res.original2, pch=16, col=1, cex=.6) 
  }

add_dplots_prep <-
  function(res1, res2, res.original1, res.original2, density.bw) {
    d1 <- density(res1, bw = density.bw)
    d1.or <- density(res.original1, bw = density.bw)
    d1$y <- d1$y[d1$x < max(res1)*1.2]
    d1$y[length(d1$y)] <- 0
    d1$x <- d1$x[d1$x < max(res1)*1.2]
    d1.or$y <- d1.or$y[d1.or$x < max(res1, res.original1)*1.2]
    d1.or$y[length(d1.or$y)] <- 0
    d1.or$x <- d1.or$x[d1.or$x < max(res1, res.original1)*1.2]
    d2 <- density(res2, bw = density.bw)
    d2.or <- density(res.original2, bw = density.bw)
    d2$y <- d2$y[d2$x < max(res2)*1.2]
    d2$y[length(d2$y)] <- 0
    d2$x <- d2$x[d2$x < max(res2)*1.2]
    d2.or$y <- d2.or$y[d2.or$x < max(res2, res.original2)*1.2]
    d2.or$y[length(d2.or$y)] <- 0
    d2.or$x <- d2.or$x[d2.or$x < max(res2, res.original2)*1.2]
    startd1 <- max(res2, res.original2)*1.2
    startd2 <- max(res1, res.original1)*1.2
    cc1 <- 1
    ratio1 <- max((max(d2$y, d2.or$y) * cc1) + startd2)/(max(res1, res.original1)*1.7)
    if(ratio1 > 1) {
      cc1 <- .1
      ratio1 <- max((max(d2$y, d2.or$y) * cc1) + startd2)/(max(res1, res.original1)*1.7)
    }
    while(ratio1 < .9) {
      ratio1 <- max((max(d2$y, d2.or$y) * cc1) + startd2)/(max(res1, res.original1)*1.7)
      cc1 <- cc1 + .01
    }
    c2 <- 1
    ratio2 <- max((max(d1$y, d1.or$y) * c2) + startd1)/(max(res2, res.original2)*1.7)
    if(ratio2 > 1) {
      c2 <- .1
      ratio2 <- max((max(d1$y, d1.or$y) * c2) + startd1)/(max(res2, res.original2)*1.7)
    }
    while(ratio2 < .9) {
      ratio2 <- max((max(d1$y, d1.or$y) * c2) + startd1)/(max(res2, res.original2)*1.7)
      c2 <- c2 + .01
    }
    return(list("d1"=d1, "d2"=d2, "cc1"=cc1, "c2"=c2, "startd1"=startd1, "startd2"=startd2, "d1.or"=d1.or, "d2.or"=d2.or))
  }

chp_perpoint <-
  function(df.point, col.polygon, res.original,
           reduce.polygon, pch, cex, conf, ppoly, ...) {
    res1 <- df.point[,1]
    res2 <- df.point[,2]
    res.original1 <- res.original[1]
    res.original2 <- res.original[2]
    chp.xy <- chull(res1, res2)
    dxy <- data.frame(res1, res2)
    pol <- dxy[chp.xy,]
    pol.area <- polygon_area(pol)$area
    
    if(reduce.polygon == "peel") {
      dif <- 1
      while(dif > conf) {
        dxy <- dxy[-chp.xy,]
        chp.xy <- with(dxy, chull(res1, res2))
        pol <- dxy[chp.xy,]
        pol.area2 <- polygon_area(pol)$area
        dif <- pol.area2/pol.area
      }
    } else if(reduce.polygon == "bag") {
      pol <- get_reduced_bag(x = res1, y = res2, conf = conf)  
    } else {
      pol <- get_newpolygon(conf, pol, method = reduce.polygon)
    }
    
    pol.centre <- polygon_area(pol)$centre
    if(missing(col.polygon)) col.polygon <- "#C0C0C055"
    if(ppoly) polygon(pol, col=col.polygon, border="#00000055", lty=1)
    pol.inside <- is_point_inside(point=c(res.original1, res.original2), polyg=pol)
    ret <- list("pol.centre"=pol.centre, "pol.inside"=pol.inside)
    return(ret)
  }

sort_theta <-
  function(df.obj, reference, counter.clockwise) {
    n <- nrow(df.obj)
    x <- df.obj[,1]
    y <- df.obj[,2]
    theta <- atan2(y, x)
    if(counter.clockwise) {
      ord <- order(theta)
    } else {
      ord <- order(theta, decreasing = TRUE)
    }
    df <- data.frame(x, y, theta)
    df.ord <- df[ord,]
    if(counter.clockwise) {
      df.ord2 <- rbind(df.ord[df.ord$theta >= reference,],
                       df.ord[df.ord$theta < reference,])
    } else {
      df.ord2 <- rbind(df.ord[df.ord$theta < reference,],
                       df.ord[df.ord$theta >= reference,])
    }
    df.ret <- as.data.frame(df.ord2[,-3])
    #ft1 <- abs(df.ord2$theta[1] - reference)
    #ft2 <- abs(df.ord2$theta[n] - reference)
    #if(ft1 < ft2) df.ret <- df.ret[c(n,1:(n-1)),]
    return(df.ret)
  }
