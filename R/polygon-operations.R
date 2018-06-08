polygon.area <-
function(P) {
  if(dim(P)[2]!=2) stop("supply data.frame or matrix with two columns")
  x <- P[,1]
  x <- c(x, x[1])
  y <- P[,2]
  y <- c(y, y[1])
  xy <- x[-length(x)] * y[-1]
  yx <- y[-length(y)] * x[-1]
  x12 <- x[-length(x)] + x[-1]
  y12 <- y[-length(y)] + y[-1]
  area <- sum(xy - yx)/2
  centre.x <- 1/(6*area) * sum(x12 * (xy - yx))
  centre.y <- 1/(6*area) * sum(y12 * (xy - yx))
  return(list("area"=abs(area), "centre"=c(centre.x, centre.y)))
}

get.k <-
  function(P, conf) {
    AP <- polygon.area(P)$area
    C <- polygon.area(P)$centre
    di <- apply(P, 1, function(x) dist(rbind(x, C)))
    di <- c(di, di[1])
    xtil <- P[,1] - C[1]
    ytil <- P[,2] - C[2]
    xtil <- c(xtil, xtil[1])
    ytil <- c(ytil, ytil[1])
    xi <- P[,1]
    xi <- c(xi, xi[1])
    yi <- P[,2]
    yi <- c(yi, yi[1])
    
    a <- sum((xtil[-length(di)]*ytil[-1] - xtil[-1]*ytil[-length(di)])/(di[-length(di)]*di[-1]))
    b <- sum((di[-1]*xi[-1]*ytil[-length(di)] +
                di[-length(di)]*xtil[-1]*yi[-length(di)] -
                di[-1]*xtil[-length(di)]*yi[-1] -
                di[-length(di)]*xi[-length(di)]*ytil[-1])/(di[-length(di)]*di[-1]))
    c1 <- sum(xi[-length(di)]*yi[-1] - xi[-1]*yi[-length(di)]) + 2*conf*AP
    c2 <- sum(xi[-length(di)]*yi[-1] - xi[-1]*yi[-length(di)]) - 2*conf*AP
    
    delta1 <- b^2 - 4*a*c1
    delta2 <- b^2 - 4*a*c2
    
    options(warn=-1)
    k1 <- (-b+sqrt(delta1))/(2*a)
    k2 <- (-b-sqrt(delta1))/(2*a)
    k3 <- (-b+sqrt(delta2))/(2*a)
    k4 <- (-b-sqrt(delta2))/(2*a)
    options(warn=0)
    
    K <- c(k1, k2, k3, k4)
    return(min(K[!is.na(K)]))
  }

get.newpolygon <-
  function(conf, P, method = c("proportional", "get.k")) {
    C <- polygon.area(P)$centre
    di <- apply(P, 1, function(x) dist(rbind(x, C)))
    dil <- switch(method, "proportional" = di * sqrt(conf),
                          "get.k" = di - get.k(P, conf))
    yil <- dil/di*(P[,2] - C[2]) + C[2]
    xil <- dil/di*(P[,1] - C[1]) + C[1]
    newPoly <- data.frame(xil, yil)
    return(newPoly)
  }

#pol <- data.frame(x=c(2,1,3,4.5,5), y=c(-10,3,5,4.5,2))
#np1 <- get.newpolygon(conf = .7, P = pol, method = "get.k")
#np2 <- get.newpolygon(conf = .7, P = pol, method = "proportional")
#pC <- polygon.area(pol)$centre
#apply(pol, 1, function(x) arrows(x[1], x[2], pC[1], pC[2], length = 0))
#points(pC[1], pC[2], pch = 16)

#plot(pol, asp = T)
#polygon(pol)
#polygon(np1, lty = 2, col = "#00FF0055")
#polygon(np2, lty = 2, col = "#FF000055")
#polygon.area(np1)
#polygon.area(np2)
#points(polygon.area(pol)$centre[1], polygon.area(pol)$centre[2],
#       pch = 21, bg = 1, cex = 1.5)
#points(polygon.area(np1)$centre[1], polygon.area(np1)$centre[2],
#       pch = 21, bg = 3)
#points(polygon.area(np2)$centre[1], polygon.area(np2)$centre[2],
#       pch = 21, bg = 2)