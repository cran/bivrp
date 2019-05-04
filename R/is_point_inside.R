is_point_inside <-
function(point, polyg) {
  p <- as.numeric(point)
  #library(mgcv)
  #return(in.out(as.matrix(polyg), p))
  is.vertex <- sum(apply(polyg, 1, function(x) all(x == p)))
  if(is.vertex==1) return(FALSE)
  sx1 <- sum(! p[1] <= polyg[,1])
  sx2 <- sum(! p[1] >= polyg[,1])
  if(sx1==0 | sx2==0) return(FALSE)
  sy1 <- sum(! p[2] <= polyg[,2])
  sy2 <- sum(! p[2] >= polyg[,2])
  if(sy1==0 | sy2==0) return(FALSE)
  px <- polyg[,1]
  py <- polyg[,2]
  px <- c(px, px[1])
  py <- c(py, py[1])
  segments <- (p[2] - py)/(p[1] - px)
  condition <- segments[-1] == segments[-length(segments)]
  is.segment <- sum(condition)
  if(is.segment>0) {
    c1 <- px[-1][condition]
    c2 <- px[-length(px)][condition]
    if(c1 < c2) {
      if(c1 < p[1] & p[1] < c2) return(FALSE)
    } else {
      if(c2 < p[1] & p[1] < c1) return(FALSE)
    }
  }
  xcross <- polyg[,1] + (px[-1] - px[-length(px)])*(p[2] - polyg[,2])/(py[-1] - py[-length(py)])  
  xcross <- c(xcross[length(xcross)], xcross[-length(xcross)])
  px2 <- polyg[,1]
  px2 <- c(px2[length(px2)], px2[-length(px2)])  
  px1 <- px[-length(px)]
  xcross2 <- xcross[p[1] <= xcross]
  px1 <- px1[p[1] <= xcross]
  px2 <- px2[p[1] <= xcross]
  n.intersections1 <- sum(px1 <= xcross2 & xcross2 <= px2)
  n.intersections2 <- sum(px2 <= xcross2 & xcross2 <= px1)
  n.int <- n.intersections1 + n.intersections2
  if(any(p[2] == polyg[,2][p[1] < polyg[,1]])) n.int <- n.int - 1
  if(n.int==1) return(TRUE) else return(FALSE)
}
