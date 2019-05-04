\name{polygon-operations}
\alias{polygon_area}
\alias{get_k}
\alias{get_newpolygon}
\alias{get_reduced_bag}
\alias{compute_bagplot}
\encoding{UTF-8}
\title{Polygon operations}
\description{Convex polygon operations - determination of area, centre of mass, and area reduction}
\usage{
polygon_area(P)
get_k(P, conf)
get_newpolygon(conf, P, method)
get_reduced_bag(x, y, conf)
compute_bagplot(x, y, conf)
}
\arguments{
  \item{P}{2-column \code{matrix} or \code{data.frame} with the coordinates of the vertices of the convex polygon}
  \item{conf}{proportion of the area of polygon \code{P}}
  \item{method}{method used to reduce the area of the polygon. Use \code{method = "proportional"} to scale the distances between the centroid and the vertices by \code{sqrt(conf)}; use \code{method = "get_k"} to subtract the same distance \code{k} from the centroid to each vertex.}
  \item{x}{x coordinates (of raw data) used to obtain the reduced bag}
  \item{y}{y coordinates (of raw data) used to obtain the reduced bag}
}

\details{
The function \code{compute_bagplot} uses an adapted version of the code written by P. Segaert to obtain the bagplot, that uses the Fortran subroutine written by P.J. Rousseeuw, I. Ruts and A. Struyf.
}

\author{
Rafael A. Moral <rafael.deandrademoral@mu.ie> and John Hinde
}

\references{
Rousseeuw P.J., Ruts I., Tukey J.W. (1999). The bagplot: A bivariate boxplot. The American Statistician, 53, 382–387.
}

\seealso{
\code{\link{is_point_inside}}
\code{\link{polygon}}
}

\examples{
  oldPolygon <- data.frame(x=c(2,1,3,4.5,5), y=c(1,3,5,4.5,2))
  
  # area
  polygon_area(oldPolygon)$area
  # centre of mass
  polygon_area(oldPolygon)$centre
  
  # get a new polygon with 50% of the area of the old one
  newPolygon <- get_newpolygon(conf=.5, P=oldPolygon, method="proportional")
  polygon_area(newPolygon)$area/polygon_area(oldPolygon)$area
  
  # second method
  newPolygon2 <- get_newpolygon(conf=.5, P=oldPolygon, method="get.k")
  polygon_area(newPolygon2)$area/polygon_area(oldPolygon)$area
  
  # illustration
  plot(oldPolygon, xlim=c(0,6), ylim=c(0,6), main="(a)", pch=16)
  polygon(oldPolygon, lwd=2, col="#00000033")
  text(oldPolygon, c(expression(P[1]), expression(P[2]),
                     expression(P[3]), expression(P[4]),
                     expression(P[5])), pos=c(1,2,3,4,4), cex=2)
  polygon(newPolygon, border=4, lwd=2, col="#52A3E199")
  points(newPolygon, pch=16, col=4)
  text(newPolygon, c(expression(paste(P[1],minute)), expression(paste(P[2],minute)),
             expression(paste(P[3],minute)), expression(paste(P[4],minute)),
             expression(paste(P[5],minute))), pos=c(1,3,2,4,4), col=4, cex=2)

  C <- polygon_area(oldPolygon)$centre
  text(C[1], C[2], "C", pos=4, cex=2)
  for(i in 1:5) lines(c(C[1], oldPolygon[i,1]), 
                      c(C[2], oldPolygon[i,2]), lty=2, lwd=2, type="b")
}

\keyword{ polygon }