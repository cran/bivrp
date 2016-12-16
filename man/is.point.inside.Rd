\name{is.point.inside}
\alias{is.point.inside}
\encoding{UTF-8}

\title{Determine if point is inside or outside a simple polygon area}
\description{Returns whether a point is inside or outside the convex polygon formed with the coordinates in a data frame or matrix}
\usage{
is.point.inside(point, polyg)
}

\arguments{
  \item{point}{vector of two values for a point in the Cartesian plane}
  \item{polyg}{data frame or matrix with the coordinates forming the convex polygon}
}
\details{
The algorithm used here draws a ray from the point and counts the number of intersections made with the polygon. If the number of intersections is only one, then this means the point is inside the convex polygon.
}
\value{
This function returns \code{TRUE}, if the point is inside and \code{FALSE}, otherwise.
}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, John Hinde and Clarice G. B. \enc{Demétrio}{Demetrio}
}

\examples{
my.polygon <- data.frame(c(1, 2, 3, 4, 3),
                         c(1, 0, .5, 3, 4))
points.to.test <- list(c(0, 0), c(2.5, 1), c(3.5, 4))

unlist(lapply(points.to.test, is.point.inside, my.polygon))

}
