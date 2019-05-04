\name{plot.bivrp}
\alias{plot.bivrp}
\encoding{UTF-8}

\title{Plot Method for bivrp Objects}
\description{
Plots the bivariate residual plot with simulation polygons from a \code{bivrp} object
}
\usage{\method{plot}{bivrp}(x, kernel, superpose.points, chp, add.dplots, 
     theta.sort, add.polygon, reduce.polygon, one.dim, pch = 16, cex = 0.8, 
     conf, xlab, ylab, main, point.col, point.pch, transparent.colors,
     density.bw, ...)
}

\arguments{
  \item{x}{object of class \code{bivrp}}
  \item{kernel}{logical. If TRUE, instead of using polygons for each point, computes 2d kernels and plots the contours}
  \item{superpose.points}{only used if \code{kernel} or \code{chp} is TRUE. Logical argument, if \code{TRUE}, plots all simulated bivariate diagnostics}
  \item{chp}{logical. If \code{TRUE}, instead of using polygons for each point, performs convex hull peeling over all simulated points}
  \item{add.dplots}{logical. If \code{TRUE}, adds the marginal density plots}
  \item{theta.sort}{logical. If \code{TRUE}, produces a simulated polygon for each point}
  \item{add.polygon}{logical. If \code{TRUE}, plots the simulated polygons as well}
  \item{reduce.polygon}{method used to reduce the polygon area. Defaults to \code{proportional}, see \code{\link{get_newpolygon}} for details. If \code{reduce.polygon = "peel"}, performs convex hull peeling to reduce the area; if \code{reduce.polygon = "bag"}, computes a (\code{conf * 100})\% bagplot of the points}
  \item{one.dim}{logical. If \code{TRUE}, plots only the marginal density plots (only works with \code{theta.sort = FALSE})}
  \item{pch}{argument passed to \code{par}}
  \item{cex}{argument passed to \code{par}}
  \item{conf}{confidence level of the simulated polygons. Default is 0.95}
  \item{xlab}{argument passed to \code{par}}
  \item{ylab}{argument passed to \code{par}}
  \item{main}{argument passed to \code{par}}
  \item{point.col}{a vector of length 2 with the colors of the points that are inside and outside of the simulated polygons}
  \item{point.pch}{a vector of length 2 with the point characters of the points that are inside and outside of the simulated polygons}
  \item{transparent.colors}{logical. If \code{TRUE}, adds transparency to the marginal density plots; if \code{FALSE}, only the border lines are drawn}
  \item{density.bw}{the smoothing bandwidth to be used for the marginal densities. Defaults to \code{"SJ"} (see \code{\link{density}})}
  \item{\dots}{further arguments passed to \code{par}}
}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br> and John Hinde
}

\seealso{
\code{\link{bivrp}}
}