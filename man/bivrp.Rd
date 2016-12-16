\name{bivrp}
\alias{bivrp}
\encoding{UTF-8}

\title{Bivariate Residual Plots with Simulation Polygons}
\description{
description
}
\usage{
bivrp(obj, sim = 99, conf = 0.95, diagfun, simfun, fitfun, 
      verb = F, add.dplots = T, theta.sort = T, add.polygon = F, 
      reduce.polygon = T, kernel = F, chp = F, superpose.points = F, 
      one.dim = F, xlab, ylab, main, clear.device = F, ...)
}

\arguments{
  \item{obj}{fitted model object}
  \item{sim}{number of simulations used to compute envelope. Default is 99}
  \item{conf}{confidence level of the simulated polygons. Default is 0.95}
  \item{diagfun}{user-defined function used to obtain the diagnostic measures from the fitted model object}
  \item{simfun}{user-defined function used to simulate a random sample from the model estimated parameters}
  \item{fitfun}{user-defined function used to re-fit the model to simulated data}
  \item{verb}{logical. If \code{TRUE}, prints each step of the simulation procedure}
  \item{add.dplots}{logical. If \code{TRUE}, adds the marginal density plots}
  \item{theta.sort}{logical. If \code{TRUE}, produces a simulated polygon for each point}
  \item{add.polygon}{logical. If \code{TRUE}, plots the simulated polygons as well}
  \item{reduce.polygon}{logical. If \code{TRUE}, reduces the polygon area using the algorithm in \code{get.k} and \code{get.newpolygon}. If \code{FALSE}, performs convex hull peeling to reduce the area}
  \item{kernel}{logical. If TRUE, instead of using polygons for each point, computes 2d kernels and plots the contours}
  \item{chp}{logical. If \code{TRUE}, instead of using polygons for each point, performs convex hull peeling over all simulated points}
  \item{superpose.points}{only used if \code{kernel} or \code{chp} is TRUE. Logical argument, if \code{TRUE}, plots all simulated bivariate diagnostics}
  \item{one.dim}{logical. If \code{TRUE}, plots only the marginal density plots}
  \item{xlab}{argument passed to \code{par}}
  \item{ylab}{argument passed to \code{par}}
  \item{main}{argument passed to \code{par}}
  \item{clear.device}{logical. If \code{TRUE}, clears the plotting device after producing the bivariate residual plot with simulation polygons}
  \item{\dots}{further arguments passed to \code{par}}
}
\details{
This approach relies on the same strategy used for producing half-normal plots with simulation envelopes. Given a vector of bivariate model diagnostics, the angle each point makes with the origin is calculated to order them. Then, by default 99 bivariate response variables are simulated from the fitted model, using the same model matrices, error distribution and fitted parameters. The model is refitted to each simulated sample, obtaining the same type of model diagnostics, again ordered by the angle they form with the origin. We have, for each bivariate diagnostic, 99 simulated bivariate diagnostics forming the whole cloud of simulated diagnostics.  By default, we then obtain the convex hulls of each set of the $s$ sets of points and obtain a reduced polygon whose area is 95\% of the original convex hull's area, forming the simulated polygon. The points are then connected to the centroids of their respective simulated polygons and, if they lie outside the polygons, they are drawn in red. For the final display, the polygons are erased so as to ease visualization.

There is no automatic implementation of a bivariate model in this function, and hence users must provide three functions for \code{bivrp}. The first function, \code{diagfun}, must extract the desired model diagnostics from a model fit object. The second function, \code{simfun}, must return the response variable, simulated using the same error distributions and estimated parameters from the fitted model. The third and final function, \code{fitfun}, must return a fitted model object. See the Examples section.
}
\value{
The function returns an object of class "bivrp", which is a list containing the following components:

\item{reslist.ord }{list of ordered diagnostics from model refitting to each simulated dataset}
\item{res.original.ord }{original model diagnostics}
\item{res1 }{diagnostics from variable 1}
\item{res2 }{diagnostics from variable 2}
\item{add.polygon }{logical. Equals TRUE if \code{add.polygon=TRUE} in the bivrp call}
\item{res.original1 }{original model diagnostics for variable 1}
\item{res.original2 }{original model diagnostics for variable 2}
\item{theta.sort }{logical. Equals TRUE if \code{theta.sort=TRUE} in the bivrp call}
\item{conf }{confidence level of the simulated polygons}
\item{superpose.points }{logical. Equals TRUE if \code{superpose.points=TRUE} in the bivrp call}
\item{kernel }{logical. Equals TRUE if \code{kernel=TRUE} in the bivrp call}
\item{one.dim }{logical. Equals TRUE if \code{one.dim=TRUE} in the bivrp call}
\item{chp }{logical. Equals TRUE if \code{chp=TRUE} in the bivrp call}
\item{add.dplots }{logical. Equals TRUE if \code{add.dplots=TRUE} in the bivrp call}
\item{reduce.polygon }{logical. Equals TRUE if \code{reduce.polygon=TRUE} in the bivrp call}
}

\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, John Hinde and Clarice G. B. \enc{Demétrio}{Demetrio}
}

\examples{
## simulating a bivariate normal response variable

require(mvtnorm)
n <- 80
beta1 <- c(2, .4)
beta2 <- c(.2, .2)
x <- seq(1, 10, length=n)
X <- model.matrix(~ x)
mu1 <- X\%*\%beta1
mu2 <- X\%*\%beta2
sig1 <- 2
sig2 <- 3
sig12 <- -1.7
Sig1 <- diag(rep(sig1), n)
Sig2 <- diag(rep(sig2), n)
Sig12 <- diag(rep(sig12), n)
V <- rbind(cbind(Sig1, Sig12),
           cbind(Sig12, Sig2))

set.seed(2016)
Y <- as.numeric(rmvnorm(1, c(mu1, mu2), V))

## code for fitting the model estimating covariance or not
loglik.bn <- function(theta, Y, X, covariance=TRUE) {
  n <- length(Y)/2
  beta1 <- theta[1:2]
  beta2 <- theta[3:4]
  sig1 <- exp(theta[5])
  sig2 <- exp(theta[6])
  if(covariance) sig12 <- theta[7] else sig12 <- 0
  mu1 <- X\%*\%beta1
  mu2 <- X\%*\%beta2
  Sig1 <- diag(rep(sig1), n)
  Sig2 <- diag(rep(sig2), n)
  Sig12 <- diag(rep(sig12), n)
  V <- rbind(cbind(Sig1, Sig12),
             cbind(Sig12, Sig2))
  llik <- dmvnorm(Y, c(mu1, mu2), V, log=TRUE)
  return(-llik)
}

bivnormfit <- function(Y, X, covariance, init) {
  if(covariance) ni <- 7 else ni <- 6
  if(missing(init)) init <- rep(0, ni)
  fit <- optim(init, loglik.bn, Y=Y, X=X, covariance=covariance, method="BFGS")   
  coefs <- fit$par
  coefs[5:6] <- exp(coefs[5:6])
  fitted <- c(X\%*\%coefs[1:2], X\%*\%coefs[3:4])
  resid <- Y - fitted
  ret <- list("coefs"=coefs, "covariance"=covariance, "n"=length(Y)/2, 
              "X"=X, "fitted"=fitted, "resid"=resid, "loglik"=fit$value, "Y"=Y)
  class(ret) <- "bivnormfit"
  return(ret)
}

## fitting bivariate models with and without estimating covariance
fit0 <- bivnormfit(Y, X, covariance=FALSE)
fit1 <- bivnormfit(Y, X, covariance=TRUE)
## likelihood-ratio test
2*(fit0$loglik - fit1$loglik)
pchisq(43.33, 1, lower=FALSE)

## function for extracting diagnostics (raw residuals)
dfun <- function(obj) {
  r <- obj$resid
  n <- obj$n
  return(list(r[1:n], r[(n+1):(2*n)]))
}

## function for simulating new response variables
sfun <- function(obj) {
  n <- obj$n
  fitted <- obj$fitted
  sig1 <- obj$coefs[5]
  sig2 <- obj$coefs[6]
  if(obj$covariance) sig12 <- obj$coefs[7] else sig12 <- 0
  Sig1 <- diag(rep(sig1), n)
  Sig2 <- diag(rep(sig2), n)
  Sig12 <- diag(rep(sig12), n)
  V <- rbind(cbind(Sig1, Sig12),
             cbind(Sig12, Sig2))
  beta1 <- obj$coefs[1:2]
  beta2 <- obj$coefs[3:4]
  mu1 <- obj$X\%*\%beta1
  mu2 <- obj$X\%*\%beta2
  Y <- as.numeric(rmvnorm(1, c(mu1, mu2), V))
  return(list(Y[1:n], Y[(n+1):(2*n)], "X"=obj$X, "covariance"=obj$covariance))
}

## function for refitting the model to simulated data
ffun <- function(new.obj) {
  Ynew <- c(new.obj[[1]], new.obj[[2]])
  bivnormfit(Ynew, new.obj$X, new.obj$covariance)
}

\dontrun{
## Bivariate residual plot for model 1 (without estimating covariance)
plot1 <- bivrp(fit0, diagfun=dfun, simfun=sfun, fitfun=ffun, verb=T)
## without polygon area reduction
plot(plot1, conf=1)
## drawing polygons
plot(plot1, add.polygon=T)
## without ordering
plot(plot1, theta.sort=F, kernel=T, add.dplots=T, superpose=T)

## Bivariate residual plot for model 2 (estimating covariance)
plot2 <- bivrp(fit1, diagfun=dfun, simfun=sfun, fitfun=ffun, verb=T)
## without polygon area reduction
plot(plot2, conf=1)
## drawing polygons
plot(plot2, add.polygon=T, conf=1)
## without ordering
plot(plot2, theta.sort=F, kernel=T, add.dplots=T, superpose=T)
}
}