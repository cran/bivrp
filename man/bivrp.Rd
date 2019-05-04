\name{bivrp}
\alias{bivrp}
\alias{print.bivrp}
\encoding{UTF-8}

\title{Bivariate Residual Plots with Simulation Polygons}
\description{
Produces a bivariate residual plot with simulation polygons to assess goodness-of-fit of bivariate statistical models, provided the user supplies three functions: one to obtain model diagnostics, one to simulate data from a fitted model object, and one to refit the model to simulated data.
}
\usage{
bivrp(obj, sim = 99, conf = .95, diagfun, simfun, fitfun, verb = FALSE,
      sort.res = TRUE, closest.angle = TRUE, angle.ref = - pi,
      counter.clockwise = TRUE, xlab, ylab, main,
      clear.device = FALSE, point.col, point.pch, ...)
      
\method{print}{bivrp}(x, ...)
}

\arguments{
  \item{obj}{fitted model object}
  \item{sim}{number of simulations used to compute envelope. Default is 99}
  \item{conf}{confidence level of the simulated polygons. Default is 0.95}
  \item{diagfun}{user-defined function used to obtain the diagnostic measures from the fitted model object}
  \item{simfun}{user-defined function used to simulate a random sample from the model estimated parameters}
  \item{fitfun}{user-defined function used to re-fit the model to simulated data}
  \item{verb}{logical. If \code{TRUE}, prints each step of the simulation procedure}
  \item{sort.res}{logical. If \code{TRUE}, points will be sorted using angles formed with the origin (type of ordering can be fine-tuned with arguments \code{closest.angle}, \code{angle.ref} and \code{counter.clockwise}).}
  \item{closest.angle}{logical. If \code{FALSE}, points will be sorted starting from the angle defined in \code{angle.ref}, if \code{TRUE}, points will be sorted starting from the closest angle to the observed bivariate sample ranked as first}
  \item{angle.ref}{the reference angle from which points will be sorted starting from the closest angle to the input (in radians). Defaults to \code{- pi}}
  \item{counter.clockwise}{logical. Should the points be ordered counter-clockwise or clockwise from the reference angle?}
  \item{xlab}{argument passed to \code{par}}
  \item{ylab}{argument passed to \code{par}}
  \item{main}{argument passed to \code{par}}
  \item{clear.device}{logical. If \code{TRUE}, clears the plotting device after producing the bivariate residual plot with simulation polygons}
  \item{point.col}{a vector of length 2 with the colors of the points that are inside and outside of the simulated polygons}
  \item{point.pch}{a vector of length 2 with the point characters of the points that are inside and outside of the simulated polygons}
  \item{\dots}{further arguments passed to \code{\link{plot.bivrp}}}
  \item{x}{an object of class \code{bivrp}}
}
\details{
This approach relies on the same strategy used for producing half-normal plots with simulation envelopes. Given a vector of bivariate model diagnostics, the angle each point makes with the origin is calculated to order them. This can be fine-tuned using the logical arguments \code{closest.angle}, \code{angle.ref}, and \code{counter.clockwise}, see the Arguments section above.

Then, \code{sim} bivariate response variables are simulated from the fitted model, using the same model matrices, error distribution and fitted parameters, using the function defined as \code{simfun}. The model is refitted to each simulated sample, using the function defined as \code{fitfun}. Next, we obtain the same type of model diagnostics, using \code{diagfun}, again ordered the same way the original bivariate sample was. We have, for each bivariate diagnostic, \code{sim} simulated bivariate diagnostics forming the whole cloud of simulated diagnostics.

By default, we then obtain the convex hulls of each set of the $s$ sets of points and obtain a reduced polygon whose area is (\code{conf * 100})\% of the original convex hull's area, forming the simulated polygon. This is equivalent to passing the argument \code{reduce.polygon = "proportional"} to \code{\link{plot.bivrp}}. The argument \code{reduce.polygon = "bag"} can be used to obtain a (\code{conf * 100})\% bagplot as the simulated polygon instead of a convex hull. The points are then connected to the centroids of their respective simulated polygons and, if they lie outside the polygons, they are drawn in red. For the final display, the polygons are erased so as to ease visualization.

There is no automatic implementation of a bivariate model in this function, and hence users must provide three functions for \code{bivrp}. The first function, \code{diagfun}, must extract the desired model diagnostics from a model fit object. The second function, \code{simfun}, must return the response variable, simulated using the same error distributions and estimated parameters from the fitted model. The third and final function, \code{fitfun}, must return a fitted model object. See the Examples section.

This function produces a plot by passing the computed object to \code{\link{plot.bivrp}}. The \code{print} method returns a \code{data.frame} containing all ordered simulated bivariate diagnostics. 
}
\value{
The function returns an object of class "bivrp", which is a list containing the following components:

\item{reslist.ord }{list of ordered diagnostics from model refitting to each simulated dataset}
\item{res.original.ord }{original model diagnostics}
\item{res1 }{diagnostics from variable 1}
\item{res2 }{diagnostics from variable 2}
\item{res.original1 }{original model diagnostics for variable 1}
\item{res.original2 }{original model diagnostics for variable 2}
\item{conf }{confidence level of the simulated polygons}
}

\author{
Rafael A. Moral <rafael.deandrademoral@mu.ie> and John Hinde
}

\seealso{
\code{\link{plot.bivrp}}
}

\examples{
## simulating a bivariate normal response variable

require(mvtnorm)

n <- 40
beta1 <- c(2, .4)
beta2 <- c(.2, .2)
x <- seq(1, 10, length = n)
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
bivnormfit <- function(Y, X, covariance) {
  n <- nrow(X)
  p <- ncol(X)
  y <- cbind(Y[1:n],Y[(n+1):(2*n)])
  XtXinv <- solve(crossprod(X, X))
  beta.hat <- XtXinv \%*\% crossprod(X, y)
  mu.hat <- X\%*\%beta.hat
  sigma.hat <- 1/n * t(y - mu.hat) \%*\% (y - mu.hat)
  if(!covariance) sigma.hat <- diag(diag(sigma.hat))
  cov.betas <- sigma.hat \%x\% XtXinv
  se.s1 <- sqrt(2*sigma.hat[1]^2/(n-p+1))
  se.s2 <- sqrt(2*sigma.hat[4]^2/(n-p+1))
  if(!covariance) se.s12 <- NA else {
    rho <- sigma.hat[2]/sqrt(sigma.hat[1]*sigma.hat[4])
    se.s12 <- sqrt((1+rho^2)*sigma.hat[1]*sigma.hat[4]/(n-p+1))
  }
  se.betas <- sqrt(diag(cov.betas))
  se.sigma <- c(se.s1, se.s2, se.s12)
  coefs <- c(beta.hat, sigma.hat[1], sigma.hat[4], sigma.hat[2])
  names(coefs) <- c("beta1.0", "beta1.1", "beta2.0", "beta2.1", "sig1", "sig2", "sig12")
  fitted <- c(mu.hat)
  resid <- Y - fitted
  Sig1 <- diag(rep(sigma.hat[1]), n)
  Sig2 <- diag(rep(sigma.hat[4]), n)
  Sig12 <- diag(rep(sigma.hat[2]), n)
  V <- rbind(cbind(Sig1, Sig12),
             cbind(Sig12, Sig2))
  llik <- dmvnorm(Y, c(mu.hat), V, log = TRUE)
  ret <- list("coefs" = coefs, "covariance" = covariance, "n" = n, 
              "X" = X, "fitted" = fitted, "resid" = resid, "loglik" = llik,
              "Y" = Y, "se" = c(se.betas, se.sigma))
  class(ret) <- "bivnormfit"
  return(ret)
}

## fitting bivariate models with and without estimating covariance
fit0 <- bivnormfit(Y, X, covariance=FALSE)
fit1 <- bivnormfit(Y, X, covariance=TRUE)
## likelihood-ratio test
2*(fit0$loglik - fit1$loglik)
pchisq(54.24, 1, lower=FALSE)

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
  Y <- as.numeric(rmvnorm(1, c(mu1, mu2), V))
  return(list(Y[1:n], Y[(n+1):(2*n)], "X" = obj$X, 
              "covariance" = obj$covariance))
}

## function for refitting the model to simulated data
ffun <- function(new.obj) {
  Ynew <- c(new.obj[[1]], new.obj[[2]])
  bivnormfit(Ynew, new.obj$X, new.obj$covariance)
}

## Bivariate residual plot for model 1 (without estimating covariance)
plot1 <- bivrp(fit0, diagfun=dfun, simfun=sfun, fitfun=ffun, verb=TRUE)
## without polygon area reduction
plot(plot1, conf=1)
## drawing polygons
plot(plot1, add.polygon=TRUE)
## without ordering
plot(plot1, theta.sort=FALSE, kernel=TRUE, add.dplots=TRUE, superpose=TRUE)

## Bivariate residual plot for model 2 (estimating covariance)
plot2 <- bivrp(fit1, diagfun=dfun, simfun=sfun, fitfun=ffun, verb=TRUE)
## without polygon area reduction
plot(plot2, conf=1)
## drawing polygons
plot(plot2, add.polygon=TRUE, conf=1)
## without ordering
plot(plot2, theta.sort=FALSE, kernel=TRUE, add.dplots=TRUE, superpose=TRUE)
}