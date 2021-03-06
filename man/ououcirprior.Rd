% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ouxy.r
\name{ououcirprior}
\alias{ououcirprior}
\title{Draw prior samples for OUOUCIR model}
\usage{
ououcirprior(
  prior.model.params = prior.model.params,
  prior.reg.params = prior.reg.params
)
}
\arguments{
\item{prior.model.params}{A vectors of hyper parameters for model parameter containing (alpha.y.min, alpha.y.max) for alpha.y, (alpha.x.min, alpha.x.max) for alpha.x, (theta.x.min, theta.x.max) for theta.x, (sigmasq.x.min, sigmasq.x.max) for sigmasq.x, (alpha.tau.min, alpha.tau.max) for alpha.tau, (theta.tau.min, theta.tau.max) for theta_tau, (sigmasq.tau.min, sigmasq.tau.max) for rate parameter of tau}

\item{prior.reg.params}{A vector of hyper paramter for regression parameters. (b0.min,b0.max) for b0, (b1.min,b1.max) for b1,  (b2.min,b2.max) for b2}
}
\value{
Returns the  samples of model parameters and regression parameter
}
\description{
Simulate sample for parameters in OUOUCIR model given a set of hyper parameters
}
\details{
The function requires user to input hyper parameters for \eqn{\alpha_y,\alpha_x,\theta_x, \sigma^2_x,\alpha_\tau, \theta_\tau,\sigma^2_\tau} and hyper parameters for regression parameters \eqn{b_0,b_1,b_2} from uniform distribution with its minimum and maximum values.
}
\examples{

prior.model.params<-c(0,3,0,3,-2,2,0,1,0,3,0,2,0,1.5,0,1)


names(prior.model.params)<-c(

"alpha.y.min","alpha.y.max","alpha.x.min","alpha.x.max",

"theta.x.min","theta.x.max","sigmasq.x.min","sigmasq.x.max",

"alpha.tau.min","alpha.tau.max","theta.tau.min","theta.tau.max",

"sigmasq.tau.min","sigmasq.tau.max")
prior.reg.params<-c(-3, 3, -3, 3, -3, 3)
names(prior.reg.params)<-c("b0.min", "b0.max", "b1.min", "b1.max", "b2.min", "b2.max")
ououcirprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)


}
