# library(TreeSim)
# library(EasyABC)
# library(coda)
# library(Sim.DiffProc)
# library(MCMCpack)
# library(ape)
# library(abc)
# library(nlme)
# library(adephylo)
# library(maps)
# library(phylobase)
# library(phytools)
# if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
#########
### Model
#########

# Roxygen document formulation http://r-pkgs.had.co.nz/man.html


### OUBMBM

#' Simulate traits under OUBMBM model given a set of parameters
#'
#' Simulate traits under OUBMBM model given a set of model parameters, regression parameters, tree and ancestral values.
#'
#' The model requires user to input  model parameters \eqn{\alpha_y, \sigma^2_x,\tau} and regression parameters \eqn{b_0,b_1,b_2}, a tree object, trait dataset (one response, two covariate), ancestral values(root) which is estimated by BM or OU model from \code{\link[geiger:fitContinuous]{geiger}}. Then the algorithm starts from the root and apply the tree traversal algorithm to simulate trait on each node according the OUOUBM dynamics.
#'
#' @param model.params  A vector of model parameters including alpha.y: force parameter of response trait, sigmasq.x: rate parameter of covariate and tau: rate parameter of response trait
#' @param reg.params A vector of regression paramters including  b0: intercept parameter of regression, b1: slope parameter of first covariate,  b2: slope paramter of the second covariate
#' @param root A vector of numerical values for root of species estimated from \code{OUprior} estimated from \code{OUprior}
#' @param tree  An ape: tree object stored in phylo format
#'
#' @return Returns the trait vectors \eqn{Y=(y_1,y_2,\cdots,y_n)', X_1=(x_{1,1},x_{1.2},\cdots,x_{1,n})',X_2=(x_{2,1},x_{2,2},\cdots,x_{2,n})'} simulated from the model.
#'
#'
#'
#' @references
#'  \enumerate{
#'   \item Jhwueng, D-C. (2019) Statistical modeling for adaptive trait evolution. Under review.
#'   \item Jhwueng, D-C., and Vasileios Maroulas. "Adaptive trait evolution in random environment." Journal of Applied Statistics 43.12 (2016): 2310-2324.
#'   \item Hansen, Thomas F., Jason Pienaar, and Steven Hecht Orzack. "A comparative method for studying adaptation to a randomly evolving environment." Evolution: International Journal of Organic Evolution 62.8 (2008): 1965-1977.
#' }
#'
#'
#' @examples
#' library(ape)
#' tree<-rcoal(5)
#' tree<-reorder(tree,"postorder")
#' root<-list(y.ou.sigmsq=1 ,y.ou.root=0, x1.ou.root=0, x2.ou.root=0,x1.bm.root=0, x2.bm.root=0)
#' model.params<-c(0.5,1,0.2)
#' names(model.params)<-c("alpha.y","sigmasq.x","tau")
#' reg.params<-c(0,1,1)
#' names(reg.params)<-c("b0","b1","b2")
#' oubmbmmodel(model.params,reg.params,root=root,tree=tree)
#'
#' @import stats ape coda Sim.DiffProc MCMCpack abc phytools nlme TreeSim adephylo maps geiger EasyABC
#'
#' @import utils
#'
#' @export
oubmbmmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  pos <- 1
  envir = as.environment(pos)
  assign("alpha.y", alpha.y, envir = envir)
  sigmasq.x<-model.params[2]
  tau<-model.params[3]
  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  n<-Ntip(tree)
  x1nodestates<-array(NA,c(2*n-1))
  x1nodestates[n+1]<- root$x1.bm.root
  x2nodestates<-array(NA,c(2*n-1))
  x2nodestates[n+1]<- root$x2.bm.root
  optimnodestates<-array(NA,c(2*n-1))
  optimnodestates[n+1]<-  b0 + b1*root$x1.bm.root + b2*root$x2.bm.root
  ynodestates<-array(NA,c(2*n-1))
  ynodestates[n+1]<-root$y.ou.root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length
  for(index in N:1){
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1nodestates[anc[index]],sd= sqrt(sigmasq.x*treelength[index]))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2nodestates[anc[index]],sd= sqrt(sigmasq.x*treelength[index]))
    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    sigmasq.theta<- b1^2*sigmasq.x + b2^2*sigmasq.x

    ###
    INT1var<-sigmasq.theta*(exp(2*alpha.y*treelength[index])-1)/(2*alpha.y)
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(INT1var))

    ###
    #fexpr<-expression(exp(alpha.y*t)*w)
    fexpr<-expression(exp(t/2)*w) #assume alpha.y = 0.5
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-tau*exp(-alpha.y*treelength[index])*median(res$X)

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
  }
  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
}

#' Draw prior samples for OUBMBM model
#'
#' Simulate sample for parameters in OUBMBM model given a set of hyper parameters
#'
#' The function requires user to input hyper parameters for \eqn{\alpha_y, \sigma^2_x,\tau} and hyper parameters for regression parameters \eqn{b_0,b_1,b_2} from uniform distribution with its minimum and maximum values.
#'
#' @param prior.model.params  A vectors of hyper parameters for model parameter containing (alpha.y.min, alpha.y.max) for alpha.y, (sigmasq.x.min, sigmasq.x.max) for sigmasq.x, (tau.y.min, tau.max) for rate parameter of tau
#' @param prior.reg.params A vector of hyper paramter for regression parameters. (b0.min,b0.max) for b0, (b1.min,b1.max) for b1,  (b2.min,b2.max) for b2
#'
#' @return Returns the  samples of model parameters and regression parameter
#'
#' @export
#'
#' @examples
#'
#' prior.model.params<-c(0,3,0,3,0,1)
#'
#' names(prior.model.params)<-c("alpha.y.min","alpha.y.max",
#'
#' "tau.min","tau.max","sigmasq.x.min","sigmasq.x.max")
#'
#' prior.reg.params<-c(-3, 3, -3, 3, -3, 3)
#' names(prior.reg.params)<-c("b0.min", "b0.max", "b1.min", "b1.max", "b2.min", "b2.max")
#' oubmbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
#'
#'
oubmbmprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
  alpha.y.min <-prior.model.params["alpha.y.min"]
  alpha.y.max <-prior.model.params["alpha.y.max"]
  sigmasq.x.min <-prior.model.params["sigmasq.x.min"]
  sigmasq.x.max <-prior.model.params["sigmasq.x.max"]
  tau.min <- prior.model.params["tau.min"]
  tau.max <- prior.model.params["tau.max"]

  alpha.y<-runif(n=1,min=alpha.y.min,max=alpha.y.max)
  sigmasq.x<-runif(n=1,min=sigmasq.x.min,max=sigmasq.x.max)
  tau<-runif(n=1,min=tau.min,max=tau.max)

  b0.min<-prior.reg.params[1]
  b0.max<-prior.reg.params[2]
  b1.min<-prior.reg.params[3]
  b1.max<-prior.reg.params[4]
  b2.min<-prior.reg.params[5]
  b2.max<-prior.reg.params[6]
  b0<-runif(n=1, min=b0.min, max=b0.max)
  b1<-runif(n=1, min=b1.min, max=b1.max)
  b2<-runif(n=1, min=b2.min, max=b2.max)

  model.params<-c(alpha.y, sigmasq.x, tau)
  reg.params<- c(b0, b1, b2)

  return(list(model.params=model.params, reg.params=reg.params))
}


### OUOUBM
#' Simulate traits under OUOUBM model given a set of parameters
#'
#' Simulate traits under OUOUBM model given a set of model parameters, regression parameters, tree and ancestral values.
#'
#' The model requires user to input  model parameters \eqn{\alpha_y,\alpha_x,\theta_x, \sigma^2_x,\tau} and regression parameters \eqn{b_0,b_1,b_2}, a tree object, trait dataset (one response, two covariate), ancestral values(root) which is estimated by BM or OU model from \code{\link[geiger:fitContinuous]{geiger}}. Then the algorithm starts from the root and apply the tree traversal algorithm to simulate trait on each node according the OUOUBM dynamics.
#'
#'
#'
#'
#' @param model.params  A vector of model parameters including alpha.y: force parameter of response trait, alpha.x: force parameter of covariate, theta.x: optimum parameter of covariate, sigmasq.x: rate parameter of covariate and tau rate parameter of response trait
#'
#' @param reg.params A vector of regression paramters including  b0: intercept parameter of regression, b1: slope parameter of first covariate,  b2: slope paramter of the second covariate
#' @param root A vector of numerical values for root of species estimated from \code{OUprior}
#' @param tree  An ape: tree object stored in phylo format
#'
#' @return Returns the trait vectors \eqn{Y=(y_1,y_2,\cdots,y_n)', X_1=(x_{1,1},x_{1.2},\cdots,x_{1,n})',X_2=(x_{2,1},x_{2,2},\cdots,x_{2,n})'} simulated from the model.
#'
#' @export
#'
#' @references
#'  \enumerate{
#'   \item Jhwueng, D-C. (2019) Statistical modeling for adaptive trait evolution. Under review.
#'   \item Jhwueng, D-C., and Vasileios Maroulas. "Adaptive trait evolution in random environment." Journal of Applied Statistics 43.12 (2016): 2310-2324.
#'   \item Hansen, Thomas F., Jason Pienaar, and Steven Hecht Orzack. "A comparative method for studying adaptation to a randomly evolving environment." Evolution: International Journal of Organic Evolution 62.8 (2008): 1965-1977.
#' }
#'
#'
#' @examples
#'library(ape)
#'tree<-rcoal(5)
#'tree<-reorder(tree,"postorder")
#'root<-list(y.ou.sigmsq=1 ,y.ou.root=0, x1.ou.root=0, x2.ou.root=0,x1.bm.root=0, x2.bm.root=0)
#'model.params<-c(0.5,0.25,0,1,0.2)
#'names(model.params)<-c("alpha.y","alpha.x","theta.x","sigmasq.x","tau")
#'reg.params<-c(0,1,1)
#'names(reg.params)<-c("b0","b1","b2")
#'ououbmmodel(model.params,reg.params,root=root,tree=tree)
#'
#'
ououbmmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]

  pos <- 1
  envir = as.environment(pos)
  assign("alpha.y", alpha.y, envir = envir)

  alpha.x<-model.params[2]
  theta.x<-model.params[3]
  sigmasq.x<-model.params[4]
  tau<-model.params[5]

  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  n<-Ntip(tree)
  x1nodestates<-array(NA,c(2*n-1))
  x1nodestates[n+1]<- root$x1.ou.root
  x2nodestates<-array(NA,c(2*n-1))
  x2nodestates[n+1]<- root$x2.ou.root
  optimnodestates<-array(NA,c(2*n-1))
  optimnodestates[n+1]<- b0 + b1*root$x1.ou.root + b2*root$x2.ou.root
  ynodestates<-array(NA,c(2*n-1))
  ynodestates[n+1]<-root$y.ou.root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length

  for(index in N:1){
    x1ou.mean<-x1nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + theta.x*(1-exp(-alpha.x*treelength[index]))
    x1ou.sd<-sqrt((sigmasq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    if(x1ou.sd<1e-10){x1ou.sd<-sigmasq.x*treelength[index]}
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd=x1ou.sd)

    x2ou.mean<-x2nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + theta.x*(1-exp(-alpha.x*treelength[index]))
    x2ou.sd<-sqrt((sigmasq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    if(x2ou.sd<1e-10){x2ou.sd<-sigmasq.x*treelength[index]}
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd=x2ou.sd)

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    sigmasq.theta<- b1^2*sigmasq.x + b2^2*sigmasq.x

    ###
    theta0.y<-optimnodestates[des[index]]
    A<-(alpha.y*theta0.y/ (alpha.y-alpha.x)) *(exp((alpha.y-alpha.x)*treelength[index]) -1)
    tilde.theta.y<- optimnodestates[des[index]]
    B<- tilde.theta.y*(exp(alpha.y*treelength[index])-1) - (alpha.y*tilde.theta.y/(alpha.y-alpha.x))*(exp((alpha.y-alpha.x)*treelength[index]) -1)
    vs<-sigmasq.theta*alpha.y^2*exp(2*alpha.y*treelength[index])*(1-exp(-2*alpha.x*treelength[index]))/ (2*alpha.x)
    C<-pnorm(treelength[index],mean=0,sd=sqrt(vs))- 0.5#integral starts from 0
    INTtime<- A+B+C
    INT1<-exp(-alpha.y*treelength[index])*INTtime

    ###
    fexpr<-expression(exp(alpha.y*t)*w)
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-tau*exp(-alpha.y*treelength[index])*median(res$X)

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
    }

  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  }


#' Draw prior samples for OUOUBM model
#'
#' Simulate sample for parameters in OUOUBM model given a set of hyper parameters
#'
#' The function requires user to input hyper parameters for \eqn{\alpha_y,\alpha_x,\theta_x, \sigma^2_x,\tau} and hyper parameters for regression parameters \eqn{b_0,b_1,b_2} from uniform distribution with its minimum and maximum values.
#'
#' @param prior.model.params  A vectors of hyper parameters for model parameter containing (alpha.y.min, alpha.y.max) for alpha.y, (alpha.x.min, alpha.x.max) for alpha.x,(theta.y.min, theta.y.max) for theta.y, (sigmasq.x.min, sigmasq.x.max) for sigmasq.x, (tau.y.min, tau.max) for rate parameter of tau
#' @param prior.reg.params A vector of hyper paramter for regression parameters. (b0.min,b0.max) for b0, (b1.min,b1.max) for b1,  (b2.min,b2.max) for b2
#'
#' @return Returns the  samples of model parameters and regression parameter
#'
#' @export
#'
#' @examples
#'
#'prior.model.params<-c(0,3,0,3,-5,5,0,1,0,2)
#'names(prior.model.params)<-c(
#'
#'"alpha.y.min","alpha.y.max","alpha.x.min","alpha.x.max",
#'
#'"theta.x.min","theta.x.max","sigmasq.x.min","sigmasq.x.max",
#'
#'"tau.min","tau.max")
#'prior.reg.params<-c(-3, 3, -3, 3, -3, 3)
#'names(prior.reg.params)<-c("b0.min", "b0.max", "b1.min", "b1.max", "b2.min", "b2.max")
#'ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
#'
#'
  ououbmprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
    alpha.y.min <-prior.model.params["alpha.y.min"]
    alpha.y.max <-prior.model.params["alpha.y.max"]
    alpha.x.min <- prior.model.params["alpha.x.min"]
    alpha.x.max <- prior.model.params["alpha.x.max"]
    theta.x.min <-prior.model.params["theta.x.min"]
    theta.x.max <- prior.model.params["theta.x.max"]
    sigmasq.x.min <-prior.model.params["sigmasq.x.min"]
    sigmasq.x.max <- prior.model.params["sigmasq.x.max"]
    tau.min <- prior.model.params["tau.min"]
    tau.max <- prior.model.params["tau.max"]

    alpha.y<-runif(n=1,min=alpha.y.min,max=alpha.y.max)
    alpha.x<-runif(n=1,min=alpha.x.min,max=alpha.x.max)
    theta.x<-runif(n=1,min=theta.x.min,max=theta.x.max)
    sigmasq.x<-runif(n=1,min=sigmasq.x.min,max=sigmasq.x.max)
    tau<-runif(n=1,min=tau.min,max=tau.max)

    b0.min<-prior.reg.params[1]
    b0.max<-prior.reg.params[2]
    b1.min<-prior.reg.params[3]
    b1.max<-prior.reg.params[4]
    b2.min<-prior.reg.params[5]
    b2.max<-prior.reg.params[6]
    b0<-runif(n=1, min=b0.min, max=b0.max)
    b1<-runif(n=1, min=b1.min, max=b1.max)
    b2<-runif(n=1, min=b2.min, max=b2.max)

    model.params<-c(alpha.y, alpha.x, theta.x, sigmasq.x, tau)
    reg.params<- c(b0, b1, b2)

    return(list(model.params=model.params, reg.params=reg.params))
  }


### OUBMCIR
  #' Simulate traits under OUBMCIR model given a set of parameters
  #'
  #' Simulate traits under OUBMCIR model given a set of model parameters, regression parameters, tree and ancestral values.
  #'
  #' The model requires user to input  model parameters \eqn{\alpha_y,\sigma^2_x,\alpha_\tau,\theta_\tau,\sigma^2_\tau} and regression parameters \eqn{b_0,b_1,b_2}, a tree object, trait dataset (one response, two covariate), ancestral values(root) which is estimated by BM or OU model from \code{\link[geiger:fitContinuous]{geiger}}. Then the algorithm starts from the root and apply the tree traversal algorithm to simulate trait on each node according the OUOUBM dynamics.
  #'
  #' @param model.params  A vector of model parameters including alpha.y: force parameter of response trait, sigmasq.x: rate parameter of covariate, alpha.tau: force parameter of rate, theta.tau optimum parameter of rate, sigmasq.tau: rate parameter of rate
  #'
  #' @param reg.params A vector of regression paramters including  b0: intercept parameter of regression, b1: slope parameter of first covariate,  b2: slope paramter of the second covariate
  #' @param root A vector of numerical values for root of species estimated from \code{OUprior}
  #' @param tree  An ape: tree object stored in phylo format
  #'
  #' @return Returns the trait vectors \eqn{Y=(y_1,y_2,\cdots,y_n)', X_1=(x_{1,1},x_{1.2},\cdots,x_{1,n})',X_2=(x_{2,1},x_{2,2},\cdots,x_{2,n})'} simulated from the model.
  #'
  #' @export
  #'
  #' @references
  #'  \enumerate{
  #'   \item Jhwueng, D-C. (2019) Statistical modeling for adaptive trait evolution. Under review.
  #'   \item Jhwueng, D-C., and Vasileios Maroulas. "Adaptive trait evolution in random environment." Journal of Applied Statistics 43.12 (2016): 2310-2324.
  #'   \item Hansen, Thomas F., Jason Pienaar, and Steven Hecht Orzack. "A comparative method for studying adaptation to a randomly evolving environment." Evolution: International Journal of Organic Evolution 62.8 (2008): 1965-1977.
  #' }
  #'
  #' @examples
  #' library(ape)
  #' tree<-rcoal(5)
  #' tree<-reorder(tree,"postorder")
  #' root<-list(y.ou.sigmsq=1 ,y.ou.root=0, x1.ou.root=0, x2.ou.root=0,x1.bm.root=0, x2.bm.root=0)
  #' model.params<-c(0.5,0.25,0.3,1,0.2)
  #' names(model.params)<-c("alpha.y","alpha.x","theta.x","sigmasq.x","sigmasq.tau")
  #' reg.params<-c(0,1,1)
  #' names(reg.params)<-c("b0","b1","b2")
  #' oubmcirmodel(model.params,reg.params,root=root,tree=tree)
  #'
  oubmcirmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  sigmasq.x<-model.params[2]
  alpha.tau<-model.params[3]
  theta.tau<-model.params[4]
  sigmasq.tau<-model.params[5]
  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  n<-Ntip(tree)
  x1nodestates<-array(NA,c(2*n-1))
  x1nodestates[n+1]<- root$x1.bm.root
  x2nodestates<-array(NA,c(2*n-1))
  x2nodestates[n+1]<- root$x2.bm.root
  optimnodestates<-array(NA,c(2*n-1))
  optimnodestates[n+1]<-  b0 + b1*root$x1.bm.root + b2*root$x2.bm.root
  sigmasqnodestates<-array(NA,c(2*n-1))
  sigmasqnodestates[n+1]<- root$y.ou.sigmsq
  ynodestates<-array(NA,c(2*n-1))
  ynodestates[n+1]<-root$y.ou.root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length
  for(index in N:1){

    x1nodestates[des[index]]<-rnorm(n=1,mean=x1nodestates[anc[index]],sd= sqrt(sigmasq.x*treelength[index]))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2nodestates[anc[index]],sd= sqrt(sigmasq.x*treelength[index]))

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    sigmasq.theta<- b1^2*sigmasq.x + b2^2*sigmasq.x

    c<- sigmasq.tau*(1-exp(-alpha.tau*treelength[index]))/(4*alpha.tau)
    k<- (4*theta.tau*alpha.tau)/sigmasq.tau
    tau0<-sigmasqnodestates[anc[index]]
    lambda<-4*alpha.tau*exp(-alpha.tau*treelength[index])
    lambda<-lambda/(sigmasq.tau*(1-exp(-alpha.tau*treelength[index])))*tau0
    tmp = rchisq(n=1, df=k, ncp = lambda)
    sig_u <- c*tmp
    sigmasqnodestates[des[index]]<-sig_u

    ###
    INT1var<-sigmasq.theta*(exp(2*alpha.y*treelength[index])-1)/(2*alpha.y)
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(INT1var))

    ###
    a.var<-(1-exp(-2*alpha.y*treelength[index]))/(2*alpha.y)
    a <- rnorm(n=1, mean=0, sd=sqrt(a.var))
    b.var<- (sigmasqnodestates[anc[index]]-theta.tau)^2/(2*(alpha.y-alpha.tau))
    b.var<-b.var*(exp(-2*alpha.tau*treelength[index])-exp(-2*alpha.y*treelength[index]))
    b <- rnorm(n=1, mean=0, sd=sqrt(b.var))

    n_t <- 1
    n_s <- 1
    outer.int.sum=0
    for(outer.index in 1:n_t){
      inner.int.sum = 0
      for(inner.index in 1:n_s){
        c<- sigmasq.tau*(1-exp(-alpha.tau*(inner.index/n_s)))/(4*alpha.tau)
        k<- (4*theta.tau*alpha.tau)/sigmasq.tau
        lambda<- 4*alpha.tau*exp(-alpha.tau*(inner.index/n_s))
        tau0<-sigmasqnodestates[anc[index]]
        lambda<-lambda/(sigmasq.tau*(1-exp(-alpha.tau*(inner.index/n_s))))*tau0
        tmp = rchisq(n=1, df=k, ncp = lambda)
        sig_u <- c*tmp
        inner.int.sum  <-  inner.int.sum + exp(alpha.tau*(inner.index/n_s))*rnorm(n=1,mean=0, sd=sqrt(1/n_s))*sqrt(sig_u)
      }
      outer.int.sum <- outer.int.sum + exp(-(alpha.y+alpha.tau)*(outer.index/n_t))*inner.int.sum*rnorm(n=1,mean=0, sd=sqrt(1/n_t))
    }
    c <- sqrt(sigmasq.tau)*outer.int.sum
    INT2 <- (a + b + c)
    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
  }
  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  }

  #' Draw prior samples for OUBMCIR model
  #'
  #' Simulate sample for parameters in OUBMCIR model given a set of hyper parameters
  #'
  #' The function requires user to input hyper parameters for \eqn{\alpha_y,\sigma^2_x,\alpha_\tau, \theta_\tau,\sigma^2_\tau} and hyper parameters for regression parameters \eqn{b_0,b_1,b_2} from uniform distribution with its minimum and maximum values.
  #'
  #' @param prior.model.params  A vectors of hyper parameters for model parameter containing (alpha.y.min, alpha.y.max) for alpha.y, (sigmasq.x.min, sigmasq.x.max) for sigmasq.x, (alpha.tau.min, alpha.tau.max) for alpha.tau, (theta.tau.min, theta.tau.max) for theta_tau, (sigmasq.tau.min, sigmasq.tau.max) for rate parameter of tau
  #' @param prior.reg.params A vector of hyper paramter for regression parameters. (b0.min,b0.max) for b0, (b1.min,b1.max) for b1,  (b2.min,b2.max) for b2
  #'
  #' @return Returns the  samples of model parameters and regression parameter
  #'
  #' @export
  #'
  #' @examples
  #'
  #'prior.model.params<-c(0,3,0,3,0,1,0,3,0,2,0,1.5)
  #'
  #'
  #'names(prior.model.params)<-c(
  #'
  #'"alpha.y.min","alpha.y.max","sigmasq.x.min","sigmasq.x.max",
  #'
  #'"alpha.tau.min","alpha.tau.max","theta.tau.min","theta.tau.max",
  #'
  #'"sigmasq.tau.min","sigmasq.tau.max")
  #'prior.reg.params<-c(-3, 3, -3, 3, -3, 3)
  #'names(prior.reg.params)<-c("b0.min", "b0.max", "b1.min", "b1.max", "b2.min", "b2.max")
  #'oubmcirprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
  #'
  #'
  oubmcirprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
    alpha.y.min <-prior.model.params["alpha.y.min"]
    alpha.y.max <-prior.model.params["alpha.y.max"]
    sigmasq.x.min <-prior.model.params["sigmasq.x.min"]
    sigmasq.x.max <- prior.model.params["sigmasq.x.max"]
    alpha.tau.min <- prior.model.params["alpha.tau.min"]
    alpha.tau.max <- prior.model.params["alpha.tau.max"]
    theta.tau.min <- prior.model.params["theta.tau.min"]
    theta.tau.max <- prior.model.params["theta.tau.max"]
    sigmasq.tau.min <- prior.model.params["sigmasq.tau.min"]
    sigmasq.tau.max <- prior.model.params["sigmasq.tau.max"]

    alpha.y<-runif(n=1,min=alpha.y.min,max=alpha.y.max)
    sigmasq.x<-runif(n=1,min=sigmasq.x.min,max=sigmasq.x.max)
    alpha.tau <- runif(n=1,min=alpha.tau.min,max=alpha.tau.max)
    theta.tau <- runif(n=1,min=theta.tau.min,max=theta.tau.max)
    sigmasq.tau<-runif(n=1,min=sigmasq.tau.min,max=sigmasq.tau.max)

    b0.min<-prior.reg.params[1]
    b0.max<-prior.reg.params[2]
    b1.min<-prior.reg.params[3]
    b1.max<-prior.reg.params[4]
    b2.min<-prior.reg.params[5]
    b2.max<-prior.reg.params[6]
    b0<-runif(n=1, min=b0.min, max=b0.max)
    b1<-runif(n=1, min=b1.min, max=b1.max)
    b2<-runif(n=1, min=b2.min, max=b2.max)

    model.params<-c(alpha.y, sigmasq.x, alpha.tau, theta.tau, sigmasq.tau)
    reg.params<- c(b0, b1, b2)

    return(list(model.params=model.params, reg.params=reg.params))
  }


### OUOUCIR
  #' Simulate traits under OUBMCIR model given a set of parameters
  #'
  #' Simulate traits under OUBMCIR model given a set of model parameters, regression parameters, tree and ancestral values.
  #'
  #' The model requires user to input  model parameters \eqn{\alpha_y,\alpha_x,\theta_x,\sigma^2_x,\alpha_\tau,\theta_\tau,\sigma^2_\tau} and regression parameters \eqn{b_0,b_1,b_2}, a tree object, trait dataset (one response, two covariate), ancestral values(root) which is estimated by BM or OU model from \code{\link[geiger:fitContinuous]{geiger}}. Then the algorithm starts from the root and apply the tree traversal algorithm to simulate trait on each node according the OUOUBM dynamics.
  #'
  #' @param model.params  A vector of model parameters including alpha.y: force parameter of response trait, alpha.x: force parameter of covariate, theta.x" optimum parameter of covariate, sigmasq.x: rate parameter of covariate, alpha.tau: force parameter of rate, theta.tau optimum parameter of rate, sigmasq.tau: rate parameter of rate
  #'
  #' @param reg.params A vector of regression paramters including  b0: intercept parameter of regression, b1: slope parameter of first covariate,  b2: slope paramter of the second covariate
  #' @param root A vector of numerical values for root of species estimated from \code{OUprior}
  #' @param tree  An ape: tree object stored in phylo format
  #'
  #'
  #' @return Returns the trait vectors \eqn{Y=(y_1,y_2,\cdots,y_n)', X_1=(x_{1,1},x_{1.2},\cdots,x_{1,n})',X_2=(x_{2,1},x_{2,2},\cdots,x_{2,n})'} simulated from the model.
  #'
  #' @export
  #'
  #' @references
  #'  \enumerate{
  #'   \item Jhwueng, D-C. (2019) Statistical modeling for adaptive trait evolution. Under review.
  #'   \item Jhwueng, D-C., and Vasileios Maroulas. "Adaptive trait evolution in random environment." Journal of Applied Statistics 43.12 (2016): 2310-2324.
  #'   \item Hansen, Thomas F., Jason Pienaar, and Steven Hecht Orzack. "A comparative method for studying adaptation to a randomly evolving environment." Evolution: International Journal of Organic Evolution 62.8 (2008): 1965-1977.
  #' }
  #'
  #' @examples
  #' library(ape)
  #' tree<-rcoal(5)
  #' tree<-reorder(tree,"postorder")
  #' root<-list(y.ou.sigmsq=1 ,y.ou.root=0, x1.ou.root=0, x2.ou.root=0,x1.bm.root=0, x2.bm.root=0)
  #' model.params<-c(0.5,0.25,0,1,0.125,0.15,0.2)
  #' names(model.params)<-c("alpha.y","alpha.x","theta.x","sigmasq.x"
  #'
  #',"alpha.tau","theta.tau","sigmasq.tau")
  #' reg.params<-c(0,1,1)
  #' names(reg.params)<-c("b0","b1","b2")
  #' ououcirmodel(model.params,reg.params,root=root,tree=tree)
  #'
  ououcirmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  alpha.x<-model.params[2]
  theta.x<-model.params[3]
  sigmasq.x<-model.params[4]
  alpha.tau<-model.params[5]
  theta.tau<-model.params[6]
  sigmasq.tau<-model.params[7]
  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  n<-Ntip(tree)
  x1nodestates<-array(NA,c(2*n-1))
  x1nodestates[n+1]<-root$x1.ou.root
  x2nodestates<-array(NA,c(2*n-1))
  x2nodestates[n+1]<-root$x2.ou.root
  optimnodestates<-array(NA,c(2*n-1))
  optimnodestates[n+1]<-  b0 + b1*root$x1.ou.root + b2*root$x2.ou.root
  sigmasqnodestates<-array(NA,c(2*n-1))
  sigmasqnodestates[n+1]<- root$y.ou.sigmsq
  ynodestates<-array(NA,c(2*n-1))
  ynodestates[n+1]<-root$y.ou.root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length

  for(index in N:1){
    x1ou.mean<-x1nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + theta.x*(1-exp(-alpha.x*treelength[index]))
    x1ou.sd<-sqrt((sigmasq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    if(x1ou.sd<1e-10){x1ou.sd<-sigmasq.x*treelength[index]}
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd=x1ou.sd)

    x2ou.mean<-x2nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + theta.x*(1-exp(-alpha.x*treelength[index]))
    x2ou.sd<-sqrt((sigmasq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    if(x2ou.sd<1e-10){x2ou.sd<-sigmasq.x*treelength[index]}
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd=x2ou.sd)

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    sigmasq.theta<- b1^2*sigmasq.x + b2^2*sigmasq.x

    c<- sigmasq.tau*(1-exp(-alpha.tau*treelength[index]))/(4*alpha.tau)
    k<- (4*theta.tau*alpha.tau)/sigmasq.tau
    tau0<-sigmasqnodestates[anc[index]]
    lambda<-4*alpha.tau*exp(-alpha.tau*treelength[index])
    lambda<-lambda/(sigmasq.tau*(1-exp(-alpha.tau*treelength[index])))*tau0
    tmp = rchisq(n=1, df=k, ncp = lambda)
    sig_u <- c*tmp
    sigmasqnodestates[des[index]]<-sig_u

    #####
    theta0.y<-optimnodestates[des[index]]
    A<-(alpha.y*theta0.y/ (alpha.y-alpha.x)) *(exp((alpha.y-alpha.x)*treelength[index]) -1)
    tilde.theta.y<- optimnodestates[des[index]]
    B<- tilde.theta.y*(exp(alpha.y*treelength[index])-1) - (alpha.y*tilde.theta.y/(alpha.y-alpha.x))*(exp((alpha.y-alpha.x)*treelength[index]) -1)
    vs<-sigmasq.theta*alpha.y^2*exp(2*alpha.y*treelength[index])*(1-exp(-2*alpha.x*treelength[index]))/ (2*alpha.x)
    C<-pnorm(treelength[index],mean=0,sd=sqrt(vs))- 0.5#integral starts from 0
    INTtime<- A+B+C
    INT1<-exp(-alpha.y*treelength[index])*INTtime

    #####
    a.var<-(1-exp(-2*alpha.y*treelength[index]))/(2*alpha.y)
    a <- rnorm(n=1, mean=0, sd=sqrt(a.var))
    b.var<- (sigmasqnodestates[anc[index]]-theta.tau)^2/(2*(alpha.y-alpha.tau))
    b.var<-b.var*(exp(-2*alpha.tau*treelength[index])-exp(-2*alpha.y*treelength[index]))
    b <- rnorm(n=1, mean=0, sd=sqrt(b.var))

    n_t <- 1
    n_s <- 1
    outer.int.sum=0
    for(outer.index in 1:n_t){
      inner.int.sum = 0
      for(inner.index in 1:n_s){
        c<- sigmasq.tau*(1-exp(-alpha.tau*(inner.index/n_s)))/(4*alpha.tau)
        k<- (4*theta.tau*alpha.tau)/sigmasq.tau
        lambda<- 4*alpha.tau*exp(-alpha.tau*(inner.index/n_s))
        tau0<-sigmasqnodestates[anc[index]]
        lambda<-lambda/(sigmasq.tau*(1-exp(-alpha.tau*(inner.index/n_s))))*tau0
        tmp = rchisq(n=1, df=k, ncp = lambda)
        sig_u <- c*tmp
        inner.int.sum  <-  inner.int.sum + exp(alpha.tau*(inner.index/n_s))*rnorm(n=1,mean=0, sd=sqrt(1/n_s))*sqrt(sig_u)
      }
      outer.int.sum <- outer.int.sum + exp(-(alpha.y+alpha.tau)*(outer.index/n_t))*inner.int.sum*rnorm(n=1,mean=0, sd=sqrt(1/n_t))
    }
    c <- sqrt(sigmasq.tau)*outer.int.sum
    INT2 <- (a + b + c)
    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
    }
  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  }

  #' Draw prior samples for OUOUCIR model
  #'
  #' Simulate sample for parameters in OUOUCIR model given a set of hyper parameters
  #'
  #' The function requires user to input hyper parameters for \eqn{\alpha_y,\alpha_x,\theta_x, \sigma^2_x,\alpha_\tau, \theta_\tau,\sigma^2_\tau} and hyper parameters for regression parameters \eqn{b_0,b_1,b_2} from uniform distribution with its minimum and maximum values.
  #'
  #' @param prior.model.params  A vectors of hyper parameters for model parameter containing (alpha.y.min, alpha.y.max) for alpha.y, (alpha.x.min, alpha.x.max) for alpha.x, (theta.x.min, theta.x.max) for theta.x, (sigmasq.x.min, sigmasq.x.max) for sigmasq.x, (alpha.tau.min, alpha.tau.max) for alpha.tau, (theta.tau.min, theta.tau.max) for theta_tau, (sigmasq.tau.min, sigmasq.tau.max) for rate parameter of tau
  #' @param prior.reg.params A vector of hyper paramter for regression parameters. (b0.min,b0.max) for b0, (b1.min,b1.max) for b1,  (b2.min,b2.max) for b2
  #'
  #' @return Returns the  samples of model parameters and regression parameter
  #'
  #' @export
  #'
  #' @examples
  #'
  #'prior.model.params<-c(0,3,0,3,-2,2,0,1,0,3,0,2,0,1.5,0,1)
  #'
  #'
  #'names(prior.model.params)<-c(
  #'
  #'"alpha.y.min","alpha.y.max","alpha.x.min","alpha.x.max",
  #'
  #'"theta.x.min","theta.x.max","sigmasq.x.min","sigmasq.x.max",
  #'
  #'"alpha.tau.min","alpha.tau.max","theta.tau.min","theta.tau.max",
  #'
  #'"sigmasq.tau.min","sigmasq.tau.max")
  #'prior.reg.params<-c(-3, 3, -3, 3, -3, 3)
  #'names(prior.reg.params)<-c("b0.min", "b0.max", "b1.min", "b1.max", "b2.min", "b2.max")
  #'ououcirprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
  #'
  #'
  ououcirprior <- function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
    alpha.y.min <-prior.model.params["alpha.y.min"]
    alpha.y.max <-prior.model.params["alpha.y.max"]
    alpha.x.min <-prior.model.params["alpha.x.min"]
    alpha.x.max <-prior.model.params["alpha.x.max"]
    theta.x.min <-prior.model.params["theta.x.min"]
    theta.x.max <- prior.model.params["theta.x.max"]
    sigmasq.x.min <-prior.model.params["sigmasq.x.min"]
    sigmasq.x.max <- prior.model.params["sigmasq.x.max"]
    alpha.tau.min <- prior.model.params["alpha.tau.min"]
    alpha.tau.max <- prior.model.params["alpha.tau.max"]
    theta.tau.min <- prior.model.params["theta.tau.min"]
    theta.tau.max <- prior.model.params["theta.tau.max"]
    sigmasq.tau.min <- prior.model.params["sigmasq.tau.min"]
    sigmasq.tau.max <- prior.model.params["sigmasq.tau.max"]

    alpha.y<-runif(n=1,min=alpha.y.min,max=alpha.y.max)
    alpha.x<-runif(n=1,min=alpha.x.min,max=alpha.x.max)
    theta.x<-runif(n=1,min=theta.x.min,max=theta.x.max)
    sigmasq.x<-runif(n=1,min=sigmasq.x.min,max=sigmasq.x.max)
    alpha.tau <- runif(n=1,min=alpha.tau.min,max=alpha.tau.max)
    theta.tau <- runif(n=1,min=theta.tau.min,max=theta.tau.max)
    sigmasq.tau<-runif(n=1,min=sigmasq.tau.min,max=sigmasq.tau.max)

    b0.min<-prior.reg.params[1]
    b0.max<-prior.reg.params[2]
    b1.min<-prior.reg.params[3]
    b1.max<-prior.reg.params[4]
    b2.min<-prior.reg.params[5]
    b2.max<-prior.reg.params[6]
    b0<-runif(n=1, min=b0.min, max=b0.max)
    b1<-runif(n=1, min=b1.min, max=b1.max)
    b2<-runif(n=1, min=b2.min, max=b2.max)

    model.params<-c(alpha.y, alpha.x, theta.x, sigmasq.x, alpha.tau, theta.tau, sigmasq.tau)
    reg.params<- c(b0, b1, b2)

    return(list(model.params=model.params, reg.params=reg.params))

  }


######
###Empirical Analysis
#######

#' Fit OU model for univariate data
#'
#' Fit OU model given tree and trait
#'
#' Parameter estimates \eqn{\alpha, \theta, \sigma^2} are estimated by BM (when \eqn{\alpha = 0}) or OU model from \code{\link[geiger:fitContinuous]{geiger}} for the next step analysis with function \code{HyperParam} to get the reasonable range of the hyper parameter as well as the ancestral value.
#'
#' @param tree  An ape: tree object stored in phylo format
#' @param trait a univaraite trait
#' @param model specified model preset "OU".
#'
#'
#' @return MLE parameter estimates \eqn{\alpha, \theta, \sigma^2}.
#'
#' @export
#'
#' @references
#'  \enumerate{
#'   \item Jhwueng, D-C. (2019) Statistical modeling for adaptive trait evolution. Under review.
#'   \item Jhwueng, D-C., and Vasileios Maroulas. "Adaptive trait evolution in random environment." Journal of Applied Statistics 43.12 (2016): 2310-2324.
#'   \item Hansen, Thomas F., Jason Pienaar, and Steven Hecht Orzack. "A comparative method for studying adaptation to a randomly evolving environment." Evolution: International Journal of Organic Evolution 62.8 (2008): 1965-1977.
#' }
#'
#' @examples
#'
#'\donttest{
#'library(ape)
#'tree<-rcoal(3)
#'trait<-rnorm(3)
#'names(trait)<-tree$tip.label
#'model <- "OU"
#'OUprior(tree=tree,trait=trait,model=model)
#'}
#'

OUprior<-function(tree=tree,trait=trait,model=model){
  options(warn=-1)
  ou.result<-geiger::fitContinuous(tree,trait,model=model)
  alpha <- ou.result$opt$alpha#rexp(n=1, rate=1/alpha.prior.mean)
  if(model=="OU"){
      if(alpha < 1e-5){alpha <-0.001}
    }
  theta <- ou.result$opt$z0#rnorm(n=1, mean=theta.prior.mean, sd=theta.prior.sd)
  sigmasq<-ou.result$opt$sigsq#rexp(n=1, rate=1/sigmasq.prior.mean)
  return(list(alpha=alpha,theta=theta,sigmasq=sigmasq,root=theta))
}

  #' Summary statistics
  #'
  #' Calculate summary statistics given trait and tree
  #'
  #'  This function computes the 12 summary statistics using the trait and tree. \code{\link[ape:pic]{ape}} is used for computing the contrast trait from the difference between the species and its closet neighbor.
  #'  For Bloomberg K and Pagel Lambda, statsitics are computed using \code{\link[phytools:phylosig]{phytools}}.
  #'
  #' @param tree  An ape: tree object stored in phylo format
  #' @param trait A vector of numerical trait value
  #' @param ... relevent argument
  #'
  #' @return Twelve summary statistcs:  mean, sd, median, skewness, kurtosis from the raw data as well as from data with the difference between two closet neighbors, Bloomberg K and Pagel's lambda.
  #'
  #' @name sumstat
  #' @aliases sumstat
  #' @export
  #'
  #' @references
  #'  \enumerate{
  #'   \item Paradis E. & Schliep K. 2018. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526-528.
  #'   \item Blomberg, Simon P., Theodore Garland Jr, and Anthony R. Ives. "Testing for phylogenetic signal in comparative data: behavioral traits are more labile." Evolution 57.4 (2003): 717-745.
  #'   \item Pagel, Mark. "Inferring the historical patterns of biological evolution." Nature 401.6756 (1999): 877.
  #'   \item Revell, Liam J. "phytools: an R package for phylogenetic comparative biology (and other things)." Methods in Ecology and Evolution 3.2 (2012): 217-223.
  #' }
  #'
  #' @examples
  #' library(ape)
  #' tree<-rcoal(5)
  #' trait <- rnorm(5)
  #' names(trait)<-tree$tip.label
  #' sumstat(trait=trait,tree=tree)
  #'
  #'
  #'
 sumstat<-function(trait=trait,tree=tree, ...){
  names(trait)<-tree$tip.label
  pic.trait<-pic(x=trait,phy=tree)
  sumstatvalue<-c(mean(trait),sd(trait),median(trait),skewness(trait),kurtosis(trait),mean(pic.trait),sd(pic.trait),median(pic.trait),skewness(pic.trait),kurtosis(pic.trait),phylosig(tree,x=trait,method = "K",test=T)$K,phylosig(tree,x=trait,method = "lambda",test=T)$lambda)
  names(sumstatvalue)<-c("raw mean","raw sd","raw median","raw skewness","raw kurtosis","contrast mean","contrast sd","contrast median","contrast skewness","contrast kurtosis","K","lambda")
  return(sumstatvalue)
}


#' range for regression parameters
#'
#' Set up range for regression parameters
#'
#' An ordinary least square analysis is performed on regression \eqn{y \sim x_1 + x_2}. Parameter estimates \eqn{\hat{b}} and standard errors \eqn{sd(\hat{b})} are used to construct the bound using formula \eqn{\hat{b}\pm 3sd(\hat{b})}.
#'
#'
#' @param olssum summary statistics from ordinary least square performed by \code{\link[stats:lm]{lm}}.
#'
#'
#' @return A vectors of values containing the range of regression parameters
#'
#' @export
#'
#'
#' @examples
#' resptrait<-rnorm(10)
#' predtrait1<-rnorm(10)
#' predtrait2<-rnorm(10)
#' olssum <- base::summary(lm(resptrait~predtrait1+predtrait2))
#' regboundfcn(olssum=olssum)
#'
#'
regboundfcn<-function(olssum=olssum){
  bdd<-  c(olssum$coefficients[,1]-3*olssum$coefficients[,2],olssum$coefficients[,1]+3*olssum$coefficients[,2])
  bdd<-bdd[c(1,4,2,5,3,6)]
  return(bdd)
}


#WE USE UNIFORM PRIOR FOR EMPIRICAL ANALYSIS

#' The range of parameters
#'
#' Set up range for parameters for next step
#'
#' Function \code{\link{OUprior}} is called to compute the model estimate, then return the range for parameter estimate for next step analysis. The range is set to 3 times larger/smaller than the parameter estimates.
#' Function \code{\link{regboundfcn}} is called to get the bound of regression parameter.
#' The ancestral value (root) is computed for each traits in order to used for simulation in the four functions \code{\link{oubmbmTrait}},\code{\link{ououbmTrait}}, \code{\link{oubmcirTrait}} and \code{\link{ououcirTrait}}.
#'
#'
#' @param tree   An ape: tree object stored in phylo format
#' @param traitset a dataframe that contains 3 traits
#'
#' @return A list of vectors of sample of model parameters, regression parameter and ancestral values.
#'
#' @export
#'
#'
#' @examples
#'
#' ## using coral dataset (running time more > 5 sec)
#' \donttest{
#' data(coral)
#' tree<-coral$tree
#' traitset<-coral$traitset
#' HyperParam(tree=tree,traitset=traitset)
#'}

HyperParam<-function(tree=tree,traitset=traitset){
  resptrait<-traitset$resptrait
  predtrait1<-traitset$predtrait1
  predtrait2 <-traitset$predtrait2
  names(resptrait)<-tree$tip.label
  names(predtrait1)<-tree$tip.label
  names(predtrait2)<-tree$tip.label

  y.ou.result <- OUprior(tree=tree,trait=resptrait,model="OU")
  x1.ou.result<- OUprior(tree=tree,trait=predtrait1,model="OU")
  x2.ou.result <- OUprior(tree=tree,trait=predtrait2,model="OU")
  x1.bm.result<- OUprior(tree=tree,trait=predtrait1,model="BM")
  x2.bm.result <- OUprior(tree=tree,trait=predtrait2,model="BM")

  root<-list(y.ou.sigmsq=y.ou.result$sigmasq ,y.ou.root=y.ou.result$root, x1.ou.root=x1.ou.result$root, x2.ou.root=x2.ou.result$root,x1.bm.root=x1.bm.result$root, x2.bm.root=x2.bm.result$root )# Here we have root

  alpha.y.min <-  0
  alpha.y.max <-  3*y.ou.result$alpha

  alpha.x.min <-  0
  alpha.x.max <-  3*mean(x1.ou.result$alpha,x2.ou.result$alpha)#0.126#3*mean(x1.ou.result$alpha,x2.ou.result$alpha)

  theta.x.min <-  -3*abs(mean(x1.ou.result$theta,x2.ou.result$theta))#-0.01#-3*abs(mean(x1.ou.result$theta,x2.ou.result$theta))
  theta.x.max <-   3*abs(mean(x1.ou.result$theta,x2.ou.result$theta))#0.01#3*abs(mean(x1.ou.result$theta,x2.ou.result$theta))

  tau.min <-       0
  tau.max <-       3*mean(y.ou.result$sigma,x1.ou.result$sigma,x2.ou.result$sigma)

  sigmasq.x.min <- 0
  sigmasq.x.max <- 3*mean(x1.ou.result$sigmasq,x2.ou.result$sigmasq)

  alpha.tau.min <- 0
  alpha.tau.max <- 3*mean(y.ou.result$alpha, x1.ou.result$alpha,x2.ou.result$alpha)

  theta.tau.min <- 0
  theta.tau.max <- 3*mean(abs(y.ou.result$theta), abs(x1.ou.result$theta),abs(x2.ou.result$theta))

  sigmasq.tau.min <- 0
  sigmasq.tau.max <- 3*mean(y.ou.result$sigmasq,x1.ou.result$sigmasq,x2.ou.result$sigmasq)

  # For regression
  olssum <- summary(stats::lm(resptrait~predtrait1+predtrait2))
  regbound<-regboundfcn(olssum=olssum)

  b0.min=regbound[1]
  b0.max=regbound[2]
  b1.min=regbound[3]
  b1.max=regbound[4]
  b2.min=regbound[5]
  b2.max=regbound[6]

  prior.model.params<-c(alpha.y.min,alpha.y.max,alpha.x.min,alpha.x.max,theta.x.min,theta.x.max,tau.min,tau.max,sigmasq.x.min,sigmasq.x.max,alpha.tau.min,alpha.tau.max,theta.tau.min,theta.tau.max,sigmasq.tau.min,sigmasq.tau.max)
  names(prior.model.params)<-c("alpha.y.min","alpha.y.max","alpha.x.min","alpha.x.max","theta.x.min","theta.x.max","tau.min","tau.max","sigmasq.x.min","sigmasq.x.max","alpha.tau.min","alpha.tau.max","theta.tau.min","theta.tau.max","sigmasq.tau.min","sigmasq.tau.max")
  prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
  return(list(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params,root=root))
}


###############
### oubmbmTrait ###
##############
#' Parameter samples and summary statistics
#'
#' Draw sample for parameters, simulate trait and compute the summary statistics for OUBMBM model
#'
#' Given tree, trait sets, function \code{\link{HyperParam}} is called to yield the range of parameters, then function \code{\link{oubmbmprior}} is called to draw sample for parameter, then the function \code{\link{oubmbmmodel}} is applied to simulate traits through post order tree traversal algorithm, finally the summary statistics is computed by function \code{\link{sumstat}}.
#'
#' @param tree   An ape: tree object stored in phylo format
#' @param traitset a dataframe that contains 3 traits
#' @param sims number of trait replicate
#'
#' @return A list of vectors containing a dataframe of summary statistics, and a dataframe of parameter samples
#'
#' @export
#'
#'
#' @examples
#'
#' ## using bat dataset (running time more > 5 sec)
#' \donttest{
#' data(bat)
#' tree<-bat$tree
#' traitset<-bat$traitset
#' sims<-10
#' oubmbmTrait(tree=tree,traitset=traitset,sims=sims)
#'}
oubmbmTrait<-function(tree=tree,traitset=traitset,sims=sims){
  prior.params<-HyperParam(tree=tree,traitset=traitset)
  prior.model.params<-prior.params$prior.model.params
  prior.reg.params<-prior.params$prior.reg.params
  prior.params.oubmbm <- oubmbmprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  root<-prior.params$root
  sim.oubmbm.trait<-array(NA,c(dim(traitset)[1],3,sims))
  model.params.oubmbm<-array(NA,c(3,sims))
  rownames(model.params.oubmbm)<-c("alpha.y","sigmasq.x","tau")
  reg.params.oubmbm<-array(NA,c(3,sims))
  row.names(reg.params.oubmbm)<-c("b0", "b1", "b2")
  y.sumstat.oubmbm<-array(NA,c(12,sims))
  rownames(y.sumstat.oubmbm)<-c("y.trait.mean","y.trait.sd","y.trait.median","y.trait.skewness","y.trait.kurtosis","y.pic.trait.mean","y.pic.trait.sd","y.pic.trait.mediam","y.pic.trait.skewness","y.pic.trait.kurtosis","y.pic.trait.K","y.pic.trait.lambda")
  x1.sumstat.oubmbm<-array(NA,c(12,sims))
  rownames(x1.sumstat.oubmbm)<-c("x1.trait.mean","x1.trait.sd","x1.trait.median","x1.trait.skewness","x1.trait.kurtosis","x1.pic.trait.mean","x1.pic.trait.sd","x1.pic.trait.mediam","x1.pic.trait.skewness","x1.pic.trait.kurtosis","x1.pic.trait.K","x1.pic.trait.lambda")
  x2.sumstat.oubmbm<-array(NA,c(12,sims))
  rownames(x2.sumstat.oubmbm)<-c("x2.trait.mean","x2.trait.sd","x2.trait.median","x2.trait.skewness","x2.trait.kurtosis","x2.pic.trait.mean","x2.pic.trait.sd","x2.pic.trait.mediam","x2.pic.trait.skewness","x2.pic.trait.kurtosis","x2.pic.trait.K","x2.pic.trait.lambda")

  for(simIndex in 1:sims){
    if(simIndex %%100==0){print(paste("oubmbm",simIndex,sep=""))}
    #print(paste("oubmbm simIndex= ",simIndex,sep=""))
    prior.params<-oubmbmprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
    model.params.oubmbm[,simIndex]<-prior.params$model.params#for record only
    reg.params.oubmbm[,simIndex]<-prior.params$reg.params#for record only

    sim.trait <-oubmbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
    sim.oubmbm.trait[,1,simIndex]<-sim.trait$y
    sim.oubmbm.trait[,2,simIndex]<-sim.trait$x1
    sim.oubmbm.trait[,3,simIndex]<-sim.trait$x2
    y.sumstat.oubmbm[,simIndex]<- sumstat(trait=sim.trait$y,tree=tree)
    x1.sumstat.oubmbm[,simIndex]<- sumstat(trait=sim.trait$x1,tree=tree)
    x2.sumstat.oubmbm[,simIndex]<- sumstat(trait=sim.trait$x2,tree=tree)
  }# end of loop

  sumstat.oubmbm <- cbind(t(y.sumstat.oubmbm),t(x1.sumstat.oubmbm),t(x2.sumstat.oubmbm))
  oubmbm.par.sim <- cbind(t(model.params.oubmbm),t(reg.params.oubmbm))

  return(list(sumstat.oubmbm=sumstat.oubmbm,oubmbm.par.sim=oubmbm.par.sim))
}



###############
### ououbm ###
##############
#' Parameter samples and summary statistics
#'
#' Draw sample for parameters, simulate trait and compute the summary statistics for OUOUBM model
#'
#' Given tree, trait sets, function \code{\link{HyperParam}} is called to yield the range of parameters, then function \code{\link{oubmbmprior}} is called to draw sample for parameter, then the function \code{\link{oubmbmmodel}} is applied to simulate traits through post order tree traversal algorithm, finally the summary statistics is computed by function \code{\link{sumstat}}.
#'
#' @param tree   An ape: tree object stored in phylo format
#' @param traitset a dataframe that contains 3 traits
#' @param sims number of trait replicate
#'
#' @return A list of vectors containing a dataframe of summary statistics, and a dataframe of parameter samples
#'
#' @export
#'
#'
#' @examples
#'
#' ## using lizard dataset (running time more > 5 sec)
#' \donttest{
#' data(lizard)
#' tree<-lizard$tree
#' traitset<-lizard$traitset
#' sims<-10
#' ououbmTrait(tree=tree,traitset=traitset,sims=sims)
#'}

ououbmTrait<-function(tree=tree,traitset=traitset,sims=sims){
  prior.params<-HyperParam(tree=tree,traitset=traitset)
  prior.model.params<-prior.params$prior.model.params
  prior.reg.params<-prior.params$prior.reg.params
  prior.params.ououbm <- ououbmprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  sim.ououbm.trait<-array(NA,c(dim(traitset)[1],5,sims))
  root<-prior.params$root
  model.params.ououbm<-array(NA,c(5,sims))
  rownames(model.params.ououbm)<-c("alpha.y","alpha.x","theta.x","sigmasq.x","tau")
  reg.params.ououbm<-array(NA,c(3,sims))
  row.names(reg.params.ououbm)<-c("b0", "b1", "b2")
  y.sumstat.ououbm<-array(NA,c(12,sims))
  rownames(y.sumstat.ououbm)<-c("y.trait.mean","y.trait.sd","y.trait.median","y.trait.skewness","y.trait.kurtosis","y.pic.trait.mean","y.pic.trait.sd","y.pic.trait.mediam","y.pic.trait.skewness","y.pic.trait.kurtosis","y.pic.trait.K","y.pic.trait.lambda")
  x1.sumstat.ououbm<-array(NA,c(12,sims))
  rownames(x1.sumstat.ououbm)<-c("x1.trait.mean","x1.trait.sd","x1.trait.median","x1.trait.skewness","x1.trait.kurtosis","x1.pic.trait.mean","x1.pic.trait.sd","x1.pic.trait.mediam","x1.pic.trait.skewness","x1.pic.trait.kurtosis","x1.pic.trait.K","x1.pic.trait.lambda")
  x2.sumstat.ououbm<-array(NA,c(12,sims))
  rownames(x2.sumstat.ououbm)<-c("x2.trait.mean","x2.trait.sd","x2.trait.median","x2.trait.skewness","x2.trait.kurtosis","x2.pic.trait.mean","x2.pic.trait.sd","x2.pic.trait.mediam","x2.pic.trait.skewness","x2.pic.trait.kurtosis","x2.pic.trait.K","x2.pic.trait.lambda")
  # rownames(post.model.params.ououbm)<-c("alpha.y","alpha.x","theta.x","sigma.x","tau")
  # row.names(post.reg.params.ououbm)<-c("b0", "b1", "b2")
  for(simIndex in 1:sims){
    if(simIndex %%100==0){print(paste("ououbm",simIndex,sep=""))}
    #simIndex<-2
    #print(paste("ououbm simIndex= ",simIndex,sep=""))
    prior.params <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
    model.params.ououbm[,simIndex]<-prior.params$model.params#for record only
    reg.params.ououbm[,simIndex]<-prior.params$reg.params#for record only

    sim.trait <-ououbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
    sim.ououbm.trait[,1,simIndex]<-sim.trait$y
    sim.ououbm.trait[,2,simIndex]<-sim.trait$x1
    sim.ououbm.trait[,3,simIndex]<-sim.trait$x2
    y.sumstat.ououbm[,simIndex]<- sumstat(trait=sim.trait$y,tree=tree)
    x1.sumstat.ououbm[,simIndex]<- sumstat(trait=sim.trait$x1,tree=tree)
    x2.sumstat.ououbm[,simIndex]<- sumstat(trait=sim.trait$x2,tree=tree)
  }#end of loop
  sumstat.ououbm <- cbind(t(y.sumstat.ououbm),t(x1.sumstat.ououbm),t(x2.sumstat.ououbm))
  ououbm.par.sim <- cbind(t(model.params.ououbm),t(reg.params.ououbm))
  return(list(sumstat.ououbm=sumstat.ououbm,ououbm.par.sim=ououbm.par.sim))
}

###############
### oubmcir ###
##############
#' Parameter samples and summary statistics
#'
#' Draw sample for parameters, simulate trait and compute the summary statistics for OUBMCIR model
#'
#' Given tree, trait sets, function \code{\link{HyperParam}} is called to yield the range of parameters, then function \code{\link{oubmbmprior}} is called to draw sample for parameter, then the function \code{\link{oubmbmmodel}} is applied to simulate traits through post order tree traversal algorithm, finally the summary statistics is computed by function \code{\link{sumstat}}.
#'
#' @param tree   An ape: tree object stored in phylo format
#' @param traitset a dataframe that contains 3 traits
#' @param sims number of trait replicate
#'
#' @return A list of vectors containing a dataframe of summary statistics, and a dataframe of parameter samples
#'
#' @export
#'
#'
#' @examples
#'
#' ## using coral dataset (running time more > 5 sec)
#' \donttest{
#' data(coral)
#' tree<-coral$tree
#' traitset<-coral$traitset
#' sims<-10
#' oubmcirTrait(tree=tree,traitset=traitset,sims=sims)
#'}
oubmcirTrait<-function(tree=tree,traitset=traitset,sims=sims){
  prior.params<-HyperParam(tree=tree,traitset=traitset)
  prior.model.params<-prior.params$prior.model.params
  prior.reg.params<-prior.params$prior.reg.params

  prior.params.oubmcir <- oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  root<-prior.params$root
  sim.oubmcir.trait<-array(NA,c(dim(traitset)[1],3,sims))
  model.params.oubmcir<-array(NA,c(5,sims))
  rownames(model.params.oubmcir)<-c("alpha.y","sigmasq.x","alpha.tau","theta.tau","sigmasq.tau")
  reg.params.oubmcir<-array(NA,c(3,sims))
  row.names(reg.params.oubmcir)<-c("b0", "b1", "b2")
  y.sumstat.oubmcir<-array(NA,c(12,sims))
  rownames(y.sumstat.oubmcir)<-c("y.trait.mean","y.trait.sd","y.trait.median","y.trait.skewness","y.trait.kurtosis","y.pic.trait.mean","y.pic.trait.sd","y.pic.trait.mediam","y.pic.trait.skewness","y.pic.trait.kurtosis","y.pic.trait.K","y.pic.trait.lambda")
  x1.sumstat.oubmcir<-array(NA,c(12,sims))
  rownames(x1.sumstat.oubmcir)<-c("x1.trait.mean","x1.trait.sd","x1.trait.median","x1.trait.skewness","x1.trait.kurtosis","x1.pic.trait.mean","x1.pic.trait.sd","x1.pic.trait.mediam","x1.pic.trait.skewness","x1.pic.trait.kurtosis","x1.pic.trait.K","x1.pic.trait.lambda")
  x2.sumstat.oubmcir<-array(NA,c(12,sims))
  rownames(x2.sumstat.oubmcir)<-c("x2.trait.mean","x2.trait.sd","x2.trait.median","x2.trait.skewness","x2.trait.kurtosis","x2.pic.trait.mean","x2.pic.trait.sd","x2.pic.trait.mediam","x2.pic.trait.skewness","x2.pic.trait.kurtosis","x2.pic.trait.K","x2.pic.trait.lambda")

  for(simIndex in 1:sims){
    if(simIndex %%100==0){print(paste("oubmcir",simIndex,sep=""))}

    prior.params<-oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
    model.params.oubmcir[,simIndex]<-prior.params$model.params#for record only
    reg.params.oubmcir[,simIndex]<-prior.params$reg.params#for record only

    sim.trait <-oubmcirmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
    sim.oubmcir.trait[,1,simIndex]<-sim.trait$y
    sim.oubmcir.trait[,2,simIndex]<-sim.trait$x1
    sim.oubmcir.trait[,3,simIndex]<-sim.trait$x2
    y.sumstat.oubmcir[,simIndex]<- sumstat(trait=sim.trait$y,tree=tree)
    x1.sumstat.oubmcir[,simIndex]<- sumstat(trait=sim.trait$x1,tree=tree)
    x2.sumstat.oubmcir[,simIndex]<- sumstat(trait=sim.trait$x2,tree=tree)
  }# end of loop

  sumstat.oubmcir <- cbind(t(y.sumstat.oubmcir),t(x1.sumstat.oubmcir),t(x2.sumstat.oubmcir))
  oubmcir.par.sim <- cbind(t(model.params.oubmcir),t(reg.params.oubmcir))
  return(list(sumstat.oubmcir=sumstat.oubmcir,oubmcir.par.sim=oubmcir.par.sim))
}
###############
### ououcir ###
##############
#' Parameter samples and summary statistics
#'
#' Draw sample for parameters, simulate trait and compute the summary statistics for OUOUCIR model
#'
#' Given tree, trait sets, function \code{\link{HyperParam}} is called to yield the range of parameters, then function \code{\link{oubmbmprior}} is called to draw sample for parameter, then the function \code{\link{oubmbmmodel}} is applied to simulate traits through post order tree traversal algorithm, finally the summary statistics is computed by function \code{\link{sumstat}}.
#'
#' @param tree   An ape: tree object stored in phylo format
#' @param traitset a dataframe that contains 3 traits
#' @param sims number of trait replicate
#'
#' @return A list of vectors containing a dataframe of summary statistics, and a dataframe of parameter samples
#'
#' @export
#'
#'
#' @examples
#'
#' ## using bat dataset (running time more > 5 sec)
#' \donttest{
#' data(bat)
#' tree<-bat$tree
#' traitset<-bat$traitset
#' sims<-10
#' ououcirTrait(tree=tree,traitset=traitset,sims=sims)
#'}

ououcirTrait<-function(tree=tree,traitset=traitset,sims=sims){
  prior.params<-HyperParam(tree=tree,traitset=traitset)
  prior.model.params<-prior.params$prior.model.params
  prior.reg.params<-prior.params$prior.reg.params

  prior.params.ououcir <- ououcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  root<-prior.params$root
  sim.ououcir.trait<-array(NA,c(dim(traitset)[1],3,sims))
  model.params.ououcir<-array(NA,c(7,sims))
  rownames(model.params.ououcir)<-c("alpha.y", "alpha.x", "theta.x", "sigmasq.x","alpha.tau","theta.tau","sigmasq.tau")
  reg.params.ououcir<-array(NA,c(3,sims))
  row.names(reg.params.ououcir)<-c("b0", "b1", "b2")
  y.sumstat.ououcir<-array(NA,c(12,sims))
  rownames(y.sumstat.ououcir)<-c("y.trait.mean","y.trait.sd","y.trait.median","y.trait.skewness","y.trait.kurtosis","y.pic.trait.mean","y.pic.trait.sd","y.pic.trait.mediam","y.pic.trait.skewness","y.pic.trait.kurtosis","y.pic.trait.K","y.pic.trait.lambda")
  x1.sumstat.ououcir<-array(NA,c(12,sims))
  rownames(x1.sumstat.ououcir)<-c("x1.trait.mean","x1.trait.sd","x1.trait.median","x1.trait.skewness","x1.trait.kurtosis","x1.pic.trait.mean","x1.pic.trait.sd","x1.pic.trait.mediam","x1.pic.trait.skewness","x1.pic.trait.kurtosis","x1.pic.trait.K","x1.pic.trait.lambda")
  x2.sumstat.ououcir<-array(NA,c(12,sims))
  rownames(x2.sumstat.ououcir)<-c("x2.trait.mean","x2.trait.sd","x2.trait.median","x2.trait.skewness","x2.trait.kurtosis","x2.pic.trait.mean","x2.pic.trait.sd","x2.pic.trait.mediam","x2.pic.trait.skewness","x2.pic.trait.kurtosis","x2.pic.trait.K","x2.pic.trait.lambda")

  for(simIndex in 1:sims){
    if(simIndex %%100==0){print(paste("ououcir",simIndex,sep=""))}
    prior.params<-ououcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
    model.params.ououcir[,simIndex]<-prior.params$model.params#for record only
    reg.params.ououcir[,simIndex]<-prior.params$reg.params#for record only

    sim.trait <-ououcirmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
    sim.ououcir.trait[,1,simIndex]<-sim.trait$y
    sim.ououcir.trait[,2,simIndex]<-sim.trait$x1
    sim.ououcir.trait[,3,simIndex]<-sim.trait$x2
    y.sumstat.ououcir[,simIndex]<- sumstat(trait=sim.trait$y,tree=tree)
    x1.sumstat.ououcir[,simIndex]<- sumstat(trait=sim.trait$x1,tree=tree)
    x2.sumstat.ououcir[,simIndex]<- sumstat(trait=sim.trait$x2,tree=tree)
  }# end of loop

  sumstat.ououcir <- cbind(t(y.sumstat.ououcir),t(x1.sumstat.ououcir),t(x2.sumstat.ououcir))
  ououcir.par.sim <- cbind(t(model.params.ououcir),t(reg.params.ououcir))
  return(list(sumstat.ououcir=sumstat.ououcir,ououcir.par.sim=ououcir.par.sim))
}

###############
### Analysis ##
###############

#' main program to perform analysis
#'
#' Analyze data and report the model estimates and model selection
#'
#'\code{\link{ouxy}} performs data analaysis under Approximate Bayesian Computation(ABC) procedure.  The summary statistics for the raw traitsets are first computed by by function \code{\link{sumstat}}, and the parameters ranges are computed using the tree and tratisets under function \code{\link{HyperParam}}, and sample of prior paramters are drawn from function \code{\link{oubmbmprior}}, then the function \code{\link{oubmbmmodel}} is applied to simulate traits through post order tree traversal algorithm.
#' The ABC procedure are then performed using sample of paramters and simulated traitset.
#' Posterior sample are chosen using acceptance rate \code{sims * tol}. The posterior samples are computed using rejection method \code{\link[abc]{abc}} to median of the posterior samples are as reported parameter esitmate  and Bayes factor is computed using function \code{\link[abc]{postpr}} accordingly by the ratio of the posterior model probability under each model.
#'
#' @param tree   An ape: tree object stored in phylo format
#' @param traitset a dataframe that contains 3 traits
#' @param tol acceptance rate from ABC
#' @param sims number of trait replicate
#'
#'
#' @return A list of vectors containing a dataframe of model parameter estimate, and a dataframe of Bayes factors between a pair of models
#' \enumerate{
#'   \item \strong{table.output}: The posterior median for parameter estiamtes under each model.
#'   \item \strong{s.mnlog}: Bayes factor tables comparing a pair of models.
#' }
#'
#'
#' @export
#'
#'
#' @examples
#'
#' ## using coral dataset (It takes for a whiles)
#' \donttest{
#' data(coral)
#' tree<-coral$tree
#' traitset<-coral$traitset
#' sims<-1000
#' output<-ouxy(tree=tree,traitset=traitset,tol=0.1,sims= sims)
#'
#' ## OUTPUT THE FOLLOWING
#' ## >output$s.mnlog
#'
#' ## $mnlogistic
#' ## $mnlogistic$Prob
#' ## oubmbm    oubmcir     ououbm    ououcir
#' ## 0.03081341 0.01533086 0.40779579 0.54605995
#'
#' ## $mnlogistic$BayesF
#' ##           oubmbm     oubmcir      ououbm     ououcir
#' ## oubmbm   1.00000000  2.00989403  0.07556087  0.05642861
#' ## oubmcir  0.49753867  1.00000000  0.03759446  0.02807542
#' ## ououbm  13.23436292 26.59966708  1.00000000  0.74679673
#' ## ououcir 17.72150620 35.61834960  1.33905246  1.00000000
#' ##
#' ## > output$table.out
#' ##          alpha.y alpha.x alpha.tau theta.x theta.tau   sigma.x
#' ## OUBMBM   4.3064      NA        NA      NA        NA  7.821074
#' ## OUOUBM   4.1240  5.2119        NA -0.5759        NA 10.117253
#' ## OUBMCIR  4.3720      NA    4.0736      NA    1.2326  7.825912
#' ## OUOUCIR  3.1016  4.4269    3.9930  0.0668    1.2702  9.226803
#' ## GLS          NA      NA        NA      NA        NA        NA
#' ##
#' ##           tau sigma.tau        b0         b1        b2
#' ## OUBMBM  2.2403        NA 0.1678000 0.03850000 0.2874000
#' ## OUOUBM  2.5021        NA 0.1651000 0.03260000 0.3146000
#' ## OUBMCIR     NA  1.492548 0.1706000 0.03760000 0.3049000
#' ## OUOUCIR     NA  1.516047 0.1661000 0.03480000 0.2549000
#' ## GLS         NA        NA 0.1682413 0.03931911 0.3564761
#'}
#'
#'
#'
#'
ouxy<-function(tree=tree,traitset=traitset,tol=0.1,sims= 100){
  resptrait<-traitset$resptrait
  predtrait1<-traitset$predtrait1
  predtrait2<-traitset$predtrait2
  names(resptrait)<-tree$tip.label
  names(predtrait1)<-tree$tip.label
  names(predtrait2)<-tree$tip.label

  raw.sumstat.y <- sumstat(trait = resptrait, tree=tree)
  raw.sumstat.x1 <- sumstat(trait = predtrait1, tree=tree)
  raw.sumstat.x2 <- sumstat(trait = predtrait2, tree=tree)
  raw.sumstat<-c(raw.sumstat.y,raw.sumstat.x1,raw.sumstat.x2)

  models<-rep(c("oubmbm","ououbm","oubmcir","ououcir"), each=sims)

  sampleoubmbm<-oubmbmTrait(tree=tree,traitset=traitset,sims=sims)
  sampleououbm<-ououbmTrait(tree=tree,traitset=traitset,sims=sims)
  sampleoubmcir<-oubmcirTrait(tree=tree,traitset=traitset,sims=sims)
  sampleououcir<-ououcirTrait(tree=tree,traitset=traitset,sims=sims)


  oubmbm.par.sim<-sampleoubmbm$oubmbm.par.sim
  ououbm.par.sim<-sampleououbm$ououbm.par.sim
  oubmcir.par.sim<-sampleoubmcir$oubmcir.par.sim
  ououcir.par.sim<-sampleououcir$ououcir.par.sim

  sumstat.oubmbm <- sampleoubmbm$sumstat.oubmbm
  sumstat.ououbm <- sampleououbm$sumstat.ououbm
  sumstat.oubmcir <- sampleoubmcir$sumstat.oubmcir
  sumstat.ououcir <- sampleououcir$sumstat.ououcir

  #######################
  ### posterior sample###
  ###   mnlogistic    ###
  #######################

  rejection.result.oubmbm <- abc(target = raw.sumstat,param = data.frame(oubmbm.par.sim), sumstat = sumstat.oubmbm,tol = tol,method="rejection")
  rejection.result.ououbm <- abc(target = raw.sumstat,param = data.frame(ououbm.par.sim), sumstat = sumstat.ououbm,tol = tol,method="rejection")
  rejection.result.oubmcir <- abc(target = raw.sumstat,param = data.frame(oubmcir.par.sim), sumstat = sumstat.oubmcir,tol = tol,method="rejection")
  rejection.result.ououcir <- abc(target = raw.sumstat,param = data.frame(ououcir.par.sim), sumstat = sumstat.ououcir,tol = tol,method="rejection")

  post.oubmbm <- as.data.frame(rejection.result.oubmbm$unadj.values)
  post.ououbm <- as.data.frame(rejection.result.ououbm$unadj.values)
  post.oubmcir <- as.data.frame(rejection.result.oubmcir$unadj.values)
  post.ououcir <- as.data.frame(rejection.result.ououcir$unadj.values)

  #######################
  ### model selection ###
  #######################

  full.sumstat <- rbind(sumstat.oubmbm, sumstat.ououbm, sumstat.oubmcir, sumstat.ououcir)
  modsel.mnlog<-postpr(target=raw.sumstat,models,full.sumstat, tol=tol,method="mnlogistic")
  s.mnlog<-summary(modsel.mnlog)

  ####################################
  ### Posterior Parameter Estimate ###
  ####################################
  post.oubmbm.mean <- round(apply(post.oubmbm,2,mean),digits = 4)
  post.ououbm.mean <- round(apply(post.ououbm,2,mean),digits = 4)
  post.oubmcir.mean <- round(apply(post.oubmcir,2,mean),digits = 4)
  post.ououcir.mean <- round(apply(post.ououcir,2,mean),digits = 4)

  table.params <- data.frame(matrix(NA,4,8))
  rownames(table.params) <- c("OUBMBM","OUOUBM","OUBMCIR","OUOUCIR")
  colnames(table.params)<- c("alpha.y","sigmasq.x","tau","alpha.x","theta.x","alpha.tau","theta.tau","sigmasq.tau")

  table.params["OUBMBM",c("alpha.y","sigmasq.x","tau")]<-post.oubmbm.mean[c("alpha.y","sigmasq.x","tau")]


  table.params["OUOUBM",c("alpha.y","alpha.x","theta.x", "sigmasq.x", "tau")]<-post.ououbm.mean[c("alpha.y","alpha.x","theta.x", "sigmasq.x", "tau")]

  table.params["OUBMCIR",c("alpha.y", "sigmasq.x","alpha.tau","theta.tau", "sigmasq.tau")]<-post.oubmcir.mean[c("alpha.y", "sigmasq.x","alpha.tau","theta.tau", "sigmasq.tau")]

  table.params["OUOUCIR",c("alpha.y","alpha.x","theta.x","sigmasq.x","alpha.tau","theta.tau", "sigmasq.tau")]<-post.ououcir.mean[c("alpha.y","alpha.x","theta.x","sigmasq.x","alpha.tau","theta.tau", "sigmasq.tau")]
  table.params[,"sigmasq.x"]<-sqrt(table.params[,"sigmasq.x"])
  table.params[,"sigmasq.tau"]<-sqrt(table.params[,"sigmasq.tau"])

  colnames(table.params)<- c("alpha.y","sigma.x","tau","alpha.x","theta.x","alpha.tau","theta.tau","sigma.tau")
  table.params<-table.params[, c("alpha.y","alpha.x", "alpha.tau", "theta.x","theta.tau", "sigma.x","tau","sigma.tau")]

  table.b<- matrix(NA,4,3)
  rownames(table.b) <- c("OUBMBM","OUOUBM","OUBMCIR","OUOUCIR")
  colnames(table.b)<-c("b0","b1","b2")
  table.b[1,]<-post.oubmbm.mean[c(4:6)]
  table.b[2,]<-post.ououbm.mean[c(6:8)]
  table.b[3,]<-post.oubmcir.mean[c(6:8)]
  table.b[4,]<-post.ououcir.mean[c(8:10)]

  gls<-lm(resptrait~predtrait1+predtrait2)
  gls$coefficients
  GLS<-rep(NA,11)
  GLS[9:11]<-gls$coefficients
  table.output<-rbind(cbind(table.params,table.b),GLS)
  rownames(table.output) <- c("OUBMBM","OUOUBM","OUBMCIR","OUOUCIR","GLS")
  return(list(s.mnlog=s.mnlog,table.out=table.output))
}
