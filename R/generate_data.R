#' Generate data from the continuous BLIN model
#'
#' This function generates data from the bipartite logitudinal influence network (BLIN) model \eqn{Y_t = A^T \sum_{k=1}^{lag} Y_{t-k} + \sum_{k=1}^{lag} Y_{t-k} B + X_t \beta + \tau E_t}.
#'
#' @export
#' @keywords external
#'
#' @param S Dimension of A.
#' @param L Dimension of B.
#' @param tmax Number of observations of relational data.
#' @param lag Autoregressive lag in model, defaults to 1. 
#' @param tau Optional error standard deviatiom, defaults to 1.
#' @param sigmaY Optional standard deviation of entries in \eqn{Y_t}, defaults to 1. 
#' @param muAB Optional mean of entries in decomposition of matrices \eqn{A = UV^T} and \eqn{B = WZ^T}, defaults to 0. 
#' @param sigmaAB Optional standard deviation of entries in decomposition matrices of \eqn{A = UV^T} and \eqn{B = WZ^T}, defaults to 1. 
#' @param rankA Rank of influence network matrix \eqn{A}, defaults to full rank.
#' @param rankB Optional rank of influence network matrix \eqn{B}, defaults to full rank.
#' @param use_cov Optional logical used to indicate whether to include \eqn{X_t \beta} in the model (\code{TRUE}) or not (\code{FALSE}), defaults to \code{TRUE}.
#' @param seed Optional numeric to set seed before generating, defaults to NA (no seed set).
#' @param sparse Optional degree of sparsity in A and B, i.e. \code{sparsity=.9} means 10\% of the entries in A and B are set to zero at random. Defaults to \code{NA} (no entries set to zero).
#' 
#' @details This function generates a continuous bipartite longitudinal relational data set from the BLIN model,
#'  \eqn{Y_t = A^T \sum_{k=1}^{lag} Y_{t-k} + \sum_{k=1}^{lag} Y_{t-k} B + X_t \beta + \tau E_t}, where \eqn{ \{ Y_t \}_t } is a set of \eqn{S \times L} matrices representing the bipartite relational data at each observation \eqn{t}. 
#'  The set \eqn{\{X_t \}_t} is a set of \eqn{S \times L \times p} arrays describing the influence of the
#'  coefficient vector \eqn{beta}. Finally, each matrix \eqn{E_t} consists of iid standard normal random variables.
#'  
#'  The matrices \eqn{A} and \eqn{B} are square matrices respesenting the influence networks among \eqn{S} senders and \eqn{L} receivers, respectively. The matrix \eqn{A} has decomposition \eqn{A = UV^T}, where each of \eqn{U} and \eqn{V} is an \eqn{S \times {rankA}} matrix of iid standard normal random variables with mean \code{muAB} and standard deviation \code{sigmaAB}. 
#'  Similarly, the matrix \eqn{B} has decomposition \eqn{B = WZ^T}, where each of \eqn{W} and \eqn{Z} is an \eqn{L \times {rankB}} matrix of iid standard normal random variables with standard deviation \code{sigmaAB} and mean \code{muAB} for \eqn{W} and mean \code{-muAB} for \eqn{Z}. 
#'  Lastly, the covariate array \eqn{X_t} has 3 covariates: the first is an intercept, the second consists of iid Bernoulli random variables, and the third consists of iid standard normal random variables. All coefficients are \eqn{\beta_i = 0} for \eqn{i = 1,2,3}.
#'
#' @return 
#' \item{fit}{An \code{blin} object containing summary information.}
#'
#' @seealso \code{\link{blin_mle}}
#' 
#'
#' @examples
#' S <- 5
#' L <- 4
#' tmax <- 10
#' data <- generate_blin(S,L,tmax, lag=2, sparse=.8)
#' names(data)
#' dim(data$X)
#' data$A
#' 
generate_blin <- function(S, L, tmax, lag=1, tau=1, sigmaY=1, 
                          muAB=0, sigmaAB=1, rankA=S, rankB=L, 
                          use_cov=TRUE, seed=NA, sparse=NA)
{
  
  binary <- FALSE
  gen_type="biten"

  # Set seed if desired
  if(is.numeric(seed)){ set.seed(seed) }

  
  # Generate initial Y
  Y <- array(rnorm(S*L*tmax, 0, sigmaY), c(S, L, tmax))
  if(tmax <= lag){
    stop("Input 'tmax' must be larger than lag.")
  }

  # Generate X
  if(use_cov){
    X1 <- array(1, c(S,L,tmax,1))
    X2 <- array(sample(c(0,1), S*L*tmax, replace=T), c(S,L,tmax,1))
    X3 <- array(rnorm(S*L*tmax), c(S,L,tmax,1))
    X <- abind::abind(X1,X2,X3)
    p <- dim(X)[4]
    beta_true <- matrix(rep(1,p), nrow=1)
    Xbeta <- drop(amprod(X, beta_true, 4))
  } else {
    X <- NULL
    Xbeta <- 0
    beta_true <- NA
  }
  
  # Generate A and B^T
  U_true <- matrix(rnorm(S*rankA, muAB, sigmaAB), S, rankA)
  V_true <- matrix(rnorm(S*rankA, muAB, sigmaAB), S, rankA)
  W_true <- matrix(rnorm(L*rankB, muAB, sigmaAB), L, rankB)
  Z_true <- matrix(rnorm(L*rankB, -muAB, sigmaAB), L, rankB)
  
  A_true <- tcrossprod(U_true, V_true)
  BT_true <- tcrossprod(Z_true, W_true)
  
  if(is.numeric(sparse)){   # set sparse % of elements to zero
    if(sparse >=0 & sparse <= 1){ 
      Aind <- matrix(sample(c(0,1), S^2, replace=T, prob=c(1-sparse, sparse)), S, S)
      Bind <- matrix(sample(c(0,1), L^2, replace=T, prob=c(1-sparse, sparse)), L, L)
      A_true <- Aind*A_true
      BT_true <- Bind*BT_true
    } else {stop("Input 'sparse' must be a numeric between zero and 1 or FALSE")}
  }
  
  # # Generate Y's
  # if(binary){
  #   # E <- array(rnorm(S*L*tmax, 0, tau*pi/sqrt(3)), c(S,L,tmax))
  #   E <- array(rlogis(S*L*tmax, 0, tau), c(S,L,tmax))
  # } else {
  #   E <- array(rnorm(S*L*tmax, 0, tau), c(S,L,tmax))
  # }
  
  # if (strtrim(gen_type,3) == "bil"){
  #   Yout <- Xbeta + amprod(amprod(D, A_true, 1), BT_true, 2) + E   # bilinear multiliplicative model
  #   
  # } else if (strtrim(gen_type,3) == "sad") {  
  #   Jl <- matrix(1, L, n)
  #   Js <- matrix(1, S, m)
  #   Yout <- Xbeta + amprod(amprod(D, A_true, 1), Jl, 2) + amprod(amprod(D, Js, 1), BT_true, 2) + E
  #   
  # } else 
  if (strtrim(gen_type,3) == "bit") { 
    A_true <- A_true*S^1.5/rankA
    BT_true <- BT_true*L^1.5/rankB   # increase variability for biten models
    
    A_true <- A_true / 2 / max(abs(A_true))
    BT_true <- BT_true / 2 / max(abs(BT_true))
    E <- tau*array(rnorm(S*L*tmax), c(S,L,tmax))
    for(t in (lag+1):tmax){
      if(lag>1){
        D <- apply(Y[,,(t-lag):(t-1),drop=FALSE], 1:2, sum) #, drop=F)
      } else {
        D <- Y[,,t-1,drop=TRUE]
      }
      if(use_cov){
        Xbt <- Xbeta[,,t, drop=TRUE]
      } else {
        Xbt <- 0
      }
      Y[,,t] <- Xbt + A_true %*%D + D %*% t(BT_true) + E[,,t,drop=TRUE]
    }
    
  } else { stop("Invalid model type for prediction")}
  
  # LLtrue <- -length(Yout)*log( sum( E^2) ) - length(Yout)   # Calc true 2*LL 
  
  # if(binary){
  #   p <- c((1+exp(-Yout + E))^{-1})   # probabilities
  #   Yout <- 1*(Yout>0)  # threshold if binary
  #   LLtrue <- 2*sum(c(Yout)*log(p) + (1-c(Yout))*log(1-p))
  # }
  
  return(list(Y=Y, X=X, E=E, beta=beta_true, A=t(A_true), B=t(BT_true), call=match.call()))
}

