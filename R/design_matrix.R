#' Build the BLIN design matrix
#'
#' @export
#' @keywords external
#' 
#' @param Y Response 3-mode array.
#' @param X Optional 4-mode array of covariates, defaults to no covariates.
#' @param lag Optional numeric specifying autoregressive lag in model, defaults to 1. 
#' @param showWarnings Optional logical whether matrix memory size should be evaluated and warning provided (see details), defaults to TRUE. 
#' 
#' @details This function takes an \eqn{S \times L \times T} array \eqn{Y} that is a representation of a longitudinal bipartite relational data set. 
#' Optional input is an \eqn{S \times L \times T \times p} array  \eqn{X} of covariates that influence the evolution of the data set in equation over time.
#' The function returns an \eqn{(SL(T - lag)) \times (S^2 + L^2 + p)} design matrix, of sparse class, upon which \code{Y[,,lag:T]} may be regressed.
#' If \code{showWarnings = TRUE}, and if the estimated size of the design matrix is greater than 1GB, a warning is thrown. 
#' 
#'
#' @return 
#' \item{}{A sparse design matrix}
#'
#' @seealso \code{\link{generate_blin}} \code{\link{blin_mle}}
#' 
#'
#' @examples
#' S <- 5
#' L <- 4
#' tmax <- 10
#' data <- generate_blin(S,L,tmax, lag=2, sparse=.8, seed=1)
#' dim(data$Y)
#' 
#' Xreg <- build_design(data$Y, data$X, lag=2)
#' dim(Xreg)
#' class(Xreg)
#' 
build_design <- function(Y, X=NULL, lag=1, showWarnings=TRUE)
{
  lag <- as.numeric(lag)
  showWarnings <- as.logical(showWarnings)
  
  if(length(dim(Y))!=3){"Y must be a 3-mode array"}
  type="biten"
  
  use_cov=!is.null(X)

  S <- L <- tmax <- NULL
  
  data <- lag_Y(lag, Y, X)
  D <- data$D
  X <- data$X
  
  
  # Size check
  if(showWarnings){
    nd <- sum(data$D != 0)
    nx <- sum(data$X != 0)
    size <- 8*nd*nx/1e9
    if(size >= 1){
      warning("Building desing matri that may be larger than 1 GB")
    }
  }
  
  # if(!sparsedata){

    # Find sizes
    S = nrow(D[,,1])
    L = ncol(D[,,1])
    tmax = dim(D)[3]
    
    # BLIN matrices
    Js <- diag(S)
    Jl <- diag(L)
    
    # for(t in 1:tmax){  # vec A then vec B
    #   Xreg[ 1:(S*L) + (t - 1)*S*L, 1:S^2] <- kronecker(t(D[,,t] %*% Jl), diag(S))   # A columns
    #   Xreg[ 1:(S*L) + (t - 1)*S*L, S^2 + 1:L^2] <- kronecker(diag(L), Js %*% D[,,t])    # B columns
    # }
    # 
    ikeep <- which(D!=0, arr.ind=TRUE)
    Darray <- cbind(ikeep, D[ikeep])
    
    if(use_cov){
      if(length(dim(X))!=4){"X must be a 4-mode array"}
      p <- dim(X)[4]
      ikeepx <-  which(X!=0, arr.ind=TRUE)
      ikeepx <- unique(ikeepx[, 1:3])
      Xarray <- ikeepx
      for(k in 1:p){
        Xarray <- cbind(Xarray, X[cbind(ikeepx, k)]) 
      }
    }
  
    
    response <- "response"
    D <- data.frame(Darray)
    names(D) <- c("i", "j", "t", response)
    X <- data.frame(Xarray)
    names(X)[1:3] <- c("i", "j", "t")

      numcols <- S^2 + L^2 
      if(use_cov){
        p <- ncol(X) - 3   # remove i,j,t values
        numcols <- numcols + p
        
        X <- X[X$i <= S & X$j <= L & X$t <= tmax,]
      }
      
      # Xreg <- sparseMatrix(i=1,j=1,x=0, dims=c(tmax*S*L, numcols))   # initialize
      onerows <- which(as.vector(D[,response]) != 0)    # rows of D that have 1s in response
      # reg1s <- matrix(0, length(onerows)*(S+L), 2)   # indices in Xreg that are 1s
      count <- 0
      Xreg <- sparseMatrix(i=1, j=1, x=0, dims=c(tmax*S*L, numcols))   # initialize
      
      for(k in onerows){
        count <- count+1
        
        i <- D$i[k]    ;   j <- D$j[k]   ;   t <- D$t[k]

        ia <- cbind(S*L*(t-1) + S*(j -1) + 1:S, S*(i-1) + 1:S)
        ib <- cbind(S*L*(t-1) + i + S*(0:(L-1)), S^2 + j + L*(0:(L-1)))
        Xreg[ia] <- Xreg[ib] <- D[k, response]
      }

      if(use_cov){
        keep <- which(!(names(X) %in% c("i","j","t")))  # columns to keep that aren't i,j, or t
        Xreg[X$i + (X$j-1)*S + (X$t-1)*S*L, S^2 + L^2 + 1:p] <- as.matrix(X[, keep])
      }
      
  
  return(Xreg)
}

