#' Check inputs to blin_fit()
#'
#' @keywords internal
check_args <- function(Y,lag,X=NULL)
{
  if(length(dim(Y))!=3){stop("Y must be a 3-mode array")}
  if(dim(Y)[3] <= lag){stop("Need more than `lag` observations of Y")}
  if(sum(is.na(Y)) > 0){warning("There are NAs in Y; omitting these observations from model fit or imputing as reasonable")}
  if(!is.null(X)){
    if(length(dim(X))!=4){stop("X must be a 4-mode array")}
    if(any(dim(X)[1:3] != dim(Y))){stop("First three dimensions of X must be equal those of Y")}
    if(sum(is.na(X)) > 0){
      warning("There are NAs in X; setting these values to zero")
      ix <- which(is.na(X))
      iy <- which(is.na(Y))
      if(length(ix) == length(iy)){
        if(any(sort(ix) != sort(iy))){
          warning("NAs in X and Y do not match")
        }
      } else {
        warning("NAs in X and Y do not match")
      }
      X[is.na(X)] <- 0
    }
  }
  
  return(list(X=X))
}


#' Check inputs for reduced rank fit of BLIN model
#'
#' @keywords internal
check_ranks <- function(rankA, rankB, S, L)
{
  if(is.null(rankA) | is.null(rankB)){stop("rankA and rankB must be specified for reduced rank fit")}
  if(rankA < 1 | rankA > S){stop("rankA must be between 1 and dim(Y)[1]")}
  if(rankB < 1 | rankB > L){stop("rankB must be between 1 and dim(Y)[2]")}
}


#' Lag Y array to get autoregressive array D
#'
#' @keywords internal
lag_Y <- function(lag,Y,X=NULL)
{
  tmax <- dim(Y)[3]
  Ynew <- Y[,,(lag+1):tmax]
  Dnew <- 0*Ynew
  if(lag > 1){
    for(t in 1:(tmax - lag)){
      Dnew[,,t] <- apply(Y[,,t + 1:lag - 1], 1:2, sum)   # sum previous lags
    }
  } else {
    Dnew <- Y[,,1:(tmax-1)]
  }
  
  if(!is.null(X)){
    X <- X[,,(lag+1):tmax,]
  }
  
  return(list(D=Dnew, Y=Ynew, X=X))
}



#' Find diagonal indices of k-mode array, where "diagonal" is i=j for first two indices of input array Y
#'
#' @keywords internal
adiag <- function(Y)
{
  if(length(dim(Y)) < 3){ stop("Y is not a 3-mode array")}
  if(dim(Y)[1] != dim(Y)[2]){ stop("Y is not square in first two dimensions")}
  
  n <- dim(Y)[1]
  rest <- as.matrix(expand.grid(lapply(dim(Y)[-c(1,2)], function(x) 1:x)))   # combinations of all dimensions beyond first two
  r <- length(dim(Y)) - 2   # number of additional dimensions
  first2 <- cbind(rep(1:n, times=nrow(rest)), rep(1:n, times=nrow(rest)))
  rest <- matrix(rep(rest, each=n), ncol=r)
  
  return(cbind(first2, rest))
}


#' Check whether there is a unique solution to the BLIN model
#'
#' @keywords internal
check_unique <- function(S,L,tmax)
{
  temp <- tmax*S*L
  if(temp < (S^2 + L^2 - 1)){warning("BLIN design matrix is rank-deficient; there is no unique solution to the least squares criterion")}
}


#' Compute generalized inverse of X^T X, for X 
#'
#' @keywords internal
ginvXX <- function(D,X=NULL)
{
  
  tmax <- dim(D)[3]
  DDT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x])))    # t(D) %*% D
  DTD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x])))    # D %*% t(D)
  
  
  
}



#' Multiply an array by a matrix along a given mode, by
#' @author Peter Hoff
#' 
#' This function multiplies a matricized array by another 
#' matrix, and then reforms the result into a new array. 
#'
#' @keywords internal
amprod<-function(A,M,k)
{
  K<-length(dim(A))
  AM<-M%*%mat(A,k)
  AMA<-array(AM, dim=c(dim(M)[1],dim(A)[-k]) )
  aperm(AMA,  match(1:K,c(k,(1:K)[-k]) ) )
}



#' Matricize an array
#'
#' This functions matricizes an array along a given mode. 
#'
#' @author Peter Hoff
#' 
#' @keywords internal
mat<-function(A,k)
{
  Ak<-t(apply(A,k,"c"))
  if(nrow(Ak)!=dim(A)[k])  { Ak<-t(Ak) }
  Ak
}


