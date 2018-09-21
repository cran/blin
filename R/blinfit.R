#' Estimate the BLIN model using maximum likelihood estimator
#'
#' This function estimates the bipartite logitudinal influence network (BLIN) model \eqn{Y_t = A^T \sum_{k=1}^{lag} Y_{t-k} + \sum_{k=1}^{lag} Y_{t-k} B + X_t \beta + \tau E_t} using maximum likelihood estimator.
#'
#' @export
#' @keywords external
#' 
#' @param Y Response 3-mode array.
#' @param X Optional 4-mode array of covariates, defaults to no covariates.
#' @param type Optional string specifying BLIN model type: full, reduced_rank, or sparse. Defaults to full. 
#' @param lag Optional numeric specifying autoregressive lag in model, defaults to 1. 
#' @param rankA Optional numeric rank of influence network matrix \eqn{A} for reduced rank model type, defaults to full rank.
#' @param rankB Optional numeric rank of influence network matrix \eqn{B}, defaults to rank of \eqn{A}.
#' @param maxit Optional numeric maximum number of iterations for full and reduced rank block coordinate descents, defaults to 1e3.
#' @param tol Optional numeric convergence tolerance for full and reduced rank block coordinate descents, defaults to 1e-8.
#' @param init Optional string specifying initialization type for full and reduced rank block coordinate descents, defaults to "I", identity for \eqn{A} and \eqn{B}. Also allows "random" for random initialization of \eqn{A} and \eqn{B}.
#' @param sigma_init Optional numeric standard deviation for random initialization of \eqn{A} and \eqn{B} in  full and reduced rank block coordinate descents, defaults to 1.
#' @param verbose Optional logical specifying whether progress should be printed out (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{FALSE}.
#' @param calcses  Optional logical specifying whether standard errors should be calculated (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{FALSE}. Only standard errors for the full BLIN model are implemented.
#' @param randseed Optional numeric specifying seed for random ininitialization of \eqn{A} and \eqn{B} in full and reduced rank block coordinate descents, defaults to \code{NA} (no seed set).
#' 
#' @details This function estimates the continuous BLIN model,
#'  \deqn{Y_t = A^T Y_{t-1} + Y_{t-1} B + X_t \beta + \tau E_t}, where \eqn{ \{ Y_t\}_t } is a set of \eqn{S \times L} matrices representing the bipartite relation data at each observation \eqn{t}. 
#'  The set \eqn{\{X_t \}_t} is a set of \eqn{S \times L \times p} arrays describing the influence of the
#'  coefficient vector \eqn{beta}. Finally, each matrix \eqn{E_t} is assumed to consist of iid standard normal random variables. The matrices \eqn{A} and \eqn{B} are square matrices respesenting the influence networks among \eqn{S} senders and \eqn{L} receivers, respectively. 
#'  
#'  This function estimates the BLIN model using maximum likelihood (and related) methods. The "full" model places no restrications on the influence networks \eqn{A} and \eqn{B}, and estimates
#'  these matrices (along with \eqn{\beta}) by block coordinate descent. In addition, if \code{calcses==TRUE}, the standard errors for each coefficient will be estimated. Note that the standard error procedure 
#'  may require large amounts of memory to build the BLIN design matrix; a warning is produced if the estimated size of the desgn is greater than 0.5GB. 
#'  
#'  The "reduced rank" BLIN model assumes that the matrix \eqn{A} has decomposition \eqn{A = UV^T}, where each of \eqn{U} and \eqn{V} is an \eqn{S \times \code{rankA}} matrix, and 
#'  the matrix \eqn{B} has decomposition \eqn{B = WZ^T}, where each of \eqn{W} and \eqn{Z} is an \eqn{L \times \code{rankB}} matrix. This model is also estimated using block coordinate descent. 
#'  
#'  Finally, the "sparse" BLIN model assumes that \eqn{A} and \eqn{B} matrices have many entries that are small or zero. The \code{cv.glmnet(.)} function from the \code{glmnet} package is used 
#'  to estimate the entries in \eqn{A}, \eqn{B}, and \eqn{beta}. The object resuling from \code{cv.glmnet(.)} is returned in this case. 
#'  
#'  Notice that the diagonals of \eqn{A} and \eqn{B} are not identifiable. However, the sum of each diagonal entry in \eqn{A} and \eqn{B}, i.e. \eqn{a_{ii} + b_{jj}}, is identifiable. Thus, 
#'  the diagonal sums are broken out as separate estimates under the name \code{diagAB}.
#'  
#'  If \code{calcses = TRUE} and \code{type = full}, then standard errors will be returned. These standard errors are based on the assumption that each \eqn{E_t} consists of iid standard normal random variables. 
#'  In this case, the full design matrix is built, which we call \eqn{W} here. Then, the variance-covariance matrix of the estimated coefficients is formed by \eqn{\hat{\tau}^2 (W^T W)^{-1}}, where \eqn{\hat{\tau}^2} is the usual unbiased estimator of the error variance.
#'
#' @return 
#' \item{fit}{A \code{blin} object containing summary information.}
#'
#' @seealso \code{\link{generate_blin}} \code{\link{build_design}}
#' 
#'
#' @examples
#' S <- 5
#' L <- 4
#' tmax <- 10
#' data <- generate_blin(S,L,tmax, lag=2, sparse=.8, seed=1)
#' 
#' fit <- blin_mle(data$Y, data$X, lag=2, calcses=TRUE)
#' summary(fit)
#' 
blin_mle <- function(Y, X=NULL, type="full", lag=1, rankA=NULL, rankB=rankA, 
                     maxit=1e3, tol=1e-8, init="I", sigma_init=1, verbose=FALSE, 
                     calcses=FALSE, randseed=NA)
{
  link <- "gaussian"
  
  if(strtrim(link,1) == "g"){
    # check inputs
    r <- check_args(Y,lag,X)
    X <- Xin <- r$X
    use_cov <- !is.null(X)  # flag on whether accounting for additional covariate information
    Yin <- Y
    
    # Lag and process data for fitting
    r <- lag_Y(lag,Y,X)
    D <- r$D
    Y <- r$Y
    X <- r$X
    rm(r)
    gc()
    
    # dimensions
    S <- dim(Y)[1]
    L <- dim(Y)[2]
    tmax <- dim(Y)[3]
    check_unique(S,L,tmax)
    if(use_cov){p <- ncol(X)} else {p <- 0}
    n <- sum(!is.na(Y))

    if(strtrim(type,1) == "f" | strtrim(type,1) == "F"){   # full fit
      type <- "full"
      D[is.na(D)] <- 0
      out <- fit_MLE_array_additive(Y, D, X, lag,
                                    tol=tol, maxit=maxit, use_cov=use_cov, 
                                    init=init, sigma_init=sigma_init, 
                                    verbose=verbose, Yold=Y[,,1:lag, drop=FALSE],
                                    type="biten", randseed=randseed)
      glmnetobj <- NA
      
    } else if (strtrim(type,1) == "R" | strtrim(type,1) == "r"){  # reduced rank fit
      type <- "reduced_rank"
      D[is.na(D)] <- 0
      check_ranks(rankA, rankB, S, L)
      if(rankA == S & rankB == L){   # do full rank fit
        warning("Reduced rank model specified but full ranks input; fitting full rank BLIN model.")
        type <- "full"
        out <- fit_MLE_array_additive(Y, D, X, lag,
                                      tol=tol, maxit=maxit, use_cov=use_cov, 
                                      init=init, sigma_init=sigma_init, Yold=Y[,,1:lag, drop=FALSE],
                                      verbose=verbose, type="biten", randseed=randseed)
        
      } else {
        out <- fit_MLE_array(Y, D, X, lag,
                             rankA=rankA, rankB=rankB, Yold=Y[,,1:lag, drop=FALSE],
                             tol=tol, init=init, sigma_init=sigma_init, use_cov=use_cov, 
                             maxit=maxit, randseed=randseed, verbose=verbose)
      }
      glmnetobj <- NA
      
    } else if (strtrim(type,1) == "S" | strtrim(type,1) == "s"){  # sparse fit
      type <- "sparse"
      D[is.na(D)] <- 0
      out <- additive_regression(Yin, Xin, lag, type="biten", use_cov=use_cov)   
      glmnetobj <- out$cv.glmnet.object
      
    } else { stop("invalid model type input")}
    
  } else {
    stop("only implemented for Gaussian link at this point")
  }
  
  diag_ests <- outer(diag(out$A), diag(out$B), "+")
  nits <- NA
  if(type != "sparse"){
    nits <- out$nit
  }
  out$A <- t(out$A)   # legacy transpose
  
  fit <- list(call=match.call(), 
              A=out$A, B=out$B, beta=out$beta, diagAB=diag_ests,
              fitted.values=out$Yhat, 
              residuals=Y-out$Yhat, 
              D=D, X=X, Y=Y, 
              Xin=Xin, Yin=Yin, lag=lag,
              use_cov=use_cov, type=type, nit=nits,
              cv.glmnet.object=glmnetobj)
  class(fit) <- "blin"
  fit$R2 <- 1 - mean( (fit$residuals)^2, na.rm=TRUE) / mean(fit$Y^2, na.rm=TRUE)
  
  if(calcses & type=="full"){
    theta <- c(c(out$A), c(out$B)) 
    if(use_cov){
      theta <- c(theta, out$beta)
    }
    fit$vcov <- NA
    temp <- vcov(fit)
    fit$vcov <- temp
    fit$pval <- 1 - pnorm(abs(theta) / sqrt(diag(temp)))
    
    adiag <- seq(1, S^2, by=S+1)
    bdiag <- seq(1, L^2, by=L+1) + S^2
    se_diag <- matrix(0, S, L)
    for(i in 1:length(adiag)){
      for(j in 1:length(bdiag)){
        se_diag[i,j] <- sqrt(temp[adiag[i],adiag[i]] + temp[bdiag[j],bdiag[j]] + 2*temp[adiag[i],bdiag[j]])
      }
    }
    
    fit$se_diag <- se_diag
    fit$pval_diag <- 1 - pnorm(abs(diag_ests) / se_diag)
    
  } else {
    if(type != "full" & calcses){
      warning("Standard error computation for non-full BLIN model fits not implemented; returning NA")
    }
    fit$vcov <- NA
    fit$se_diag <- NA
    fit$pval <- fit$pval_diag <- NA
  }
  
  return(fit)
}


##################
###  Generics  ###
##################


#' vcov S3 generic for class blin
#' @export
#' @param object blin object
#' @param ... ignored
vcov.blin <- function(object, ...)
{
  if(is.na(object$vcov[1]) & object$type == "full"){
    Xreg <- model.matrix(object)
    S <- nrow(object$A)
    L <- nrow(object$B)
    if(object$use_cov){
      p <- dim(object$X)[4]
    } else {
      p <- 0
    }
    df <- S^2 + L^2 - 1 + p
    n <- sum(!is.na(object$Y))
    mse <-  sum(object$residuals^2, na.rm=T) / (n - df)
    vout <- mse*ginv( as.matrix( Matrix::crossprod(Xreg) ) )
  } else if (!is.na(object$vcov[1]) & object$type == "full"){
    vout <- object$vcov
  } else if (object$type != "full"){
    warning("Variance computation for non-full BLIN model fits not implemented; returning NA")
    vout <- NA
  }
  return(vout)
}


#' Coef S3 generic for class blin
#' @export
#' @param object blin object
#' @param whichcoef optional string (or NULL) indicating which coefficient to retrun, i.e. A, B, beta, or diagAB. If NULL, returns list of all coefficients.
#' @param ... ignored
coef.blin <- function(object, whichcoef=NULL, ...)
{
  if(is.null(whichcoef)){
    out <- list(A=object$A, B=object$B, beta=object$beta, diagAB=object$diagAB)
  } else {
    if(!is.character(whichcoef)){stop("Must specify character or NULL for 'whichcoef' ")}
    if (strtrim(whichcoef, 1) == "A"){
      out <- (object$A)
    } else if (strtrim(whichcoef, 1) == "B"){
      out <- (object$B)
    } else if (strtrim(whichcoef, 2) == "be"){
      out <- (object$beta)
    } else if (strtrim(whichcoef, 1) == "d"){
      out <- (object$diagAB)
    }
  }
  return(out)
}


#' Print S3 generic for class blin
#' @export
#' @param x blin object
#' @param hn optional numeric length of each coefficient printed
#' @param ... ignored
print.blin <- function(x, hn=10,...)
{
  if(!is.numeric(hn)){stop("Numeric input required for hn")}
  cat("\nCall: \n")
  print(x$call)
  cat("\nbeta coefficients:\n")
  print(head(x$beta, hn))
  cat("\n ... ")
  cat("\nA coefficients:\n")
  print(head(x$A, hn))
  cat("\n ... ")
  cat("\nB coefficients:\n")
  print(head(x$B, hn))
  cat("\n ... ")
  cat("\ndiag coefficients:\n")
  print(head(x$diagAB, hn))
  cat("\n ... ")
}


#' Plot S3 generic for class blin
#' @export
#' @param x blin object
#' @param ... ignored
plot.blin <- function(x, ...)
{
  hist(scale(resid(x)), freq=F, xlab="standardized residuals", main="")
  
  plot(fitted.values(x), scale(resid(x)), xlab="fitted values", ylab="standardized residuals", main="")
  
  qqnorm(scale(resid(x)), main="Normal Q-Q Plot for residuals")
  abline(0,1, col="red")
}


#' model.matrix S3 generic for class blin
#' @export
#' @param object blin object
#' @param ... ignored
model.matrix.blin <- function(object, ...)
{
  S <- dim(object$D)[1]
  L <- dim(object$D)[2]
  if(!is.null(object$X)){
    p <- dim(object$X)[4]
  } else {
    p <- 0
  }
  sizeinGB <- prod(dim(object$D))*(S^2 + L^2 + p)*8/1e9
  if(sizeinGB > 1){
    warning("model.matrix.blin() may take awhile / lots of memory for large S, L")
    cat("\n estimated model matrix size ", sizeinGB, "GB \n")
  }
  Xreg <- build_design(Y=object$Yin, X=object$Xin, lag=object$lag)  # build design matrix
  return(Xreg)
}


#' Summary S3 generic for class blin
#' @export
#' @param object blin object
#' @param whichcoef optional string (or NULL) indicating which coefficient to retrun, i.e. A, B, beta, or diagAB. If NULL, returns list of all coefficients.
#' @param ... ignored
summary.blin <- function(object, whichcoef=NULL, ...)
{
  secompute <- !is.na(object$vcov[1])
  if(!secompute){
    warning("Standard errors not yet computed; nothing much to summarize")
    out <- list(A=object$A, B=object$B, beta=object$beta, diagAB=object$diagAB, call=object$call)
  }
  
  if(is.null(whichcoef)){
    if(secompute){
      vc <- object$vcov
      pv <- object$pval
      S <- nrow(object$A)
      L <- nrow(object$B)
      p <- length(object$beta)
      
      ia <- 1:S^2
      acoef <- cbind(c(object$A), sqrt(diag(vc)[ia]), c(object$A) / sqrt(diag(vc)[ia]), pv[ia])
      ib <- S^2 + 1:(L^2)
      bcoef <- cbind(c(object$B), sqrt(diag(vc)[ib]), c(object$B) / sqrt(diag(vc)[ib]), pv[ib])
      ibeta <- S^2 + (L^2) + 1:p
      betacoef <- cbind(c(object$beta), sqrt(diag(vc)[ibeta]), c(object$beta) / sqrt(diag(vc)[ibeta]), pv[ibeta])
      diagcoef <- cbind(c(object$diagAB), c(object$se_diag), c(object$diagAB) / c(object$se_diag), 
                        c(object$pval_diag))
      
      colnames(acoef) <- colnames(bcoef) <- colnames(betacoef) <- 
        colnames(diagcoef) <- c("Estimate", "Std. Error", "t value", "Pr(|t| > 0)")             
      
      out <- list(A=acoef, B=bcoef, beta=betacoef, diagAB=diagcoef, call=object$call)
    } 
    
  } else {
    if(!is.character(whichcoef)){stop("Must specify character or NULL for 'whichcoef' ")}
    
    if (strtrim(whichcoef, 1) == "A"){
      ia <- 1:S^2
      acoef <- cbind(c(object$A), sqrt(diag(vc)[ia]), c(object$A) / sqrt(diag(vc)[ia]), pv[ia])
      colnames(acoef) <-  c("Estimate", "Std. Error", "t value", "Pr(|t| > 0)")
      out <- list(A=acoef, call=object$call)
      
    } else if (strtrim(whichcoef, 1) == "B"){
      ib <- S^2 + 1:(L^2)
      bcoef <- cbind(c(object$B), sqrt(diag(vc)[ib]), c(object$B) / sqrt(diag(vc)[ib]), pv[ib])
      colnames(bcoef) <-  c("Estimate", "Std. Error", "t value", "Pr(|t| > 0)")
      out <- list(B=bcoef, call=object$call)
      

    } else if (strtrim(whichcoef, 2) == "be"){
      ibeta <- S^2 + (L^2) + 1:p
      betacoef <- cbind(c(object$beta), sqrt(diag(vc)[ibeta]), c(object$beta) / sqrt(diag(vc)[ibeta]), pv[ibeta])
      colnames(betacoef) <-  c("Estimate", "Std. Error", "t value", "Pr(|t| > 0)")
      out <- list(beta=betacoef, call=object$call)
    
    } else if (strtrim(whichcoef, 1) == "d"){
      diagcoef <- cbind(c(object$diagAB), object$se_diag, c(object$diagAB) / object$sediag, 
                        object$pval_diag)
      colnames(diagcoef) <-  c("Estimate", "Std. Error", "t value", "Pr(|t| > 0)")
      out <- list(diagAB=betacoef, call=object$call)    
    }
  }
  
  class(out) <- "summary.blin"
  return(out)
}


#' Print S3 generic for class summary.blin
#' @export
#' @param x summary.blin object
#' @param hn optional numeric length of each coefficient printed
#' @param ... ignored
print.summary.blin <- function(x, hn=10,...)
{
  if(!is.numeric(hn)){stop("Numeric input required for hn")}
  cat("\nCall: \n")
  print(x$call)
  
  for(i in 1:(length(x) - 1)){
    cat("\n",  names(x)[[i]], "coefficients:\n")
    printCoefmat(head(x[[i]], hn))
    cat("\n ... ")
  }
}

