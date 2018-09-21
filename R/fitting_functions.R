#' Fit full BLIN model
#' D is matrix of lagged Y (covariates)
#' Y is oservations
#' X is covariates
#' for multiple time periods
#' Returns A, B, beta, betaOLS, Yhat, 2LL, 2LLinit
#'
#' @keywords internal
fit_MLE_array_additive <- function(Y, D, X, lag, type, verbose=FALSE, printout=1, init="I", sigma_init=1, tol=1e-8, use_cov=TRUE, maxit=1e3, randseed=NA, Yold=NULL)
{
  # Get sizes and check
  if(dim(Y)[3] != dim(D)[3]){  stop("Dimenions of Y and D don't match")}
  # if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  # if(length(dim(Y)) < 3){ stop("Y is not an array")}
  # 
  impute <- sum(is.na(Y)) > 0
  if(impute & is.null(Yold)){
    stop("Cannot impute without discarded entries in Y")
  }
  i_na <- which(is.na(Y))
  Yold[is.na(Yold)] <- 0
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  m <- dim(D)[1]
  n <- dim(D)[2]
  tmax <- dim(Y)[3]
  
  
  # Initialize
  if(strtrim(init,1) == "r" | strtrim(init,1) == "R") {
    if(verbose == T){
      cat("Randomly initializing A,B \n")
    }
    if(is.numeric(randseed)){set.seed(randseed)}
    
    A <- matrix(rnorm(S*m, 0, sigma_init), S, m)
    B <- matrix(rnorm(L*n, 0, sigma_init), L, n)
    
  } else if (strtrim(init, 1) == "I" | strtrim(init, 1) == "i") {
    if(verbose == T){
      cat("Initializing A, B as identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    A <- diag(S)
    B <- diag(L)
    
  } else { stop("Invalid initialization type")}
  
  
  Ainit <- A
  Binit <- B
  
  # Preliminary matrices to avoid recalculating
  if(type == "sadd"){
    Jl <- matrix(1,L,L)
    Js <- matrix(1,S,S)
    DDT <- L*Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, D[,,x] )))    # LDJD^T
    DTD <- S*Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% D[,,x])))    # SD^T J
  } else if (type == "biten") {
    Jl <- diag(L)
    Js <- diag(S)
    DDT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x])))    # D %*% t(D)
    DTD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x])))    # t(D) %*% D
  } else { stop( "Invalid type " ) }
  
  if(impute){
    Y[i_na] <- (amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2) )[i_na]   # initialize NAs
    data <- lag_Y(lag, abind(Yold, Y), X=NULL)   # update D if imputing
    D <- data$D
  }
  
  if(use_cov){
    betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
    betaOLS[which(is.na(betaOLS))] <- 0   # if not full rank, get NAs
  } else { betaOLS <- NA}  
  
  # Find optimal values
  change <- 100
  count <- 0
  # if(nper < 5){ warning(paste0("Only ", round(nper, 3), " data points per coefficient"))}
  
  
  while(change > tol & count < maxit){

    # Update for beta
    Ystar <- Y - amprod(amprod(D, A, 1), Jl, 2) - amprod(amprod(D, Js, 1), t(B), 2)
    if(use_cov){
      beta <- matrix(lm(c(Ystar) ~ -1 + t(mat(X, 4)) )$coef, nrow=1)
      beta[which(is.na(beta))] <- 0
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2, na.rm=TRUE) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Ytilde <- Y - Xbeta
    if(type == "sadd"){
      Jl <- matrix(1,L,L)
      Js <- matrix(1,S,S)
      DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, Ytilde[,,x])))    # D J Y^T
      DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% Ytilde[,,x])))    # D^T J Y
    } else if (type == "biten") {
      Jl <- diag(L)
      Js <- diag(S)
      DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], Ytilde[,,x])))    # D %*% t(Y)
      DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Ytilde[,,x])))    # t(D) %*% Y
    } else { stop( "Invalid type " ) }
    result <- update_MLE_additive(A, B, D, DDT, DTD, DYT, DTY, type)
    Anew <- result$A
    Bnew <- result$B
    changeAB <- max(abs(c(c(A - Anew), c(B-Bnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    B <- Bnew
    
    Yhat <- Y - Ytilde + amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2) 
    if(impute){
      Y[i_na] <- Yhat[i_na]   # update NAs
      data <- lag_Y(lag, abind(Yold, Y), X=NULL)   # update D if imputing
      D <- data$D
    }
    
    LLnew <- -length(Y)*log( sum( (Y - Yhat)^2, na.rm=TRUE) ) - length(Y)    #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    # cat("True 2LL*tau^2: \t", LLtrue, "\n")
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  Yhat <- Xbeta + amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2)
  
  return(list(A=A, B=B, Yhat=Yhat, beta= beta, betaOLS = betaOLS, LLt2 = LL, LLt2_init=LLinit, nit=count))
}





#' Fit reduced rank BLIN model
#' D is matrix of lagged Y (covariates)
#' Y is oservations
#' X is covariates
#' for multiple time periods
#' Returns A, B, beta, betaOLS, Yhat, 2LL, 2LLinit
#'
#' @keywords internal
fit_MLE_array <- function(Y, D, X, lag, rankA=1, rankB=1, verbose=FALSE, printout=1, tol=1e-8, init="I", sigma_init=1, use_cov=TRUE, maxit=1e3, randseed=NA, Yold=NULL)
{
  # Get sizes and check
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  # if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  # if(length(dim(Y)) < 3){ stop("Y is not an array")}
  impute <- sum(is.na(Y)) > 0
  if(impute & is.null(Yold)){
    stop("Cannot impute without discarded entries in Y")
  }
  i_na <- which(is.na(Y))  # save NAs
  Y[i_na] <- 0  # set to zero
  Yold[is.na(Yold)] <- 0
  
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  tmax <- dim(Y)[3]
  
  k <- rankA 
  m <- rankB
  nper <- length(Y)/(2*(S*k + L*m) + dim(X)[4])    # Warning for data size
  
  # Initialize
  if(strtrim(init,1) == "r" | strtrim(init,1) == "R") {
    if(verbose == TRUE){
      cat("Randomly initializing U,V,W,Z \n")
    }
    if(is.numeric(randseed)){set.seed(randseed)}
    U <- matrix(rnorm(S*k, 0, sigma_init), S, k)
    V <- matrix(rnorm(S*k, 0, sigma_init), S, k)
    W <- matrix(rnorm(L*m, 0, sigma_init), L, m)
    Z <- matrix(rnorm(L*m, 0, sigma_init), L, m)
    
    A <- tcrossprod(U,V)
    BT <- tcrossprod(Z,W)
    
  } else if (strtrim(init, 1) == "I" | strtrim(init, 1) == "i") {
    if(verbose == TRUE){
      cat("Initializing A, B as identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    U <- eigen(diag(S))$vectors[,1:rankA]
    V <- U
    W <- eigen(diag(L))$vectors[,1:rankB]
    Z <- U
    
    A <- diag(S)
    BT <- diag(L)
    
  } else { stop("Invalid initialization type")}
  
  
  
  Ainit <- A
  BTinit <- BT
  
  # Preliminary matrices to avoid recalculating
  DDT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x])))    # t(D) %*% D
  DTD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x])))    # D %*% t(D)
  
  if(use_cov){
    betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
    betaOLS[which(is.na(betaOLS))] <- 0   # if not full rank, get NAs
  } else { betaOLS <- NA}
  
  # Find optimal values
  change <- 100
  count <- 0
  # if(nper < 5){ warning(paste0("Only ", round(nper, 3), " data points per coefficient"))}
  
  
  while(change > tol & count < maxit){

    # Update for beta
    Ystar <- Y - amprod(D, A, 1) - amprod(D, BT, 2)
    if(use_cov){
      beta <- matrix(lm(c(Ystar) ~ -1 + t(mat(X, 4)) )$coef, nrow=1)
      beta[which(is.na(beta))] <- 0
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Ytilde <- Y - Xbeta
    DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], Ytilde[,,x])))    
    DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Ytilde[,,x])))    
    result <- update_MLE_asymmetric(D, U, V, W, Z, DDT, DTD, DYT, DTY)
    U <- result$U
    V <- result$V
    W <- result$W
    Z <- result$Z
    Anew <- tcrossprod(U,V)
    BTnew <- tcrossprod(Z,W)
    changeAB <- max(abs(c(c(A - Anew), c(BT-BTnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    BT <- BTnew
    
    Yhat <- Y - Ytilde + amprod(D, A, 1) + amprod(D, BT, 2)  # estimate
    if(impute){
      Y[i_na] <- Yhat[i_na]  # update NAs
      data <- lag_Y(lag, abind(Yold, Y), X=NULL)   # update D if imputing
      D <- data$D
    }
    
    # LLnew <- -sum( (Ytilde - (amprod(D, A, 1) + amprod(D, BT, 2)) )^2)
    LLnew <- -length(Y)*log( sum( (Y - Yhat )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    # cat("True 2LL*tau^2: \t", LLtrue, "\n")
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  # Post-process
  A <- tcrossprod(U,V)
  B <- tcrossprod(W,Z)
  
  Yhat <- Xbeta + amprod(D, A, 1) + amprod(D, t(B), 2)
  
  return(list(A=A, B=B, Yhat = Yhat, beta= beta, betaOLS = betaOLS, LLt2 = LL, LLt2_init=LLinit, nit=count))
}








#' Fit sparse BLIN model
#' D is matrix of lagged Y (covariates)
#' Y is oservations
#' X is covariates
#' for multiple time periods
#' Returns A, B, beta, Yhat
#'
#' @keywords internal
additive_regression <- function(Y, X, lag, type="biten", use_cov=TRUE, penalty=1, whichlambda="min")
{
  # require(glmnet)
  
  # # Find sizes, assuming dim(D) = dim(Y)
  # S <- nrow(D[,,1])
  # L <- ncol(D[,,1])
  # tmax <- dim(D)[3]
  # Find sizes, 
  S <- nrow(Y[,,1])
  L <- ncol(Y[,,1])
  tmax <- dim(Y)[3] - lag 
  
  # 
  # Build X matrix
  Xreg <- build_design(Y, X, lag=lag)  # build design matrix
  
  
  # Perform regression and pull out coefficients
  # fit <- lm(c(Y) ~ -1 + Xreg)
  if(is.numeric(penalty)){
    yfit <- c(Y[,,-c(1:lag)])
    keep <- which(!is.na(c(yfit)))   # glmnet cannot handle NAs
    cv.fit <- cv.glmnet(x=Xreg[keep,], y=yfit[keep], alpha=penalty, family='gaussian', intercept=FALSE)
    # cv.fit <- cv.glmnet(Xreg, y=as.factor(c(Y)), alpha=penalty, family='binomial', intercept=F)
    fit <- cv.fit$glmnet.fit
    if(whichlambda == "1sd"){
      col <- cv.fit$lambda == cv.fit$lambda.1se  # .min??
    } else if (whichlambda == "min"){
      col <- cv.fit$lambda == cv.fit$lambda.min  # .min??
    } else {
      warning("Invalid choice of lambda. Defaulting to lambda with minimum CV measure")
      col <- cv.fit$lambda == cv.fit$lambda.1se  # .min??
    }
    coefs <- fit$beta[, col]   # also can use <- coef(cv.fit, s="lambda.1se") # or "lambda.min"
    
    Yhat <- array((predict(fit, newx=Xreg, type="response")[, col]), c(S,L,tmax))    # need own predict method
    
  } else {
    stop("penalty input to lasso must be numeric")
    # if(S < 30 | strtrim(type, 3) == "bit"){   # fit.lm fast for small dimensions
    #   remove <- which(is.na(c(Y)))  # NAs to remove
    #   if(length(remove) > 0){
    #     x <- Xreg[-remove, ]
    #     y <- c(Y)[-remove]
    #   } else {
    #     x = Xreg  ;  y = c(Y)
    #   }
    #   fit <- lm.fit(x=x, y=y)   # faster than lm()
    #   coefs <- fit$coef
    #   Yhat <- array(NA, c(S,L,tmax))
    #   Yhat[!is.na(Y)] <- fit$fitted.values
    # } else {
    #   fit <- solve_lm(x=Xreg, y=c(Y))  # handles NAs internally
    #   coefs <- fit$coef
    #   Yhat <- array(fit$pred, c(S,L,tmax))
    # }
  }
  A <- matrix(coefs[1:S^2], nrow=S)
  B <- matrix(coefs[S^2 + 1:L^2], nrow=L)
  
  if(use_cov){
    p <- dim(X)[4]
    beta <- matrix(coefs[S^2 + L^2 + 1:p], ncol=p)
  } else { beta <- NA }
  
  return(list(A=A, B=B, beta=beta, Yhat=Yhat, Xreg=Xreg, cv.glmnet.object=cv.fit))
}






#' Update for block coordinate descent of full BLIN model
#' D is matrix of lagged Y (covariates)
#' Y is observations
#' 
#' @keywords internal
update_MLE_additive <- function(A, B, D, DDT, DTD, DYT, DTY, type="biten")
{
  S <- dim(D)[1]
  L <- dim(D)[2]
  tmax <- dim(D)[3]
  
  if(type == "sadd"){
    Jl <- matrix(1,L,L)
    Js <- matrix(1,S,S)
    DBD <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, Js %*% D[,,x] %*% B)))    # D J B^T D^T J
    A <- t(solve(DDT, DYT - DBD))
    DAD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% A %*% D[,,x] %*% Jl)))    # D^T J A D J
    B <- (solve(DTD, DTY - DAD))   # no transpose! Use B instead of B^T
    
  } else if (type == "biten") {
    DBD <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], D[,,x] %*% B)))    # D B^T D^T
    A <- t(solve(DDT, DYT - DBD))
    DAD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], A %*% D[,,x])))    # D^T A D 
    B <- (solve(DTD, DTY - DAD))  # no transpose! Use B instead of B^T
    
  } else { stop( "Invalid type " ) }
  
  return(list(A=A, B=B))
}



# Update for block coordinate descent of reduced rank BLIN model
#' D is matrix of lagged Y (covariates)
#' Y is observations
#' 
#' @keywords internal
update_MLE_asymmetric <- function(D, U, V, W, Z, DDT, DTD, DYT, DTY)
{
  # sizes
  m <- ncol(W)
  p <- nrow(W)
  k <- ncol(U)
  n <- nrow(U)
  
  t <- dim(D)[3]   # number of time slices
  
  # W and Z update
  DAD <- Reduce("+", lapply(1:t, function(x) crossprod(D[,,x ], tcrossprod(U, V) %*% D[,,x])  ))
  Z <- t( solve( crossprod(W, DTD %*% W), crossprod(W, DTY - DAD)))
  W <- solve(DTD, DTY - DAD) %*% Z %*% solve(crossprod(Z))
  
  # U and V update
  DBD <- Reduce("+", lapply(1:t, function(x) tcrossprod(D[,,x ] %*% tcrossprod(Z, W), D[,,x])  ))
  V <- solve(DDT, DYT - DBD) %*% U %*% solve(crossprod(U))
  U <- t(solve(crossprod(V, DDT %*% V), crossprod(V, DYT - DBD)))
  
  return(list(U=U, V=V, W=W, Z=Z))
}

