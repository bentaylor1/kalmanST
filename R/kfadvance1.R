

##' kfadvance_stationary_STGP function
##'
##' A function to 
##'
##' @param obs X 
##' @param oldmean X 
##' @param oldvar X 
##' @param A X 
##' @param B X 
##' @param C X 
##' @param D X 
##' @param E X 
##' @param F X 
##' @param W X 
##' @param V X 
##' @param marglik X 
##' @param log X 
##' @param na.rm X 
##' @return ...
##' @export

kfadvance_stationary_STGP <- function (obs, oldmean, oldvar, A, B, C, D, E, F, W, V, marglik = FALSE,log = TRUE, na.rm = FALSE){
    if (na.rm) {
        if (any(is.na(obs))) {
            if (all(is.na(obs))) {
                if (log) {
                  return(list(mean = diag(A)*oldmean + B, var = t(t(diag(A)*oldvar)*diag(A)) + W, mlik = 0))
                }
                else {
                  return(list(mean = diag(A)*oldmean + B, var = t(t(diag(A)*oldvar)*diag(A)) + W, mlik = 1))
                }
            }
            else {
                #M <- diag(length(obs))
                ind <- !is.na(obs)
                #M <- M[-which(is.na(obs)), ]
                obs <- obs[ind]
                #D <- M %*% D
                D <- D[ind,,drop=FALSE]
                V <- diag(V[1,1],sum(ind))
            }
        }
    }
    T <- diag(A)*oldmean
    S <- t(t(diag(A)*oldvar)*diag(A)) + W
    thing1 <- D %*% S
    tD <- t(D)
    K <- thing1 %*% tD + V
    margmean <- D %*% T
    resid <- obs - margmean
    if (marglik == TRUE) {
        if (all(dim(K) == 1)) {
            thing2 <- S %*% tD
            newmean <- T + as.numeric(1/K) * thing2 %*% resid
            newvar <- S - as.numeric(1/K) * thing2 %*% thing1
            marginal <- dnorm(as.numeric(obs), as.numeric(margmean), sqrt(as.numeric(K)), 
                log = log)
        }
        else {
            Kchol <- chol(K)
            Kcholinv <- solve(Kchol)
            logdetK <- 2*sum(log(diag(Kchol)))
            Kinv <- Kcholinv%*%t(Kcholinv)
            #Kinv <- solve(K)
            thing3 <- tD %*% Kinv
            thing4 <- S %*% thing3
            newmean <- T + thing4 %*% resid
            newvar <- S - thing4 %*% thing1
            #if(!isSymmetric(newvar)){browser()}
            marginal <- -(1/2)*logdetK + (-1/2) * t(resid) %*% Kinv %*% resid
            #marginal <- dmvnorm(as.vector(obs),as.vector(margmean),K,log=TRUE)
            if (!log) {
                marginal <- exp(marginal)
            }
        }
        return(list(mean = newmean, var = newvar, mlik = marginal))
    }
    else {
        if (all(dim(K) == 1)) {
            thing2 <- S %*% tD
            newmean <- T + as.numeric(1/K) * thing2 %*% resid
            newvar <- S - as.numeric(1/K) * thing2 %*% thing1
        }
        else {
            Kinv <- solve(K)
            thing3 <- tD %*% Kinv
            thing4 <- S %*% thing3
            newmean <- T + thing4 %*% resid
            newvar <- S - thing4 %*% thing1
        }
        return(list(mean = newmean, var = newvar))
    }
}