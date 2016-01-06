##' kfadvance function
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

kfadvance <- function (obs, oldmean, oldvar, A, B, C, D, E, F, W, V, marglik = FALSE,log = TRUE, na.rm = FALSE){
    if (na.rm) {
        if (any(is.na(obs))) {
            if (all(is.na(obs))) {
                if (log) {
                  return(list(mean = A %*% oldmean + B, var = A %*% 
                    oldvar %*% t(A) + C %*% W %*% t(C), mlik = 0))
                }
                else {
                  return(list(mean = A %*% oldmean + B, var = A %*% 
                    oldvar %*% t(A) + C %*% W %*% t(C), mlik = 1))
                }
            }
            else {
                M <- diag(length(obs))
                M <- M[-which(is.na(obs)), ]
                obs <- obs[which(!is.na(obs))]
                D <- M %*% D
                E <- M %*% E
                F <- M %*% F
            }
        }
    }
    T <- A %*% oldmean + B
    S <- A %*% oldvar %*% t(A) + C %*% W %*% t(C)
    thing1 <- D %*% S
    tD <- t(D)
    K <- thing1 %*% tD + F %*% V %*% t(F)
    margmean <- D %*% T + E
    resid <- obs - margmean
    if (marglik == TRUE) {
        if (all(dim(K) == 1)) {
            thing2 <- S %*% tD
            newmean <- T + as.numeric(1/K) * thing2 %*% resid
            newvar <- S - as.numeric(1/K) * thing2 %*% thing1
            marginal <- dnorm(obs, as.numeric(margmean), sqrt(as.numeric(K)), 
                log = log)
        }
        else {
            Kinv <- solve(K)
            thing3 <- tD %*% Kinv
            thing4 <- S %*% thing3
            newmean <- T + thing4 %*% resid
            newvar <- S - thing4 %*% thing1
            marginal <- (-1/2) * t(resid) %*% Kinv %*% resid
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
