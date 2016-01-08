##' kffit_stationary_STGP function
##'
##' A function to 
##'
##' @param form X 
##' @param param list of parameters, see ?paramSpec
##' @param pid X 
##' @param tid X 
##' @param tmax X
##' @param data X 
##' @param cov.model X 
##' @param history.means X 
##' @param history.vars X 
##' @param diagonal_only X 
##' @param prior.mean X 
##' @param prior.var X 
##' @param fit X 
##' @param se.fit X 
##' @param se.predict X 
##' @param noisy X 
##' @param na.rm X 
##' @return ...
##' @export


kffit_stationary_STGP <- function(  form,
                                    param,
                                    pid,
                                    tid,
                                    tmax,
                                    data,
                                    cov.model,
                                    history.means=TRUE,
                                    history.vars=TRUE,
                                    diagonal_only=TRUE,
                                    prior.mean=NULL,
                                    prior.var=NULL,
                                    fit=TRUE,
                                    se.fit=TRUE,
                                    se.predict=TRUE,
                                    noisy=TRUE,
                                    na.rm=TRUE){

    if(!(inherits(form,"formula") | inherits(form,"character"))){
        stop("form must be an object inheriting class formula, or character")
    }
    
    if(is.character(form)){
        Nfixed <- 0
    }
    else{
        modelmat <- model.matrix(form,data)
        Nfixed <- ncol(modelmat)
    }


    uqid <- sort(unique(data[[pid]]))
    crds <- coordinates(data)
    uqcrds <- t(sapply(uqid,function(x){crds[which(data[[pid]]==x)[1],]}))

    kf_spec <- list()
    kf_spec$Nspatial <- length(uqid)
    kf_spec$Nfixed <- Nfixed
    kf_spec$T <- tmax

    a = param$a
    sigma_S = param$sigma_S
    phi_S = param$phi_S
    sigma_Y = param$sigma_Y

    npar <- kf_spec$Nspatial + kf_spec$Nfixed 

    u <- as.vector(as.matrix(dist(uqcrds)))
    Sigma <- matrix(EvalCov(cov.model,u=u,parameters=c(sigma_S,phi_S)),kf_spec$Nspatial,kf_spec$Nspatial)

    W <- matrix(0,npar,npar)
    W[1:kf_spec$Nspatial,1:kf_spec$Nspatial] <- Sigma
    W <- (1-a^2)*W

    nobs <- sapply(1:tmax,function(x){sum(data[[tid]]==x)})

    obser <- matrix(NA,kf_spec$Nspatial,tmax)

    Xt <- list()
    if(Nfixed>0){
        mm <- model.frame(form,data,na.action=function(x){x}) # model frame because it deals with NA's in the way I want
        mm <- mm[,-1] # remove outcome column
        mm <- model.matrix(form,mm)
    }

    for(t in 1:tmax){
        if(sum(data[[tid]]==t)==0){
            next
        }
        else{
            if(Nfixed>0){
                ob <- data[[as.character(form[2])]][data[[tid]]==t]
            }
            else{
                ob <- data[[form]][data[[tid]]==t]   
            }
            id <- data[[pid]][data[[tid]]==t]
            mch <- match(id,uqid) 
            obser[mch,t] <- ob
            if(Nfixed>0){ 
                Xt[[t]] <- matrix(NA,kf_spec$Nspatial,kf_spec$Nfixed)
                Xt[[t]][mch,] <- mm[data[[tid]]==t,]
            }
        }
    }

    Amat <- diag(rep(c(a,1),c(kf_spec$Nspatial,kf_spec$Nfixed)))
    Bmat <- matrix(0,npar,1)
    Cmat <- diag(1,npar)
    Emat <- matrix(0,kf_spec$Nspatial,1)
    Fmat <- diag(1,kf_spec$Nspatial)
    Vmat <- diag(sigma_Y^2,kf_spec$Nspatial)

    kf_spec$A <- function(t){return(Amat)}
    kf_spec$B <- function(t){return(Bmat)}
    kf_spec$C <- function(t){return(Cmat)}
    kf_spec$W <- function(t){return(W)}
    if(Nfixed>0){
        kf_spec$D <- function(t){return(cbind(diag(1,kf_spec$Nspatial),Xt[[t]]))}
    }
    else{
        kf_spec$D <- function(t){return(diag(1,kf_spec$Nspatial))}   
    }
    kf_spec$E <- function(t){return(Emat)}
    kf_spec$F <- function(t){return(Fmat)}
    kf_spec$V <- function(t){return(Vmat)}
    kf_spec$Y <- function(t){return(matrix(obser[,t],kf_spec$Nspatial,1))}

    out <- run_filter(kf_spec=kf_spec,
                        optim=FALSE,        # set this to FALSE if not estimating parameters, otherwise, if parameters have already been estimated, then set to TRUE 
                        history.means=history.means,# whether to save a matrix of filtered E(Theta_t) 
                        history.vars=history.vars, # whether to save a list of filtered V(Theta_t) 
                        diagonal_only=diagonal_only,
                        prior.mean=prior.mean,    # optional prior mean. column vector.    # Normally set  
                        prior.var=prior.var,     # optional prior mean. matrix.           # these inside the code 
                        fit=fit,          # whether to return a matrix of fitted values 
                        se.fit=se.fit,       # whether to return the standard error of the fitted values 
                        se.predict=se.predict,   # whether to return the prediction standard error = se(fitted values) + observation variance V_t 
                        noisy=noisy,         # whether to print a progress bar, useful. 
                        na.rm=na.rm)

    retlist <- list()
    retlist$Nspatial <- kf_spec$Nspatial 
    retlist$Nfixed <- kf_spec$Nfixed
    retlist$T <- tmax 
    retlist$form <- form
    retlist$param <- param
    retlist$id <- id
    retlist$data <- data
    retlist$diagonal_only <- diagonal_only
    retlist$obser <- obser
    retlist <- c(retlist,out)

    return(retlist)
}

##' run_filter function
##'
##' A function to 
##'
##' @param kf_spec X
##' @param optim X 
##' @param history.means X 
##' @param history.vars X 
##' @param diagonal_only X 
##' @param prior.mean X 
##' @param prior.var X 
##' @param fit X 
##' @param se.fit X 
##' @param se.predict X 
##' @param noisy X 
##' @param na.rm X 
##' @return ...
##' @export

 
run_filter <- function(kf_spec,
                    optim=FALSE,        # set this to FALSE if not estimating parameters, otherwise, if parameters have already been estimated, then set to TRUE 
                    history.means=FALSE,# whether to save a matrix of filtered E(Theta_t) 
                    history.vars=FALSE, # whether to save a list of filtered V(Theta_t) 
                    diagonal_only=TRUE,
                    prior.mean=NULL,    # optional prior mean. column vector.    # Normally set  
                    prior.var=NULL,     # optional prior mean. matrix.           # these inside the code 
                    fit=FALSE,          # whether to return a matrix of fitted values 
                    se.fit=FALSE,       # whether to return the standard error of the fitted values 
                    se.predict=FALSE,   # whether to return the prediction standard error = se(fitted values) + observation variance V_t 
                    noisy=TRUE,         # whether to print a progress bar, useful. 
                    na.rm=FALSE){       # whether to use NA handling. set to TRUE if any Y is NA  
           
    start <- Sys.time()   
     
    if(se.predict & !se.fit){ 
        stop("Must have se.fit=TRUE in order to compute se.predict") # leads to a computational saving  
    }     

    npar <- kf_spec$Nspatial + kf_spec$Nfixed     
     
    if(is.null(prior.mean)){ 
        Xpost <- rep(0,npar)
    } 
    else{ 
        Xpost <- prior.mean 
    } 
     
    if(is.null(prior.var)){    
        Vpost <- diag(1000,npar) 
    } 
    else{ 
        Vpost <- prior.var 
    }     
     
    if (history.means){ 
        Xrec <- matrix(NA,length(Xpost),kf_spec$T+1)
        Xrec[,1] <- Xpost 
    } 
    if(history.vars){ 
        if(diagonal_only){
            Vrec <- matrix(NA,length(Xpost),kf_spec$T+1)
            Vrec[,1] <- diag(Vpost)
        }
        else{
            Vrec <- list() 
            Vrec[[1]] <- Vpost 
        } 
    } 
    loglik <- 0  
    fitmat <- matrix(NA,kf_spec$Nspatial,kf_spec$T)  
    sefitmat <- matrix(NA,kf_spec$Nspatial,kf_spec$T)
    sepredictmat <- matrix(NA,kf_spec$Nspatial,kf_spec$T)
     
    if(noisy){ 
        pb <- txtProgressBar(min=1,max=kf_spec$T,style=3) 
    } 
       
    for(t in 1:kf_spec$T){ 
     
        # delete or complete the following rows as necessary 
        A <- kf_spec$A(t)
        B <- kf_spec$B(t)
        C <- kf_spec$C(t)
        W <- kf_spec$W(t)
        D <- kf_spec$D(t)
        E <- kf_spec$E(t)
        F <- kf_spec$F(t)
        V <- kf_spec$V(t)
        obs <- kf_spec$Y(t)


         
        # this bit calls KF advance       
        new <- kfadvance_stationary_STGP(obs=obs,oldmean=Xpost,oldvar=Vpost,A=A,B=B,C=C,D=D,E=E,F=F,W=W,V=V,marglik=TRUE,log=TRUE,na.rm=na.rm) 
         
        Xpost <- new$mean 
        Vpost <- new$var 
         
        if(t==1){ # used when this function is called iteratively one step at a time 
            running.mean <- Xpost 
            running.var <- Vpost 
        } 
         
        loglik <- loglik + new$mlik 
             
         
        if (history.means){ 
            Xrec[,t+1] <- Xpost 
        } 
        if(history.vars){ 
            if(diagonal_only){
                Vrec[,t+1] <- diag(Vpost)
            }
            else{
                Vrec[[t+1]] <- Vpost # since first entry is the prior 
            }
        } 

        if(fit){ 
            fitmat[,t] <- D%*%Xpost + E
        } 
        if(se.fit){ 
            sefitmat[,t] <- sqrt(diag(D%*%Vpost%*%t(D)))
        } 
        if(se.predict){ 
            sepredictmat[,t] <- sqrt(sefitmat[,t]^2+diag(F%*%V%*%t(F)))
        } 
         
        if(noisy){ 
            setTxtProgressBar(pb,t)  
        } 
    } 
     
    if(noisy){         
        close(pb) 
    } 
    end <- Sys.time() 
     
    if(noisy){         
        cat("Time taken:",difftime(end,start,units="secs")," seconds.\n") 
    } 
     
    if(optim){ 
        return(-loglik) # just return the -log likelihood if in parameter estimation mode 
    } 
    else{ 
        retlist <- list(mean=Xpost,var=Vpost,mlik=loglik,data=data,running.mean=running.mean,running.var=running.var) 
        if(history.means){ 
            retlist$history.means <- Xrec 
        } 
        if(history.vars){ 
            retlist$history.vars <- Vrec 
        }  
        if(fit){ 
            retlist$fit <- fitmat 
        } 
        if(se.fit){ 
            retlist$se.fit <- sefitmat 
        } 
        if(se.predict){ 
            retlist$se.predict <- sepredictmat         
        } 
        return(retlist) 
    } 
} 
 

