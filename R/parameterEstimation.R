##' paramSpec function
##'
##' A function to 
##'
##' @param a X 
##' @param sigma_S X 
##' @param phi_S X 
##' @param sigma_Y X 
##' @return ...
##' @export

paramSpec <- function(a=NULL,sigma_S=NULL,phi_S=NULL,sigma_Y=NULL){
    return(list(a=a,sigma_S=sigma_S,phi_S=phi_S,sigma_Y=sigma_Y))
}


## KFparest function
##
## A function to 
##
## @param data X 
## @param ... X 
## @return ...
## @export


# KFparest <- function(   data, # data ie Y  
#                         ...){ # delete and paste in OTHER ARGUMENTS TO BE PASSED TO KFfit 
#     start <- Sys.time() 
     
#     inits <- PUT INITIAL VALUES FOR PARAMETER VECTOR HERE 
     
#     # use optim to find optimal parameters 
#     oppars <- optim(inits, 
#                     KFfit, 
#                     data=data, 
#                     OTHER ARGUMENTS TO BE PASSED TO KFfit,  # delete and paste in OTHER ARGUMENTS  
#                                                             # TO BE PASSED TO KFfit 
#                     optim=TRUE, 
#                     control=list(trace=100))     
     
#     end <- Sys.time() 
#     cat("\n") 
#     cat("Time Taken",difftime(end,start,units="mins"),"\n") 
#     cat("\n") 
     
#     return(oppars) 
# } 