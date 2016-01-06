##' predictPlots function
##'
##' A function to 
##'
##' @param out X 
##' @param pause X
##' @param ... X
##' @return ...
##' @export

predictPlots <- function(out,pause=TRUE,...){
    n <- nrow(out$fit)
    T <- ncol(out$fit)

    for(i in 1:n){     
        upp <- out$fit[i,] + 1.96*out$se.predict[i,]
        low <- out$fit[i,] - 1.96*out$se.predict[i,]
        plot(1:T,out$fit[i,],xlim=c(1,T),ylim=c(min(low),max(upp)),type="l")
        points(1:T,out$obser[i,],pch=19)
        polygon(c(1:T,T:1),c(upp,rev(low)),col=rgb(0,0,1,alpha=0.2),border=NA)

        if(pause){
            cat ("Press a key to continue")
            line <- readline()
        }
    }
}


# ## fit function
# ##
# ## A function to 
# ##
# ## @param out X 
# ## @return ...
# ## @export

# fit <- function(out){
#     return(out$fit)
# }



# ## se.fit function
# ##
# ## A function to 
# ##
# ## @param out X 
# ## @return ...
# ## @export

# se.fit <- function(out){
#     return(out$se.fit)
# }


# ## se.predict function
# ##
# ## A function to 
# ##
# ## @param out X 
# ## @return ...
# ## @export

# se.predict <- function(out){
#     return(out$se.predict)
# }


# ## observations function
# ##
# ## A function to 
# ##
# ## @param out X 
# ## @return ...
# ## @export

# observations <- function(out){
#     return(out$obser)
# }