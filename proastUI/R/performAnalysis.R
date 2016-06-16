#' Fit all models for the continuous or quantal response
#' @param data list, as returned from the f.scan() function;
#' @param shinyInput list, all parameter values as defined in the shiny UI
#' @param continuousResponse boolean, indicating whether the response is a
#' continuous covariate or not
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with for each model a list with all results obtained during 
#' the analysis with f.proast(); attribute "modelNames" contains the full names
#' of all models that were fitted; if error occurs return error message 
#' @export
fitModels <- function(data, shinyInput, continuousResponse, track = FALSE) {
  
  fitModelSeries <- function(main.ans, modelIndices){
    
    currentShinyInput <- shinyInput
    currentShinyInput$main.ans <- main.ans
    
    result <- lapply(modelIndices, function(iModel){
          
          currentShinyInput$model.ans <- iModel
          f.proast(odt = data, shinyInput = currentShinyInput, track = track)
          
        })
    
    return(result)
    
  }
  
  tryCatch({
        
        if(continuousResponse){
          
          savedResults <- fitModelSeries(main.ans = 4, modelIndices = c(1, 11))
          savedResults <- c(savedResults, fitModelSeries(main.ans = c(4, 6), 
                  modelIndices = c(13, 15, 23, 25))
          )
          
          modelNames <- c("Null", "Full", "Exp model 3",
              "Exp model 5", "Hill model 3", "Hill model 5")
          
          # TODO quantal response
#          } else {
#            
#            fitModelSeries(main.ans = 4, modelIndices = c(1, 14))
#            fitModelSeries(main.ans = c(4, 6),
#                modelIndices = c(26, 25, 18, 21, 19, 24, 16))
#            modelNames <- c("Null", "Full", "Logistic", "Probit",
#                "Log-logistic", "Log-probit", "Weibull", "Gamma", "Two-stage")
#            bmd <- c(NA, NA, unlist(sapply(storedResults[3:length(storedResults)], 
#                        function(x) x$CED.matr)))
          
        }
        
        attr(savedResults, "modelNames") <- modelNames
        
        return(savedResults)
        
      }, error = function(err) {
        
        return(err)
        
      })   
  
}


#' Summarize the results of all fitted models for continuous or quantal response
#' @param savedResults list as returned by the function fitModels() 
#' @return data frame with for each model the model = model names, npar = number 
#' of parameters in the model, loglik = the estimated log likelihood, aic = the
#' estimated aic, bmd, bmdl and bmdu, respectively the estimate, lower and upper 
#' bound of the benchmark dose, converged = whether the fitting procedure 
#' ended with convergion or not; attributes minimum and maximum of all bmdl resp
#' bmdu values, the index and name of the best (min aic) model, the name of the 
#' response 
#' @export
summaryModels <- function(savedResults){
  
  bmd <- sapply(savedResults, function(x) x$CED)
  bmdConf <- sapply(savedResults, function(x) x$conf.int)
  
  returnTable <- data.frame(
      model = attr(savedResults, "modelNames"),
      npar = sapply(savedResults, function(x) x$npar),
      loglik = sapply(savedResults, function(x) x$loglik),
      aic = sapply(savedResults, function(x) 
            2 * x$npar - 2 * x$loglik),
      bmd = bmd,
      bmdl = bmdConf[1,],
      bmdu = bmdConf[2,],
      converged = as.logical(sapply(savedResults, function(x) x$converged)),
      row.names = NULL)
  
  aicNotFull <- min(returnTable[which(returnTable$model != 'Full'), 'aic'], na.rm = TRUE)
  returnTable$accepted <- (returnTable$aic <= (aicNotFull + 2))
  
  
  attr(returnTable, "minBmdl") <- min(returnTable$bmdl[returnTable$accepted], na.rm = TRUE)
  attr(returnTable, "maxBmdu") <- max(returnTable$bmdu[returnTable$accepted], na.rm = TRUE)
  attr(returnTable, "bestModelIndex") <- which.min(returnTable$aic)
  attr(returnTable, "bestModel") <- returnTable$model[which.min(returnTable$aic)]
  attr(returnTable, "responseName") <- tail(savedResults, n = 1)[[1]]$y.leg
  
  return(returnTable)
  
} 


