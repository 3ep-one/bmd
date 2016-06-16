if (FALSE){

source("~/git/bmd/proastUI/R/bmdFunctions.R")

}

library(proastUI)
library(testthat)


dataDir <- system.file("extdata", package = "proastUI")


track <- FALSE



context("Proast - Continuous summary data (atra)")

atra <- f.scan(file.path(dataDir, "atra.txt"))

test_that("Load data", {
      
      expect_match(atra$info, "atrazine")
      expect_equal(names(atra), c("info", "nvar", "varnames", "data", "dtype"))
      
    })

shinyInput <- list(dtype = 10, quick.ans = 1, 
	main.ans = 4, 
	# continous model
    model.ans = 1, 
	# dose
	xans = 1, 
	# mean response
	yans = 2, 
	# associated sd
	sans = 3, 
	# sample size
	nans = 4, 
	# extra parameters
	sd.se = 1, 
	# critical effect size
	CES = 0.05, 
	# is model continous
	cont = TRUE)


#test_that("Run f.proast() - initial run for all continuous models", {
#      
#      result <- lapply(c(1, 11, 13, 15, 23, 25), function(x){
#            
#            shinyInput$model.ans <- x
#            f.proast(odt = atra, shinyInput = shinyInput)
#            
#          })
#      
#      expect_equal(sapply(result, function(x) signif(x$loglik.first, 5)),
#          c(-12.7490, 4.6821, 3.7955, -8.6963, -17.2330, -16.9170)) # 14/06/2016: two last values modified
#      
#    })




test_that("Run f.proast() - fit model for all continuous models", {
      
      shinyInput$main.ans <- 4
      result <- lapply(c(1, 11, 13, 15, 23, 25), function(x){
            
            shinyInput$model.ans <- x
            f.proast(odt = atra, shinyInput = shinyInput)
            
          })
      
      expect_equal(sapply(result, function(x) x$npar), 
          c(2, 6, 4, 5, 4, 5))
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-11.16, 4.73, 4.15, 4.27, 4.16, 4.27))
      
      expect_equal(sapply(result, function(x) signif(x$MLE, 5)),
          list(c(0.068009, 437.790000),
              c(0.05495, 485.11000, 459.57000, 468.54000, 388.52000, 352.22000),
              c(0.055379, 480.020000, 3.545400, 0.688400),
              c(0.055286, 473.730000, 8.348500, 0.732930, 1.535500),
              c(0.055367, 479.100000, 3.984300, 0.771080),
              c(0.055288, 473.450000, 8.251600, 0.687830, 1.763600)))
      
    })


test_that("Run f.proast() - calculate CED with CI for Exp and Hill models", {
      
      shinyInput$main.ans <- c(4, 6)
      shinyInput$conf.lev <- 0.9
      
      result <- lapply(c(13, 15, 23, 25), function(x){
            
            shinyInput$model.ans <- x
            f.proast(odt = atra, shinyInput = shinyInput)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$CED, 5)),
          c(3.5454, 8.3485, 3.9843, 8.2516))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(0.17464, 0.20767, 0.22209, 0.26154,
                  14.665, 26.745, 15.061, 22.382), nrow = 2, byrow = TRUE))
      
      
      tmp <- f.plot.all(result[[4]])
      
      
    })




if(FALSE){

  ## Using summary function in performAnalysis.R
  shinyInput <- list(dtype = 10, xans = 1, yans = 2, nans = 4, sans = 3, 
      sd.se = 1, CES = 0.05, cont = TRUE, conf.lev = 0.9)
  
  fittedModels <- fitModels(data = atra, shinyInput = shinyInput, 
      continuousResponse = TRUE)
  toPrint <- summaryModels(fittedModels)
  f.plot.all(fittedModels[[attr(toPrint, "bestModelIndex")]])
  
}
