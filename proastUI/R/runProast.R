
#' Run the Proast Application
#' @return no return value
#' @import shiny
#' @importFrom devtools install_github
#' @export
runProast <- function(){
  
  if(!requireNamespace(assertive)){
    
    install.packages("assertive")
    
  }
  
  if(!requireNamespace(shinysky)){
    
    devtools::install_github("AnalytixWare/ShinySky")
    
  }
  
  
  proastTmpDir <- tempdir()
  
  setwd(proastTmpDir)
  
  # Copy server.R and ui.R (not folder www)
  proastUiDir <- system.file("ui", package = "proastUI")
  proastUiFiles <- list.files(path = proastUiDir, full.names = TRUE)
  proastUiFiles <- proastUiFiles[!grepl("www", proastUiFiles)]
  
  sapply(proastUiFiles, function(x){
        file.copy(from = x, to = file.path(proastTmpDir, basename(x)),
            overwrite = TRUE)}
  )
  
  # Make www directory and copy its files
  if (!dir.exists(file.path(proastTmpDir, "www"))) {
    
    dir.create(path = file.path(proastTmpDir, "www"))
    
  }
  
  wwwFiles <- list.files(path = file.path(proastUiDir, "www"), full.names = TRUE)
  
  sapply(wwwFiles, function(x){
        file.copy(from = x, to = file.path(proastTmpDir, "www", basename(x)),
            overwrite = TRUE)}
  )
  
  
  runApp(appDir = proastTmpDir)
  
}
