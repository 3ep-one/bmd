if(FALSE){
  
  library(shiny)
  runApp()
  
}

if (TRUE){

#  library(proastUI)
	tmp <- sapply(list.files("~/git/bmd/proastUI/R/", full.names = TRUE), source)
  
}

library(shinysky) # for hotable

`%then%` <- shiny:::`%OR%`
track <- FALSE

serverFunction <- function(input, output, session){
  
  # To debug
  output$print <- renderPrint({
        
        
        
      })
  
  
  
  results <- reactiveValues()
  
  
  # Load data
  results$dataLoaded <- reactive({
        
        if(input$example > 0){
          
          dataDir <- system.file("extdata", package = "proastUI")
          return(f.scan(file.path(dataDir, "atra.txt")))
          
#            load(path.expand(file.path(dataDir, "das11.rda")))
#            return(das11)
          
          
        } else {
          
          inFile <- input$dataLoaded
          
          if (is.null(inFile))
            return(NULL)
          
          tryCatch({ 
                
                if(input$dataType == "rawData"){
                  
                  rawData <- read.table(inFile$datapath, header = TRUE, sep = input$sep, 
                      quote = '"', stringsAsFactors = FALSE)
                  
                  covariatesName <- colnames(rawData)
                  covariatesType <- c()
                  
                  for(iName in covariatesName){
                    
                    covariatesType <- c(covariatesType, input[[paste(iName, "Type")]])
                    
                  }
                  
                  return(list(info = input$dataName, 
                          nvar = ncol(rawData),
                          varnames = covariatesName, 
                          data = rawData, 
                          dtype = as.numeric(covariatesType)))
                  
                  
                } else if(input$dataType == "proastData") {
                  
                  f.scan(inFile$datapath, separator = input$sep) 
                  # also returns list with info (title), nvar, varnames, data, dtype
                  
                } 
                
                
              },  error = function(err) {
                
                return(err)
                
              })
          
        }
        
      })
  
  
  
  output$dataLoaded <- DT::renderDataTable({
        
        validate(need(results$dataLoaded(), "No data loaded") %then%
                need( !any(attr(results$dataLoaded(), "class") == "error"), 
                    paste("Data were not loaded correctly: \n", results$dataLoaded()$message)))
        
#        if(is.null(results$dataModified())){
#          
          DT::datatable(results$dataLoaded()$data, options = list(
                  lengthMenu = list(c(15, 100, -1), c('15', '100', 'All')),
                  pageLength = 15), rownames = FALSE)
#          
#        } else {
#          
#          DT::datatable(results$dataModified()$data, options = list(
#                  lengthMenu = list(c(15, 100, -1), c('15', '100', 'All')),
#                  pageLength = 15), rownames = FALSE)
#          
#        }
        
      })
  
  
  observe({
        
        if(input$dataType == 'rawData'){
          
          inFile <- input$dataLoaded
          
          if(!is.null(inFile) & input$dataName == ""){
            
            newName <- inFile$name
            updateTextInput(session, "dataName", 
                value = substring(newName, first = 1, last = nchar(newName) - 4))
            
          }
          
        }
        
        if(input$example > 0){
          
          updateSelectInput(session, inputId = "dtype", selected = 10)
		  updateRadioButtons(session, inputId = "sep", 
				selected = '\t')
          
        }
        
      })
  
  
  
  output$covariatesType <- renderUI({ 
        
        covariatesType <- list()
        covariatesName <- (results$dataLoaded())$varnames
        
        if(length(covariatesName) < 1)
          return(NULL)
        
        
        for (i in seq_along(covariatesName)) {
          
          covariatesType[[i]] <- selectInput(paste(covariatesName[i], "Type"),
              label = covariatesName[i],
              choices = c("0: Non-response" = 0,
                  "1: Continuous response" = 1,
                  "2: Binary response" = 2,
                  "3: Ordinal response" = 3,
                  "4: Quantal response" = 4,
                  "5: Nested continuous response" = 5,
                  "6: Nested quantal response" = 6,
                  "10: Mean (continuous) response" = 10)
          )
          
        }
        
        list(
            h4("Specify data type for all variables"),
            covariatesType
        )
        
      }) 
  
  # Choices are named numbers
  results$varnamesChoices <- reactive({
        
        newNames <- (results$dataLoaded())$varnames
        values <- 1:length(newNames)
        names(values) <- newNames
        return(values)
        
      })
  
  
  
  output$CES <- renderUI({
        
        # TODO shinyInput$ces.ans for quantal response only
        if(input$dtype == '4'){
          
          selectInput("ces.ans", label = "What type of Benchmark response do you want to consider?",  
              choices = c("ED50" = 1, 
                  "Additional risk, i.e. P[BMD] - P[0]" = 2,
                  "Extra risk, i.e. (P[BMD]-P[0])/(1-P[0])" = 3, 
                  "CED for latent variable" = 4))
          
          #Switch text depending on shinyInput$ces.ans:
#       text <- switch(shinyInput$ces.ans, "", "Give value for the BMR, in terms of additional risk ",
#        "\nGive value for the BMR, in terms of extra risk", "Give value for the CES, defined for latent variable ")
          
        }
        
      })
  
  
  output$parameterQuestions <- renderUI({
        
        list(
            
            selectInput("xans",
                label = "Q1: Which variable do you want to consider as independent variable?",
                choices = results$varnamesChoices()),
            
            selectInput("Vyans",
                label = "Q2: Which response variable(s) do you want to analyse?",
                choices = results$varnamesChoices(), 
                selected = results$varnamesChoices()[2], multiple = TRUE)
        
#            conditionalPanel("input.dtype == 'continuous, clustered'",
#                selectInput("nest.no",
#                    label = "Select the nested factor",
#                    choices = c("<none>" = 0, results$varnamesChoices()))
#            ),
        # TODO this can be more than one factor: Adapt data when loading, not in proast
#            selectInput("select.no", label = "Q11a: For which factor do you want to select data",
#                choices = results$varnamesChoices())
        # TODO shinyInput$covar.no
#            selectInput("covar.no", 
#                label = "Which variable do you want to consider as potential covariate?",
#                choices = c("<none>" = "0", results$varnamesChoices()))
        
        )
      })
  
  
  output$parameterQuestions2 <- renderUI({
        
        createVnans <- function(iSelected){          
          
          yNames <- names(results$varnamesChoices())[as.numeric(input$Vyans)]
          Vnans <- list()
          
          for (i in 1:length(input$Vyans)) {
            
            Vnans[[i]] <- column(11, selectizeInput(paste0("nans", i),
                    label = yNames[i],
                    choices = results$varnamesChoices(),
                    selected = results$varnamesChoices()[iSelected]), offset = 1)
            
          }
          
          return( list(strong("Q3c: Give the associated sample size(s)"),
                  Vnans) )
          
        }
        
        validate(need(results$dtype(), "") %then%
                need(input$Vyans, ""))
        
        
        if(results$dtype() %in% c(10, 250, 260)){
          
          list(
              selectInput("sans", 
                  label = "Q3a: The variation statistic associated",
                  choices = results$varnamesChoices(), 
                  selected = results$varnamesChoices()[3]),
              radioButtons("sd.se",
                  label = "Q3b: Type of variation statistic associated with the means?",
                  choices = c("standard deviations" = 1, "standard errors" = 2)),
              createVnans(iSelected = 4)
          
          )
          
          
        } else if(!results$dtype() %in% c(1, 25, 26)){
          
          createVnans(iSelected = 3)
          
        }
        
      })
  
  output$transformation <- renderUI({
        
        if(input$dtype %in% c('1', '10')){
          
          selectInput("transformation", 
              "Transformation for the response variable(s)",
              choices = c("log", "square root", "<none>"))
          
        } else {
          
          NULL
          
        }
        
      })
  
  
  results$dataModified <- reactive({
        
        validate(need(results$dataLoaded(), "No data loaded"))
        
        
        # Remove missing values
        
        dataComplete <- f.remove.NAs(variableIndices = as.numeric(c(input$xans, 
                    input$Vyans, input$sans, input$covar.no, input$nans)), 
            originalData = results$dataLoaded()$data, track = track)
        
        # index of the variables used for parameters (a, b, variance (theta), c, d)
        allFactors <- c(input$fct1.no, input$fct2.no, input$fct3.no, 
            input$fct4.no, input$fct5.no)
        
        if(length(allFactors) != 0){
          
          dataComplete <- f.remove.NAs(variableIndices = allFactors, 
              originalData = dataComplete, track = track)
          
        }
        
        
        # Check whether x and y values are numeric vectors
        
        validate(need(input$xans, "Please provide answer to Q1") %then%
                need(input$Vyans, "Please provide answer to Q2"))
        
        f.check.nonneg.num(dataFrame = dataComplete[, as.numeric(input$xans)], track = track)
        isProblem <- f.check.nonneg.num(dataFrame = dataComplete[, as.numeric(input$Vyans)], track = track)
        
        validate(need(!any(isProblem), 
                "No log-transformation can be performed due to negative response values. 
                    \nPlease change the chosen transformation."))
        
        
        # Check for one level factors & provide warning
        
        if(length(allFactors) != 0){
          
          if (length(allFactors) == 1){
            
            oneLevel <- nlevels(dataComplete[, allFactors]) == 1
            
          } else {
            
            oneLevel <- apply(dataComplete[, allFactors], 2, 
                function(x) nlevels(x) == 1)
            
          }
          
          validate(need(!any(oneLevel), 
                  paste("The factor you chose as covariate on parameter(s)", 
                      paste(c("a", "b", "var", "c", "d")[oneLevel], collapse = ","),
                      "has only one level\n you might have selected a subgroup for this factor")))
          
        }
        
        xValues <- dataComplete[, as.numeric(input$xans)]
        xName <- names(results$varnamesChoices())[as.numeric(input$xans)]
        xRange <- range(xValues[xValues != 0])
        
        validate(need(xRange[2] <= 100, 
                paste("Warning: (Nonzero)", xName, "values range from", xRange[1],
                    "to ", xRange[2], "\nit is recommended to scale",
                    xName, "to prevent numerical problems")))
        
        
        list(info = (results$dataLoaded())$info, 
            nvar = (results$dataLoaded())$nvar,
            varnames = (results$dataLoaded())$varnames, 
            data = dataComplete, 
            dtype = (results$dataLoaded())$dtype)        
        
      })
  

  
      # TODO Change default settings
        
#  output$parameterConstraints <- renderHotable({
#        
#        results$parameterConstraints()
#        
#      }, readOnly = c(TRUE, FALSE, FALSE))
  
  
  
  # Define input parameters for f.proast()
  results$dtype <- reactive({
        
        if(is.null(input$dtype))
          return(NULL)
        
        currentDtype <- as.numeric(input$dtype)
        
        
        #No log-transform of response
        if(!is.null(input$transformation)){
          
          switch(input$transformation,
              
              "<none>" = {
                
                currentDtype[currentDtype == 1] <- 25
                currentDtype[currentDtype == 10] <- 250
                
              },
              
              "square root" = {
                
                currentDtype[currentDtype == 1] <- 26
                currentDtype[currentDtype == 10] <- 260
              })
          
        }
        
        return(currentDtype)
        
      })
  
  
  results$cont <- reactive({
        
        if(results$dtype() %in% c(1, 5, 10, 15, 25, 250, 26, 260)){
          
          TRUE
          
        } else {
          
          FALSE
          
        }
        
        
      })
  
  
  
  results$shinyInput <- reactive({
        
        if(is.null(input$sf.x)){
          
          sf.x <- 1
          
        } else {
          
          sf.x <- input$sf.x
          
        } 
        
        Vnans <- c()
        
        for (i in 1:length(input$Vyans)){
          
          Vnans[i] <- as.numeric(input[[paste0("nans", i)]])
          
        } 
        
        
        inputValues <- list(dtype = results$dtype(), 
            xans = as.numeric(input$xans), 
            Vyans = as.numeric(input$Vyans),
            Vnans = Vnans,
            sans = as.numeric(input$sans), 
            sd.se = as.numeric(input$sd.se),
            sf.x = sf.x,
            CES = input$CES,
            cont = results$cont(),
            conf.lev = input$conf.lev
        
        )
        
    
        inputValues
        
      })
  
  
 
  
  output$summaryAnalysis <- renderUI({
        
        input$submit
        
        validate(need(input$submit, ""))
        
        if(input$submit == 0){
          
          return(NULL)
          
        }
        
        isolate({
              
              lapply(1:length(results$shinyInput()$Vyans), function(iResponse) {
                    
                    currentShinyInput <- results$shinyInput()
                    currentShinyInput$yans <- currentShinyInput$Vyans[iResponse]
                    currentShinyInput$Vyans <- NULL
                    currentShinyInput$nans <- currentShinyInput$Vnans[iResponse]
                    currentShinyInput$Vnans <- NULL
                    
                    fittedModels <- fitModels(data = results$dataModified(), 
                        shinyInput = currentShinyInput,
                        continuousResponse = results$cont())
                    
                    validate(need( !any(attr(fittedModels, "class") == "error"), 
                            paste("Error in calculation: \n", fittedModels$message))
                    )
                    
                    
                    tmpTable <- summaryModels(savedResults = fittedModels)
                    
                    colnames(tmpTable) <- c("Model", "Number of parameters", 
                        "Log-likelihood", "AIC", "BMD", "BMDL", "BMDU", "Converged",
                        "Accepted AIC")
                    
                    
                    fluidRow(
                        column(6,
                            
                            h4(p("Response: ", attr(tmpTable, "responseName"))),
                            renderTable(tmpTable),
                            p("Lowest BMDL: ", round(attr(tmpTable, "minBmdl")), 2),
                            p("Highest BMDU: ", round(attr(tmpTable, "maxBmdu")), 2),
                            p("Best model(s): ", attr(tmpTable, "bestModel"))
                        
                        ),
                        
                        column(6, 
                            
                            renderPlot(
                                
                                f.plot.all(fittedModels[[attr(tmpTable, "bestModelIndex")]], track = track)
                            
                            )
                        
                        )
                    
                    )
                    
                  })
              
            })
        
      })
  
  
  
  output$warnings <- renderUI({
        
        list(
            
            busyIndicator("In progress", wait = 0),
            actionButton(inputId = "submit", label = "Fit model", 
                styleclass = "primary", size = "mini"),
            h6("Output generated with Proast, version 61.3", align = "right")
        
        )
        
      })
  
  
}