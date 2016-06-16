library(shinysky)

fluidPage(
    
    helpText(
        h5(a(href="https://efsa-models.openanalytics.eu/projects/bmd/wiki", 
                target="_blank", "About"), align = "right"),
        h5(a(href="https://efsa-models.openanalytics.eu/projects/bmd/issues/new/", 
                target="_blank", "Report new issue"), align = "right")
    ),
    
    titlePanel(title = div(img(src = "EFSA_logo.JPG", 
                float = "top", height = "60px", hspace = "50px"),
            "Benchmark Dose Modelling"), 
        windowTitle = "Benchmark Dose Modelling"),
    
    tags$br(),
    
    fluidRow(
        
        column(6, 
            
            wellPanel(
                
                tabsetPanel(id = "settings",
                    
                    tabPanel("Load Data",
                        
                        fluidRow(
                            column(6, selectInput("dataType", 'Type of data',
                                    choices = c("Proast data format" = "proastData", 
                                        "Raw data" = "rawData"))),
                            column(6, fileInput('dataLoaded', '',
                                    accept=c('text/csv', 
                                        'text/comma-separated-values,text/plain', 
                                        '.csv')))
                        ),
                        
                        conditionalPanel("input.dataType == 'proastData'",
                            
                            actionLink(inputId = "helpProast", label = "Proast data format",
                                icon = icon("info-circle")),
                            
                            conditionalPanel("input.helpProast % 2 == 1",
                                p(strong("Line 1:"), "one-word title, used in plots (no minus signs or spaces allowed)",
                                    br(), strong("Line 2:"), "the number of columns of the data matrix",
                                    br(), strong("Line 3:"), "a code for the data type in that particular column, choose between",
                                    tags$ul(
                                        tags$li("0: Non-response"),
                                        tags$li("1: Continuous response"),
                                        tags$li("2: Binary response"),
                                        tags$li("3: Ordinal response"),
                                        tags$li("4: Quantal response"),
                                        tags$li("5: Nested continuous response"),
                                        tags$li("6: Nested quantal response"),
                                        tags$li("10: Mean (continuous) response")
                                    ),
                                    strong("Line 4:"), "one-word title for each column (no minus signs or spaces allowed)",
                                    br(), strong("Line 5 onwards:"), "the data matrix")
                            ),
                            
                            conditionalPanel("input.helpProast % 2 == 0",
                                tags$br()
                            )
                        
                        ),
                        
                        radioButtons('sep', 'Separator for data',
                            c("Comma (.csv file)"=',',
                                "Semicolon (.csv file)"=';',
                                "Tab (.txt file)"='\t',
                                "Space (.txt file)" = ' ')
                        ),
                        
                        conditionalPanel("input.dataType == 'rawData'",
                            
                            textInput('dataName', label = 'Name for data', value = ""),
                            
                            uiOutput('covariatesType')
                        
                        )
                    ),
                    
                    tabPanel("General Parameters",
                        
                        # Main questions
                        selectInput("dtype", 
                            label = "What type of response data do you want to consider?",
                            choices = c("continuous" = 1, "quantal" = 4, "continuous, summary data" = 10)),
#                        choices = c("continuous" = 1, "binary" = 2, "ordinal" = 3, "quantal" = 4, 
#                                        "continuous, clustered" = 5, "quantal, clustered" = 6, 
#                                        "continuous, summary data" = 10, 
#                                        "continuous, summary data, clustered" = 15, 
#                                        "quantal, CxT" = 84)),
                        
                        
                        uiOutput("parameterQuestions"),
                        uiOutput("parameterQuestions2"),
                        
                        numericInput("CES", label = "Value for CES (positive)", value = 0.05),
                        
                        actionLink(inputId = "helpCES", label = "CES",
                            icon = icon("info-circle")),
                        
                        conditionalPanel("input.helpCES % 2 == 1",
                            p(em("Critical Effect Size, i.e. the Benchmark Response (BMR) 
                                        defined as a percent change in average response compared to the response in the controls"))
                        ),
                        
                        conditionalPanel("input.helpCES % 2 == 0",
                            tags$br()
                        )
                    
                    ),
                    
                    tabPanel("Other Settings",
                        
                        numericInput("sf.x", label = "Q10: Give scaling factor for the independent variable", value = 1),
                        
                        numericInput("conf.lev", "Confidence level for the confidence intervals", 
                            value = 0.9),
                        
                        uiOutput("transformation"),
                        
                        
                        
# TODO parameter constraints as input
                        
                        
# TODO factors and scale.dum as
                        
# TODO shinyInput$constr.ans for quantal response
#constr.ans <- menu(c("no", "yes (not recommended!)"), 
#    title = "\nDo you want to constrain the models to have finite slope at zero? ")
                        
                        hotable("parameterConstraints"),
                        
                        uiOutput("CES")
                    
                    )
                )
            
            
            )
        
        ),
        
        column(6,
            
            conditionalPanel("input.settings == 'Load Data'",
                
                DT::dataTableOutput('dataLoaded')
            
            ),
            
            conditionalPanel("input.settings == 'General Parameters'",
                
                p("TODO: Plot data")
            
            )
        
        )
    
    ),
    
    
    tags$hr(),
    
#    verbatimTextOutput("print"),
    
    actionButton("example", "Example", styleclass = "success"),
    tags$br(),
    tags$br(),
    uiOutput("warnings"),
    
    uiOutput("summaryAnalysis")

)