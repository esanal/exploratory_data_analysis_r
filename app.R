library(shiny)
library(DT)
library(tidyverse)
library(Cairo)
library(ggplot2)
library(gridExtra)
library(seplyr)
library(shinythemes)
library(plotly)
library(mygene)

options(shiny.maxRequestSize=30*1024^2) 
# User Interface
ui <- fluidPage(
  # App title
  titlePanel(title = h1("Exploratory data analysis", align = "center"), windowTitle = "Proteome Explorer"),
  
  tabsetPanel(
    #Tab Panel 1: Main upload ----
    tabPanel("Data", fluid = TRUE, type = "tabs",
             
             #Sidebar layout with input and output definitions  
             sidebarLayout(
               # Sidebar panel for inputs  
               sidebarPanel(
                 
                 # Input: Select a file  
                 fileInput("uploaded_file", "Choose File",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # Horizontal line   
                 tags$hr(),
                 
                 # Input: Checkbox if file has header  
                 checkboxInput("header", "File has header", TRUE),
                 
                 # Input: Select separator  
                 radioButtons("sep", "Separator",
                              choices = c(Semicolon = ";",
                                          Comma = ",",
                                          Tab = "\t"),
                              selected = "\t"),
                 
                 # Input: Select number of rows to display  
                 radioButtons("disp", "Display",
                              choices = c(All = "all",
                                          Top = "head"),
                              selected = "all"),
                 # Horizontal line  
                 tags$hr(),
                 
                 # Select variables to display  
                 uiOutput("checkbox"),
                 uiOutput("filter_columns"),
                 uiOutput("filter_from_elements")
                 ),
               
               # Main panel for displaying outputs  
               mainPanel(id = "dataset",
                         DT::dataTableOutput("rendered_file"),
                         radioButtons(inputId = "filetype",
                                      label = "Select filetype:",
                                      choices = c("csv", "tsv"),
                                      selected = "csv"),
                         downloadButton(outputId = "download_data", label = "Download data"))
             )#End of sidebarlayout  
    ), 
    #End of Tab Panel 1: Main upload ----
    
    
    #Tab Panel 2: Statistics ----
    tabPanel('Statistics', fluid = TRUE,
             
             sidebarLayout(
               
               sidebarPanel(wellPanel(uiOutput("select_ratio_column"),
                                      uiOutput("select_intensity_column")),
                            sliderInput("MAD", 'Set cut-off %', min = 0, max = 100, value = 95),
                            wellPanel(
                                      selectInput("normalize_by", 'Normalize by:', c('Median','Mean')),
                                      div(style="display:inline-block",actionButton("normalize_button", "Normalize"),width=6,style="display:center-align")
                                      ),
                            wellPanel(
                                      numericInput("log_scale_y", 'Scale Intensity to log base:', 10),
                                      numericInput("log_scale_x", 'Scale Ratio to log base:', 2),
                                      actionButton("log_scale_button", "Scale")
                                      )
               ),
               
               mainPanel(h4("Summary of selected ratio column"),
                         verbatimTextOutput("summary_text_before_norm",placeholder = TRUE),
                         h4("Summary of selected ratio column after normalization"),
                         verbatimTextOutput("summary_text_after_norm", placeholder = TRUE),
                         h4("Test parameters"),
                         verbatimTextOutput("mOut", placeholder = TRUE),
                         DT::dataTableOutput("rendered_file_normalized")
               )
               
             )#End of sidebar layout
    ),
    #End of Tab Panel 2: Statistics ----
    
    #Tab Panel 3: Scatter plot ----
    tabPanel("Scatter-Plot", fluid = TRUE,
             # Sidebar layout with input and output definitions  
             sidebarLayout(
                           # Sidebar panel for inputs  
                           sidebarPanel(
                             uiOutput("plot_y"),
                             uiOutput("plot_x")
                           ),
               
               mainPanel(
                         downloadButton("scatterplot_download",
                                        label = "Download Scatterplot"),
                 
                         plotOutput(outputId = "scatterplot", 
                                    dblclick = "scatterplot_dblclick",
                                    brush = brushOpts(
                                      id = "scatterplot_brush",
                                      resetOnNew = TRUE
                                    )
                 ),
                 # Show data table
                 dataTableOutput(outputId = "brushed_points")
               )
             ) #End of sidebarLayout
    ),
    #End of Tab Panel 3: Scatter plot ----
    
    #Tab Panel 3.1: Scatter plot2 ----
    tabPanel("Scatter-Plot2", fluid = TRUE,
             # Sidebar layout with input and output definitions  
             sidebarLayout(
               # Sidebar panel for inputs  
               sidebarPanel(
                 uiOutput("plot_y1"),
                 uiOutput("plot_x1")
               ),
               
               mainPanel(
                 htmltools::div(style = "display:inline-block", plotlyOutput("xy"))
                 # Show data table
               )
             ) #End of sidebarLayout
    ),
    #End of Tab Panel 3: Scatter plot ----
    

    #Tab Panel 4: Histogram ----
    tabPanel("Histogram", fluid = TRUE,
             # Sidebar layout with input and output definitions  
             sidebarLayout(
                           # Sidebar panel for inputs  
                           sidebarPanel(
                                       uiOutput("plot_x_histogram"),
                                       numericInput('h_bin_size', "Bin size:", value = 50)
                           ),
               
               mainPanel(downloadButton("histogram_download", label = "Download Histogram"),
                         plotOutput(outputId = "histogram")
                        )
             ) #End of sidebarLayout
    ), 
    #End of Tab Panel 4: Histogram ----
    
    #Tab Panel 5: Annotation ----
    tabPanel("Annotation", fluid = TRUE,
             # Sidebar layout with input and output definitions  
             sidebarLayout(
               # Sidebar panel for inputs  
                            sidebarPanel(
                                          uiOutput("select_gene_name"),
                                          textInput('name_separator', "Name Separator:", value = ';'),
                                          textInput('new_gene_name_column', "New 'Gene Name' column:", value = 'GeneNameEdited'),
                                          actionButton("separate", "Generate Gene Names"),
                                          actionButton("add_annotation", "Add annotation")
                                        ),
               
               mainPanel(DT::dataTableOutput("GO_results_table"))
             )
    ) 
    #End of Tab Panel 5: Annotation ----
  )
)
# End of User Interface

# Define server logic to read selected file  
server <- function(input, output, session) {
  
  ####Reactive Values####
  # Store data in reactive value
  reactive_values <- reactiveValues(df_data = NULL, df_data_to_show = NULL, GO_results = NULL)
  #Define statistics values
  stat_values <- reactiveValues(normalize_count = 0)
  ####End of reactive values####
  
  
  #Oberve data upload
  observeEvent(input$uploaded_file, 
               {req(input$uploaded_file)
                uploaded_paths <- input$uploaded_file$datapath
                uploaded_data <- lapply(uploaded_paths, function(path) {
                  read.table(path,
                             header = input$header,
                             sep = input$sep,
                             stringsAsFactors = FALSE)
                })
                reactive_values$df_data <- dplyr::bind_rows(uploaded_data)
               }
              )
  
  #Observe filter button and filter data
  observeEvent(input$filter_button,
               {outdf <- NULL
               for (i in input$column_to_filter) {
                                                   if ( is.null(outdf) ){
                                                     outdf <- dplyr::filter(reactive_values$df_data_to_show, UQ(sym(i)) != input$filter_this)
                                                   }
                 
                                                   else{outdf <- dplyr::filter(outdf, UQ(sym(i)) != input$filter_this)
                                                   }
               }
               reactive_values$df_data_to_show <- outdf
               }
              )
  
  
  observeEvent(input$select_var,
               {reactive_values$df_data_to_show <- reactive_values$df_data[,input$select_var]})

  # Print data table  
  output$rendered_file <- DT::renderDataTable(reactive_values$df_data_to_show, filter = 'top', options = list(
    pageLength = 5, autoWidth = TRUE
  ))
  
  #Download widget
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("mydata")
    },
    content = function(file) { 
      write.csv(reactive_values$df_data_to_show, file, row.names = FALSE) 
    }
  )
  
  # Dynamically generate cols to be selected input when data is uploaded  
  output$checkbox <- renderUI({
                              selectInput(inputId = "select_var", 
                                          label = "Select variables", 
                                          choices = names(reactive_values$df_data),
                                          multiple = TRUE)
  })
  
  #define cols to be filtered
  output$filter_columns <- renderUI({
                                    tagList(
                                      selectInput(inputId = "column_to_filter", 
                                                  label = "Select columns to filter", 
                                                  choices = input$select_var,
                                                  multiple = TRUE),
                                      textInput("filter_this", "Remove rows containing:", "+"),
                                      actionButton(inputId = "filter_button",
                                                   label = "Filter")
                                            )
                                    })
  
  #Define ratio column
  output$select_ratio_column <- renderUI({
                                          req(input$select_var)
                                          selectInput(inputId = "selected_ratio_column", 
                                                      label = "Select ratio column", 
                                                      choices = input$select_var,
                                                      multiple = FALSE)
                                          })
  
  #Define intensity column
  output$select_intensity_column <- renderUI({
                                              req(input$select_var)
                                              selectInput(inputId = "selected_intensity_column", 
                                                          label = "Select intensity column", 
                                                          choices = input$select_var,
                                                          multiple = FALSE)
                                              })
  
  ####Calculate normalized values####
  #Generate summary for ratio column, normalize and MAD
  output$summary_text_before_norm <- renderText({
                                                  req(input$selected_ratio_column)
                                                  summary(data.frame(stat_values$selected_ratios))
                                                 })
  
  #Generate summary for after normalization
  output$summary_text_after_norm <- renderText({
                                                req(input$normalize_button)
                                                summary(data.frame(stat_values$normalized_ratio))
                                                })
  

  
  #Observe MAD quantile selection
  observeEvent(input$MAD, 
               {stat_values$mad_cutoff <- qnorm(input$MAD*0.01)
               }
               )
  
  #Observe changes and re-eavluate statistics
  observeEvent(input$selected_ratio_column,
               {stat_values$selected_ratios <- as.numeric(unlist(reactive_values$df_data_to_show %>% dplyr::select(input$selected_ratio_column)))
                stat_values$m <- median(stat_values$selected_ratios)
               }
               )
  
  #Add intensity column to reactive values
  observeEvent(input$selected_intensity_column,
               {stat_values$selected_intensity <- as.numeric(unlist(reactive_values$df_data_to_show %>% dplyr::select(input$selected_intensity_column)))
               stat_values$m_intensity <- median(stat_values$selected_intensity)
               }
               )
  
  #Normalize data
  observeEvent(input$normalize_button,
               {req(stat_values$selected_ratios)
                normalization_center <- if (input$normalize_by == "Mean") {
                  mean(stat_values$selected_ratios, na.rm = TRUE)
                } else {
                  median(stat_values$selected_ratios, na.rm = TRUE)
                }
                stat_values$normalized_ratio <- stat_values$selected_ratios / normalization_center
                stat_values$normalized_ratio_median <- median(stat_values$normalized_ratio, na.rm = TRUE)
                stat_values$mad <- mad(stat_values$normalized_ratio, na.rm = TRUE)
                stat_values$outliers_bool <- ifelse(
                  ((abs(stat_values$normalized_ratio - stat_values$normalized_ratio_median) / stat_values$mad) >
                     stat_values$mad_cutoff),
                  "Outlier",
                  "NS"
                )
                
                reactive_values$df_data_to_show$Normalized_Ratio <- stat_values$normalized_ratio
                reactive_values$df_data_to_show$Outlier <- stat_values$outliers_bool
                stat_values$normalize_count <- stat_values$normalize_count + 1
               }
               )
  #Log Scale
  observeEvent(input$log_scale_button,
               {req(stat_values$normalized_ratio, stat_values$selected_intensity)
                stat_values$normalized_ratio_log <- log(stat_values$normalized_ratio, base= input$log_scale_x)
                stat_values$intensity_log <- log(stat_values$selected_intensity, base= input$log_scale_y)
                reactive_values$df_data_to_show$log_Normalized_Ratio <- stat_values$normalized_ratio_log
                reactive_values$df_data_to_show$log_Intensity <- stat_values$intensity_log
               }
               )
  
  #Output statistics
  output$mOut <- renderText({
    req(input$selected_ratio_column)
    paste("The median for your selection is: ", stat_values$m,
          "\nThe median for your selection after normalization is: ", stat_values$normalized_ratio_median,
          "\nThe MAD for your selection after normalization is: ", stat_values$mad,
          "\nCut off MAD for the quantile is: ", stat_values$mad_cutoff)
                             })
  
  #Plot options
  output$plot_y <- renderUI({
    selectInput(inputId = "plot_y_on", 
                label = "Y-axis:", 
                choices = names(reactive_values$df_data_to_show),
                multiple = FALSE)
                             })
  
  output$plot_x <- renderUI({
    selectInput(inputId = "plot_x_on", 
                label = "X-axis:", 
                choices = names(reactive_values$df_data_to_show),
                multiple = FALSE)
                             })

  output$plot_y1 <- renderUI({
    selectInput(inputId = "plot_y_on_1",
                label = "Y-axis:",
                choices = names(reactive_values$df_data_to_show),
                multiple = FALSE)
  })

  output$plot_x1 <- renderUI({
    selectInput(inputId = "plot_x_on_1",
                label = "X-axis:",
                choices = names(reactive_values$df_data_to_show),
                multiple = FALSE)
  })
  
  
  #plot scatter zoomable#
  ranges <- reactiveValues(x = NULL, y = NULL, regular_color = 'darkgray', outlier_color = 'darkred' )
  
  scatter_plot_f <- function(){
    ggplot(data = reactive_values$df_data_to_show, 
           aes_string(x = input$plot_x_on,
                      y = input$plot_y_on)
    ) +
      geom_point(colour = ifelse(reactive_values$df_data_to_show$Outlier=='Outlier', ranges$outlier_color, ranges$regular_color), alpha = 0.9) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y,
                      expand = TRUE) +
      theme_bw()
  }
  
  output$scatterplot <- renderPlot({
    scatter_plot_f()
  })
  
  observeEvent(input$scatterplot_dblclick, 
               {
    brush <- input$scatterplot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } 
    
    else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  }
               )
  
  ####Download plot as pdf####
  ### function to do the pdf export
  output$scatterplot_download <- downloadHandler("test.pdf", function(theFile) {
    
    makePdf <- function(filename){
      # I use Cairo instead of the the basic pdf function. It is much more complex and it supports special characters like "u"  
      Cairo(type = 'pdf', file = filename, width = 21, height = 29.7, units='cm', bg='transparent')
      p <- scatter_plot_f()
      print(p)
      dev.off()
    }
    
    makePdf(theFile)
  })
  #brushed data
  output$brushed_points <- DT::renderDataTable({
    brushedPoints(reactive_values$df_data_to_show, input$scatterplot_brush, xvar = input$plot_x_on, yvar = input$plot_y_on) 
  })
  
  ####Scatterplot 2####
  m <- list(color = toRGB("black"))
  m2 <- list(color = toRGB("black", 0.2))
  
  #output$xy <- renderPlotly({
  #  reactive_values$df_data_to_show %>% 
  #    plot_ly(x = ~input$plot_x_on, y = ~input$plot_y_on)
  #})
  
  output$xy <- renderPlotly({
    req(input$plot_x_on_1, input$plot_y_on_1)
    # use the key aesthetic/argument to help uniquely identify selected observations
    key <- row.names(reactive_values$df_data_to_show)

      plot_ly(reactive_values$df_data_to_show, x = ~get(input$plot_x_on_1), y = ~get(input$plot_y_on_1), key = ~key, type = "scatter") %>%
        layout(dragmode = "select")
  })
  
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
  })
  
  
  ####Histogram####
  ####Histogram plot####
  
  #x axis of histogram
  output$plot_x_histogram <- renderUI({
    selectInput(inputId = "plot_x_on_histogram", 
                label = "X-axis:", 
                choices = names(reactive_values$df_data_to_show),
                multiple = FALSE)
  })
  
  #generate histogram
  
  histogram_f <- function(){
    ggplot(data = reactive_values$df_data_to_show, aes_string(x = input$plot_x_on_histogram)) +
      geom_histogram(color="black", fill="white", bins = input$h_bin_size) +
      theme_bw()
  }
  
  output$histogram <- renderPlot({
    histogram_f()
  })
  
  ####Download plot as pdf####
  ### function to do the pdf export
  output$histogram_download <- downloadHandler("test.pdf", function(theFile) {
    makePdf <- function(filename){
      # I use Cairo instead of the the basic pdf function. It is much more complex and it supports special characters like "u"  
      Cairo(type = 'pdf', file = filename, width = 16, height = 16, units='cm', bg='transparent')
      
      p <- histogram_f()
      
      print(p)
      
      dev.off()
    }
    makePdf(theFile)}
                                                )
  #Gene annotation#
  #Generate gene names
  output$select_gene_name <- renderUI({
    selectInput(inputId = "select_gene_name_col", 
                label = "Gene name column:", 
                choices = names(reactive_values$df_data_to_show),
                multiple = FALSE)
  })
  #Gene name button observer
  observeEvent(input$separate,
               {
               gene_names <- reactive_values$df_data_to_show %>% dplyr::select(input$select_gene_name_col)
               print(typeof(gene_names))
               #reactive_values$df_data_to_show %>% tidyr::separate(select_gene_name_col, 
                #                      c("Gene name"), extra='drop')
               #reactive_values$df_data_to_show <- outdf
               gene_names <- unlist(gene_names, use.names = FALSE)
               gnl_l <- length(gene_names)
               gnl <- character(gnl_l)
               k <- 0
               #print(gene_names)
               for (i in gene_names){
                 k <- k + 1
                 #print(i)
                 a <- strsplit(i, ';')
                 a<-unlist(a, use.names=FALSE)
                 #print(a[1])
                 gnl[k] <- ""
                 for (l in a) {
                   if (l == toupper(l)){
                     gnl[k] <- l
                     next
                     next
                   }
                  }
               }
               reactive_values$df_data_to_show$GeneNames <- gnl
               }
  )
  #Gene name button observer
  observeEvent(input$add_annotation,
               {
                 gnl<-reactive_values$df_data_to_show$GeneNames
                 res <- queryMany(gnl, scopes='symbol', fields=c('go'), species='human', returnall=TRUE, return.as = "records")
                 annotations <- vector(mode="list", length=length(res$response))
                 for (i in 1:length(res$response)){
                   gene_name <- res$response[[i]]$query
                   bp_cur <- c()
                   mf_cur <- c()
                   cc_cur <- c()
                   if (is.list(res$response[[i]]$go$BP[[1]])){
                    for (j in 1:length(res$response[[i]]$go$BP)){
                      bp_cur <- c(bp_cur,res$response[[i]][["go"]][["BP"]][[j]][["term"]])
                   }
                   }
                   
                   if (is.list(res$response[[i]]$go$MF[[1]])){
                   for (k in 1:length(res$response[[i]]$go$MF)){
                     mf_cur <- c(mf_cur,res$response[[i]][["go"]][["MF"]][[k]][["term"]])
                   }
                   }
                   if (is.list(res$response[[i]]$go$CC[[1]])){
                   for (l in 1:length(res$response[[i]]$go$CC)){
                     cc_cur <- c(cc_cur,res$response[[i]][["go"]][["CC"]][[l]][["term"]])
                   }
                   }
                   annotations[[i]] <- list(cc_cur, bp_cur, mf_cur)
                   names(annotations[[i]]) <- c('GO_Cellular_Compartment',
                                                'GO_Biological_Process',
                                                'GO_Molecular_Function'
                                                )
                   names(annotations)[i] <- gene_name
                 }
                 # Print GO table  
                 #output$GO_results_table <- DT::renderDataTable(reactive_values$GO_results_BP, filter = 'top',
                  #                                              options = list(pageLength = 5, autoWidth = TRUE))
                 gnl_l <- length(gnl)
                 cc_add <- c()
                 mf_add <- c()
                 bp_add <- c()
                 for (m in 1:gnl_l){
                   if (gnl[m] %in% names(annotations)){
                     mf <- paste0(annotations[[gnl[m]]]$GO_Molecular_Function, collapse = ', ')
                     bp <- paste0(annotations[[gnl[m]]]$GO_Biological_Process, collapse = ', ')
                     cc <- paste0(annotations[[gnl[m]]]$GO_Cellular_Compartment, collapse = ', ')
                   }
            
                   else {
                     mf <- ""
                     bp <- ""
                     cc <- ""
                   }
                   cc_add <- c(cc_add,cc)
                   bp_add <- c(bp_add,bp)
                   mf_add <- c(mf_add,mf)
                   
                 }
                 reactive_values$df_data_to_show$GO_MF <- mf_add
                 reactive_values$df_data_to_show$GO_BP <- bp_add
                 reactive_values$df_data_to_show$GO_CC <- cc_add
               }
  )
  
  
}

# Create Shiny app  
shinyApp(ui, server)
