library(shiny)
library(DT)
library(tidyverse)
library(Cairo)
library(ggplot2)
library(gridExtra)
library(seplyr)
library(shinythemes)
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
                 
                 # Horizontal line  
                 tags$hr(),
                 
                 # Input: Select number of rows to display  
                 radioButtons("disp", "Display",
                              choices = c(All = "all",
                                          Top = "head"),
                              selected = "all"),
                 
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
               
               sidebarPanel(uiOutput("select_ratio_column"),
                            uiOutput("select_intensity_column"),
                            sliderInput("MAD", 'Set cut-off %', min = 0, max = 100, value = 95),
                            selectInput("normalize_by", 'Normalize by:', c('Median','Mean')),
                            actionButton("outlier_button", "Analyze")
               ),
               
               mainPanel(verbatimTextOutput("summary_text_before_norm",placeholder = TRUE),
                         verbatimTextOutput("mOut"),
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
                             uiOutput("plot_x"),
                             textInput("log_scale_y", 'Scale Y to log base:', 10),
                             textInput("log_scale_x", 'Scale X to log base:', 2)
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
    
    #Tab Panel 4: Histogram ----
    tabPanel("Histogram", fluid = TRUE,
             # Sidebar layout with input and output definitions  
             sidebarLayout(
                           # Sidebar panel for inputs  
                           sidebarPanel(
                                       uiOutput("plot_x_histogram"),
                                       textInput('h_bin_size', "Bin size:", value = 50)
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
               sidebarPanel(),
               
               mainPanel()
             )
    ) 
    #End of Tab Panel 5: Annotation ----
    
  )
)
# End of User Interface

# Define server logic to read selected file  
server <- function(input, output, session) {
  
  # Store data in reactive value
  reactive_values <- reactiveValues(df_data = NULL)
  
  #Oberve data upload
  observeEvent(input$uploaded_file, 
               {reactive_values$df_data <-  read.table(input$uploaded_file$datapath,
                                                       header = input$header,
                                                       sep = input$sep,
                                                       stringsAsFactors=F)
               }
              )
  
  observeEvent(input$Go,
               {temp <- values$df_data[-input$Delete, ]
    values$df_data <- temp
               }
              )
  
  observeEvent(input$soil,
               {})
  
  #Selecting the columns from uploaded file
  df_sel <- reactive({
    req(input$select_var)
    df_sel <- df()[,input$select_var]
    print(length(input$select_var))
    return(df_sel)
  })
  
  # Filtering based on selected columns to filter from selected columns
  df_f <- reactive({
    req(input$column_to_filter)
    outdf <- NULL
    for (i in input$column_to_filter) {
      
      if ( is.null(outdf) ){
        outdf <- dplyr::filter(df_sel(), UQ(sym(i)) != input$filter_this)
      }
      
      else{outdf <- dplyr::filter(outdf, UQ(sym(i)) != input$filter_this)
      }
      
    }
    return (outdf)
  })
  
  # Normalized and outlier found dataframe
  df_f_o <- reactive({
    req(input$column_to_filter)
    outdf <- NULL
    for (i in input$column_to_filter) {
      
      if ( is.null(outdf) ){
        outdf <- dplyr::filter(df_sel(), UQ(sym(i)) != input$filter_this)
      }
      
      else{outdf <- dplyr::filter(outdf, UQ(sym(i)) != input$filter_this)
      }
      
    }
    return (cbind(outdf, Significance = stat_values$outliers_bool))
  })
  
  # Print data table  
  output$rendered_file <- DT::renderDataTable(reactive_values$df_data)
  
  
  #Download widget
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("mydata")
    },
    content = function(file) { 
      write.csv(df_sel(), file, row.names = FALSE) 
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
      textInput("filter_this", "Remove rows containing:", "+")
    )
  })
  
  #Define ratio column
  output$select_ratio_column <- renderUI({
    selectInput(inputId = "selected_ratio_column", 
                label = "Select ratio column", 
                choices = input$select_var,
                multiple = FALSE)
  })
  
  #Define intensity column
  output$select_intensity_column <- renderUI({
    selectInput(inputId = "selected_intensity_column", 
                label = "Select intensity column", 
                choices = input$select_var,
                multiple = FALSE)
  })
  
  ####Calculate normalized values####
  #Generate summary for ratio column, normalize and MAD
  output$summary_text_before_norm <- renderText({
    c(summary(df_sel() %>% select(input$selected_ratio_column)))
  })
  
  #Define statistics values
  stat_values <- reactiveValues()
  
  #Observe changes and re-eavluate statistics
  observe({
    req(input$selected_ratio_column)
    stat_values$selected_ratios <- as.numeric(unlist(df_sel() %>% select(input$selected_ratio_column)))
    stat_values$m <- median(stat_values$selected_ratios)
    stat_values$mad_cutoff <- qnorm(input$MAD*0.01)
    stat_values$mad <- mad(stat_values$selected_ratios)
    #stat_values$outliers_bool <- which(((stat_values$selected_ratios - stat_values$m) / stat_values$mad) > stat_values$mad_cutoff)
    stat_values$outliers_bool <- ifelse( ( ( abs(stat_values$selected_ratios - stat_values$m) / stat_values$mad ) > stat_values$mad_cutoff), "Outlier", "NS")
    print(stat_values$outliers_bool)
    #df_sel() <- cbind(df_sel(), stat_values$outliers_bool)
  })
  
  #Output statistics
  output$mOut <- renderText({
    req(input$selected_ratio_column)
    paste("The median for your selection is: ", stat_values$m,
          "\nThe MAD for your selection is: ", stat_values$mad,
          "\nCut off MAD for the quantile is: ", stat_values$mad_cutoff)
  })
  
  #Plot options
  output$plot_y <- renderUI({
    selectInput(inputId = "plot_y_on", 
                label = "Y-axis:", 
                choices = input$select_var,
                multiple = FALSE)
  })
  
  output$plot_x <- renderUI({
    selectInput(inputId = "plot_x_on", 
                label = "X-axis:", 
                choices = input$select_var,
                multiple = FALSE)
  })
  
  
  #plot scatter zoomable#
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  scatter_plot_f <- function(){
    ggplot(data = df_sel(), 
           aes_string(x = input$plot_x_on,
                      y = input$plot_y_on),
           color=factor(season)
    ) +
      geom_point(colour = "darkgreen", alpha = 0.9) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y,
                      expand = FALSE) +
      theme_bw()
  }
  
  output$scatterplot <- renderPlot({
    scatter_plot_f()
  })
  
  observeEvent(input$scatterplot_dblclick, {
    brush <- input$scatterplot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } 
    
    else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
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
    brushedPoints(df(), input$scatterplot_brush, xvar = input$plot_x_on, yvar = input$plot_y_on) 
  })
  
  ####Histogram####
  ####Histogram plot####
  
  #x axis of histogram
  output$plot_x_histogram <- renderUI({
    selectInput(inputId = "plot_x_on_histogram", 
                label = "X-axis:", 
                choices = input$select_var,
                multiple = FALSE)
  })
  
  #generate histogram
  
  histogram_f <- function(){
    ggplot(data = df(), aes_string(x = input$plot_x_on_histogram)) +
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
    
    makePdf(theFile)
  })
  #       #
  
  
  
  
}

# Create Shiny app  
shinyApp(ui, server)