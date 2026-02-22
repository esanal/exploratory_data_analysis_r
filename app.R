library(shiny)
library(DT)
library(tidyverse)
library(Cairo)
library(ggplot2)
library(gridExtra)
#library(seplyr)
library(shinythemes)
library(plotly)
library(mygene)

options(shiny.maxRequestSize=200*1024^2) 

collapse_preview <- function(x, max_items = 8) {
  if (length(x) == 0) {
    return("none")
  }
  x <- as.character(x)
  if (length(x) > max_items) {
    paste0(paste(x[seq_len(max_items)], collapse = ", "), ", ...")
  } else {
    paste(x, collapse = ", ")
  }
}

read_uploaded_table <- function(path, sep, has_header) {
  df <- readr::read_delim(
    file = path,
    delim = sep,
    col_names = has_header,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN", "NULL"),
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  names(df) <- make.names(names(df), unique = TRUE)
  df
}

validate_uploaded_schemas <- function(uploaded_data) {
  if (length(uploaded_data) <= 1) {
    return(invisible(TRUE))
  }
  reference_cols <- names(uploaded_data[[1]])
  reference_file <- names(uploaded_data)[1]
  for (idx in seq_along(uploaded_data)[-1]) {
    current_cols <- names(uploaded_data[[idx]])
    if (!identical(reference_cols, current_cols)) {
      missing_cols <- setdiff(reference_cols, current_cols)
      extra_cols <- setdiff(current_cols, reference_cols)
      stop(
        paste0(
          "Schema mismatch between '", reference_file, "' and '", names(uploaded_data)[idx], "'. ",
          "Missing columns: ", collapse_preview(missing_cols), ". ",
          "Unexpected columns: ", collapse_preview(extra_cols), "."
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

drop_empty_columns <- function(df, ignore_cols = character()) {
  is_empty <- function(x) {
    if (is.character(x)) {
      all(is.na(x) | stringr::str_trim(x) == "")
    } else {
      all(is.na(x))
    }
  }
  keep_cols <- vapply(
    names(df),
    function(col_name) {
      col_name %in% ignore_cols || !is_empty(df[[col_name]])
    },
    logical(1)
  )
  df[, keep_cols, drop = FALSE]
}

load_uploaded_files <- function(uploaded_file, sep, has_header) {
  parsed_files <- lapply(seq_len(nrow(uploaded_file)), function(idx) {
    read_uploaded_table(uploaded_file$datapath[idx], sep = sep, has_header = has_header)
  })
  names(parsed_files) <- uploaded_file$name

  validate_uploaded_schemas(parsed_files)

  parsed_with_source <- lapply(seq_along(parsed_files), function(idx) {
    parsed_files[[idx]] %>%
      dplyr::mutate(
        source_file = uploaded_file$name[idx],
        source_row_id = dplyr::row_number()
      )
  })

  combined <- dplyr::bind_rows(parsed_with_source)
  combined_clean <- drop_empty_columns(combined, ignore_cols = c("source_file", "source_row_id"))
  dropped_empty_cols <- setdiff(names(combined), names(combined_clean))
  combined_typed <- suppressWarnings(
    readr::type_convert(
      combined_clean,
      na = c("", "NA", "NaN", "NULL"),
      guess_integer = TRUE
    )
  )
  combined_typed <- combined_typed %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::relocate(row_id, .before = 1)

  list(data = combined_typed, dropped_empty_cols = dropped_empty_cols)
}

measurement_columns_from_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(character())
  }
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  numeric_cols[stringr::str_detect(numeric_cols, "^(Ratio|Intensity|Fraction|Normalized|log)")]
}

build_measurement_long <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble::tibble())
  }
  id_cols <- intersect(c("row_id", "source_file", "source_row_id"), names(df))
  value_cols <- measurement_columns_from_df(df)
  if (length(value_cols) == 0) {
    return(tibble::tibble())
  }

  df %>%
    dplyr::select(dplyr::all_of(id_cols), dplyr::all_of(value_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(value_cols),
      names_to = "metric_column",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric_family = dplyr::case_when(
        stringr::str_starts(metric_column, "Ratio") ~ "Ratio",
        stringr::str_starts(metric_column, "Intensity") ~ "Intensity",
        stringr::str_starts(metric_column, "Fraction") ~ "Fraction",
        stringr::str_starts(metric_column, "Normalized") ~ "Normalized",
        stringr::str_starts(metric_column, "log") ~ "Log",
        TRUE ~ "Other"
      ),
      sample_id = stringr::str_extract(metric_column, "(IEF\\d+_\\d+|Fraction\\.\\d+)$"),
      metric_name = stringr::str_remove(metric_column, "\\.(IEF\\d+_\\d+|Fraction\\.\\d+)$"),
      sample_id = dplyr::if_else(is.na(sample_id), "global", sample_id),
      sample_id = stringr::str_replace(sample_id, "^Fraction\\.", "Fraction_")
    )
}

compose_display_data <- function(view_df, row_features_df) {
  if (is.null(view_df)) {
    return(NULL)
  }
  if (is.null(row_features_df) || !("row_id" %in% names(view_df)) || !("row_id" %in% names(row_features_df))) {
    return(view_df)
  }
  derived_cols <- setdiff(names(row_features_df), names(view_df))
  if (length(derived_cols) == 0) {
    return(view_df)
  }
  join_df <- row_features_df[, c("row_id", derived_cols), drop = FALSE]
  dplyr::left_join(view_df, join_df, by = "row_id")
}

upsert_row_features <- function(base_df, updates_df) {
  if (is.null(base_df) || nrow(base_df) == 0 || is.null(updates_df) || nrow(updates_df) == 0) {
    return(base_df)
  }
  if (!"row_id" %in% names(base_df) || !"row_id" %in% names(updates_df)) {
    return(base_df)
  }

  target_idx <- match(updates_df$row_id, base_df$row_id)
  valid_rows <- !is.na(target_idx)
  if (!any(valid_rows)) {
    return(base_df)
  }

  target_idx <- target_idx[valid_rows]
  updates_df <- updates_df[valid_rows, , drop = FALSE]

  update_cols <- setdiff(names(updates_df), "row_id")
  for (col_name in update_cols) {
    if (!col_name %in% names(base_df)) {
      base_df[[col_name]] <- NA
    }
    base_df[[col_name]][target_idx] <- updates_df[[col_name]]
  }
  base_df
}

# User Interface
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      table.dataTable.nowrap th,
      table.dataTable.nowrap td {
        white-space: nowrap;
      }
      table.dataTable td.dt-truncate,
      table.dataTable th.dt-truncate {
        max-width: 180px;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      .annotation-progress-track {
        width: 100%;
        height: 12px;
        background: #e9ecef;
        border-radius: 999px;
        overflow: hidden;
        margin-top: 8px;
      }
      .annotation-progress-fill {
        height: 100%;
        background: #2c7fb8;
        width: 100%;
        transform: scaleX(0);
        transform-origin: left center;
        transition: transform 0.15s linear;
      }
      .annotation-status-text {
        margin-top: 6px;
      }
      .annotation-progress-label {
        margin-top: 4px;
        font-size: 12px;
        color: #4b5563;
      }
      div.dataTables_wrapper {
        width: 100%;
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('annotation-progress-update', function(message) {
        var fill = document.getElementById('annotation-progress-fill');
        var label = document.getElementById('annotation-progress-label');
        if (fill && typeof message.progress === 'number') {
          var pct = Math.max(0, Math.min(100, message.progress));
          fill.style.transform = 'scaleX(' + (pct / 100) + ')';
          fill.style.width = pct + '%';
        }
        if (label && typeof message.progress === 'number') {
          label.textContent = Math.round(message.progress) + '%';
        }
      });
    "))
  ),
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
                                      "text/tab-separated-values",
                                      "text/comma-separated-values,text/plain",
                                      ".csv",
                                      ".tsv",
                                      ".txt")),
                 
                 # Horizontal line   
                 tags$hr(),
                 
                 # Input: Species selection used for gene-name parsing
                 radioButtons(
                   "selected_species",
                   "Species",
                   choices = c("Human", "Mouse"),
                   selected = "Human"
                 ),
                 
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
                                          actionButton("add_annotation", "Add annotation"),
                                          tags$div(
                                            class = "annotation-progress-track",
                                            tags$div(id = "annotation-progress-fill", class = "annotation-progress-fill", style = "transform:scaleX(0); width:0%;")
                                          ),
                                          div(id = "annotation-progress-label", class = "annotation-progress-label", "0%"),
                                          div(class = "annotation-status-text", textOutput("annotation_status")),
                                          div(class = "annotation-status-text", textOutput("annotation_detail"))
                                        ),
               
               mainPanel(
                 h4("Annotated data"),
                 DT::dataTableOutput("rendered_file_annotation")
               )
             )
    ) 
    #End of Tab Panel 5: Annotation ----
  )
)
# End of User Interface

# Define server logic to read selected file  
server <- function(input, output, session) {
  
  ####Reactive Values####
  # Store raw/view/derived data separately.
  data_store <- reactiveValues(
    raw = NULL,
    view = NULL,
    row_features = NULL,
    measurement_long_raw = NULL,
    measurement_long_view = NULL,
    GO_results = NULL
  )
  #Define statistics values
  stat_values <- reactiveValues(normalize_count = 0)
  annotation_status <- reactiveVal("Idle")
  annotation_detail <- reactiveVal("")
  annotation_progress <- reactiveVal(0)
  ####End of reactive values####

  initialize_data_store <- function(base_df) {
    data_store$raw <- base_df
    data_store$view <- base_df
    data_store$row_features <- tibble::tibble(row_id = base_df$row_id)
    data_store$measurement_long_raw <- build_measurement_long(base_df)
    data_store$measurement_long_view <- build_measurement_long(base_df)
  }

  refresh_view_state <- function(new_view) {
    data_store$view <- new_view
    if (is.null(data_store$row_features)) {
      return(invisible(NULL))
    }
    current_view <- compose_display_data(data_store$view, data_store$row_features)
    data_store$measurement_long_view <- build_measurement_long(current_view)
    invisible(NULL)
  }

  update_derived_rows <- function(updates_df) {
    if (is.null(data_store$row_features)) {
      return(invisible(NULL))
    }
    data_store$row_features <- upsert_row_features(data_store$row_features, updates_df)
    current_view <- compose_display_data(data_store$view, data_store$row_features)
    data_store$measurement_long_view <- build_measurement_long(current_view)
    invisible(NULL)
  }

  current_data <- reactive({
    req(data_store$view)
    compose_display_data(data_store$view, data_store$row_features)
  })
  
  
  #Oberve data upload
  observeEvent(input$uploaded_file, 
               {req(input$uploaded_file)
                upload_result <- tryCatch(
                  {
                    load_uploaded_files(
                      uploaded_file = input$uploaded_file,
                      sep = input$sep,
                      has_header = input$header
                    )
                  },
                  error = function(e) {
                    showNotification(paste("Upload failed:", e$message), type = "error", duration = NULL)
                    NULL
                  }
                )
                req(upload_result)
                initialize_data_store(upload_result$data)
                if (length(upload_result$dropped_empty_cols) > 0) {
                  showNotification(
                    paste0("Dropped empty columns: ", collapse_preview(upload_result$dropped_empty_cols)),
                    type = "message"
                  )
                }
               },
               ignoreInit = TRUE
              )
  
  #Observe filter button and filter data
  observeEvent(input$filter_button,
               {req(data_store$view)
                req(input$column_to_filter)
                if (is.null(input$filter_this) || input$filter_this == "") {
                  showNotification("Enter text to filter rows.", type = "warning")
                  return()
                }
                outdf <- data_store$view
                for (col_name in input$column_to_filter) {
                  if (!col_name %in% names(outdf)) {
                    next
                  }
                  col_values <- dplyr::coalesce(as.character(outdf[[col_name]]), "")
                  outdf <- outdf[
                    !stringr::str_detect(col_values, pattern = stringr::fixed(input$filter_this)),
                    ,
                    drop = FALSE
                  ]
                }
                refresh_view_state(outdf)
               }
              )
  
  
  observeEvent(input$select_var,
               {req(data_store$raw)
                required_cols <- intersect(c("row_id", "source_file", "source_row_id"), names(data_store$raw))
                selected_cols <- unique(c(required_cols, input$select_var))
                if (is.null(input$select_var) || length(input$select_var) == 0) {
                  refresh_view_state(data_store$raw)
                } else {
                  valid_cols <- intersect(selected_cols, names(data_store$raw))
                  refresh_view_state(data_store$raw[, valid_cols, drop = FALSE])
                }
               })

  data_for_display <- reactive({
    req(current_data())
    if (identical(input$disp, "head")) {
      return(utils::head(current_data(), 25))
    }
    current_data()
  })

  scrollable_dt_options <- function(page_length = 25, scroll_height = "60vh") {
    list(
      pageLength = page_length,
      lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All")),
      autoWidth = FALSE,
      scrollX = TRUE,
      scrollY = scroll_height,
      scrollCollapse = TRUE,
      deferRender = TRUE,
      scroller = TRUE,
      columnDefs = list(
        list(
          targets = "_all",
          className = "dt-truncate",
          createdCell = DT::JS(
            "function(td, cellData) {",
            "  if (cellData !== null && cellData !== undefined) {",
            "    td.setAttribute('title', cellData.toString());",
            "  }",
            "}"
          )
        )
      )
    )
  }

  build_scrollable_table <- function(df, page_length = 25, scroll_height = "60vh") {
    DT::datatable(
      df,
      filter = "top",
      rownames = FALSE,
      extensions = c("Scroller"),
      class = "display compact nowrap stripe hover",
      options = scrollable_dt_options(page_length = page_length, scroll_height = scroll_height)
    )
  }

  # Print data table  
  output$rendered_file <- DT::renderDataTable({
    build_scrollable_table(data_for_display(), page_length = 25, scroll_height = "60vh")
  })
  output$rendered_file_normalized <- DT::renderDataTable({
    build_scrollable_table(current_data(), page_length = 25, scroll_height = "60vh")
  })
  output$rendered_file_annotation <- DT::renderDataTable({
    build_scrollable_table(current_data(), page_length = 25, scroll_height = "60vh")
  })
  output$annotation_status <- renderText({
    annotation_status()
  })
  output$annotation_detail <- renderText({
    annotation_detail()
  })
  
  update_annotation_progress <- function(progress = NULL, status = NULL, detail = NULL) {
    if (!is.null(status)) {
      annotation_status(status)
    }
    if (!is.null(detail)) {
      annotation_detail(detail)
    }
    if (!is.null(progress)) {
      progress_value <- max(0, min(100, progress))
      annotation_progress(progress_value)
      session$sendCustomMessage("annotation-progress-update", list(progress = progress_value))
      if (exists("flushReact", where = asNamespace("shiny"), inherits = FALSE)) {
        getFromNamespace("flushReact", "shiny")()
      }
    }
  }
  
  #Download widget
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("mydata_", format(Sys.Date(), "%Y%m%d"), ".", input$filetype)
    },
    content = function(file) { 
      export_df <- current_data()
      if (identical(input$filetype, "tsv")) {
        write.table(export_df, file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
      } else {
        write.csv(export_df, file, row.names = FALSE, na = "")
      }
    }
  )
  
  preferred_selected_cols <- c(
    "Gene.names",
    "Ratio.H.L",
    "Intensity"
  )
  
  # Dynamically generate cols to be selected input when data is uploaded  
  output$checkbox <- renderUI({
                              req(data_store$raw)
                              default_selected_cols <- intersect(preferred_selected_cols, names(data_store$raw))
                              selectizeInput(
                                inputId = "select_var",
                                label = "Select variables",
                                choices = names(data_store$raw),
                                selected = default_selected_cols,
                                multiple = TRUE,
                                options = list(
                                  closeAfterSelect = FALSE,
                                  plugins = list("remove_button")
                                )
                              )
  })
  
  #define cols to be filtered
  output$filter_columns <- renderUI({
                                    tagList(
                                      selectInput(inputId = "column_to_filter", 
                                                  label = "Select columns to filter", 
                                                  choices = input$select_var,
                                                  multiple = TRUE),
                                      textInput("filter_this", "Remove rows containing text:", ""),
                                      actionButton(inputId = "filter_button",
                                                   label = "Filter")
                                            )
                                    })

  numeric_metric_choices <- reactive({
    req(current_data())
    df <- current_data()
    setdiff(
      names(df)[vapply(df, is.numeric, logical(1))],
      c("row_id", "source_row_id")
    )
  })

  ratio_metric_choices <- reactive({
    ml <- data_store$measurement_long_view
    if (is.null(ml) || nrow(ml) == 0 || !all(c("metric_family", "metric_column") %in% names(ml))) {
      return(numeric_metric_choices())
    }
    ratio_cols <- ml %>%
      dplyr::filter(metric_family == "Ratio") %>%
      dplyr::distinct(metric_column) %>%
      dplyr::pull(metric_column)
    if (length(ratio_cols) == 0) {
      return(numeric_metric_choices())
    }
    ratio_cols
  })

  intensity_metric_choices <- reactive({
    ml <- data_store$measurement_long_view
    if (is.null(ml) || nrow(ml) == 0 || !all(c("metric_family", "metric_column") %in% names(ml))) {
      return(numeric_metric_choices())
    }
    intensity_cols <- ml %>%
      dplyr::filter(metric_family == "Intensity") %>%
      dplyr::distinct(metric_column) %>%
      dplyr::pull(metric_column)
    if (length(intensity_cols) == 0) {
      return(numeric_metric_choices())
    }
    intensity_cols
  })
  
  #Define ratio column
  output$select_ratio_column <- renderUI({
                                          req(current_data())
                                          selectInput(inputId = "selected_ratio_column", 
                                                      label = "Select ratio column", 
                                                      choices = ratio_metric_choices(),
                                                      multiple = FALSE)
                                          })
  
  #Define intensity column
  output$select_intensity_column <- renderUI({
                                              req(current_data())
                                              selectInput(inputId = "selected_intensity_column", 
                                                          label = "Select intensity column", 
                                                          choices = intensity_metric_choices(),
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
               {req(current_data())
                req(input$selected_ratio_column %in% names(current_data()))
                stat_values$selected_ratios <- as.numeric(current_data()[[input$selected_ratio_column]])
                stat_values$m <- median(stat_values$selected_ratios)
               }
               )
  
  #Add intensity column to reactive values
  observeEvent(input$selected_intensity_column,
               {req(current_data())
                req(input$selected_intensity_column %in% names(current_data()))
                stat_values$selected_intensity <- as.numeric(current_data()[[input$selected_intensity_column]])
               stat_values$m_intensity <- median(stat_values$selected_intensity)
               }
               )
  
  #Normalize data
  observeEvent(input$normalize_button,
               {req(current_data())
                req(input$selected_ratio_column %in% names(current_data()))
                stat_values$selected_ratios <- as.numeric(current_data()[[input$selected_ratio_column]])
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
                updates <- tibble::tibble(
                  row_id = current_data()$row_id,
                  Normalized_Ratio = stat_values$normalized_ratio,
                  Outlier = stat_values$outliers_bool
                )
                update_derived_rows(updates)
                stat_values$normalize_count <- stat_values$normalize_count + 1
               }
               )
  #Log Scale
  observeEvent(input$log_scale_button,
               {req(current_data())
                req(input$selected_intensity_column %in% names(current_data()))
                req("Normalized_Ratio" %in% names(current_data()))
                stat_values$normalized_ratio <- as.numeric(current_data()[["Normalized_Ratio"]])
                stat_values$selected_intensity <- as.numeric(current_data()[[input$selected_intensity_column]])
                stat_values$normalized_ratio_log <- log(stat_values$normalized_ratio, base= input$log_scale_x)
                stat_values$intensity_log <- log(stat_values$selected_intensity, base= input$log_scale_y)
                updates <- tibble::tibble(
                  row_id = current_data()$row_id,
                  log_Normalized_Ratio = stat_values$normalized_ratio_log,
                  log_Intensity = stat_values$intensity_log
                )
                update_derived_rows(updates)
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
    req(current_data())
    selectInput(inputId = "plot_y_on", 
                label = "Y-axis:", 
                choices = ratio_metric_choices(),
                multiple = FALSE)
                             })
  
  output$plot_x <- renderUI({
    req(current_data())
    selectInput(inputId = "plot_x_on", 
                label = "X-axis:", 
                choices = intensity_metric_choices(),
                multiple = FALSE)
                             })

  output$plot_y1 <- renderUI({
    req(current_data())
    selectInput(inputId = "plot_y_on_1",
                label = "Y-axis:",
                choices = ratio_metric_choices(),
                multiple = FALSE)
  })

  output$plot_x1 <- renderUI({
    req(current_data())
    selectInput(inputId = "plot_x_on_1",
                label = "X-axis:",
                choices = intensity_metric_choices(),
                multiple = FALSE)
  })
  
  
  #plot scatter zoomable#
  ranges <- reactiveValues(x = NULL, y = NULL, regular_color = 'darkgray', outlier_color = 'darkred' )
  
  scatter_plot_f <- function(){
    req(current_data())
    req(input$plot_x_on %in% names(current_data()))
    req(input$plot_y_on %in% names(current_data()))
    plot_df <- current_data()
    point_colors <- if ("Outlier" %in% names(plot_df)) {
      ifelse(plot_df$Outlier == "Outlier", ranges$outlier_color, ranges$regular_color)
    } else {
      ranges$regular_color
    }
    ggplot(data = plot_df, 
           aes_string(x = input$plot_x_on,
                      y = input$plot_y_on)
    ) +
      geom_point(colour = point_colors, alpha = 0.9) +
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
    req(input$plot_x_on %in% names(current_data()))
    req(input$plot_y_on %in% names(current_data()))
    brushed_df <- brushedPoints(current_data(), input$scatterplot_brush, xvar = input$plot_x_on, yvar = input$plot_y_on)
    build_scrollable_table(brushed_df, page_length = 10, scroll_height = "35vh")
  })
  
  ####Scatterplot 2####
  m <- list(color = toRGB("black"))
  m2 <- list(color = toRGB("black", 0.2))
  
  output$xy <- renderPlotly({
    req(input$plot_x_on_1, input$plot_y_on_1)
    req(current_data())
    req(input$plot_x_on_1 %in% names(current_data()))
    req(input$plot_y_on_1 %in% names(current_data()))
    # use the key aesthetic/argument to help uniquely identify selected observations
    key <- current_data()$row_id

      plot_ly(current_data(), x = ~get(input$plot_x_on_1), y = ~get(input$plot_y_on_1), key = ~key, type = "scatter") %>%
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
    req(current_data())
    selectInput(inputId = "plot_x_on_histogram", 
                label = "X-axis:", 
                choices = names(current_data()),
                multiple = FALSE)
  })
  
  #generate histogram
  
  histogram_f <- function(){
    req(current_data())
    ggplot(data = current_data(), aes_string(x = input$plot_x_on_histogram)) +
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
    req(current_data())
    default_gene_name_col <- if ("Gene.names" %in% names(current_data())) "Gene.names" else names(current_data())[1]
    selectInput(inputId = "select_gene_name_col", 
                label = "Gene name column:", 
                choices = names(current_data()),
                selected = default_gene_name_col,
                multiple = FALSE)
  })
  #Gene name button observer
  observeEvent(input$separate,
               {
                 update_annotation_progress(
                   progress = 0,
                   status = "Generating gene names...",
                   detail = "Scanning selected column values..."
                 )
                 tryCatch(
                   {
                     req(current_data())
                     req(input$select_gene_name_col %in% names(current_data()))
                     req(input$selected_species)
                     is_human_species <- identical(input$selected_species, "Human")
                     gene_names <- as.character(current_data()[[input$select_gene_name_col]])
                     gnl_l <- length(gene_names)
                     gnl <- character(gnl_l)
                     k <- 0
                     split_sep <- ifelse(is.null(input$name_separator) || input$name_separator == "", ";", input$name_separator)
                     for (i in gene_names){
                       k <- k + 1
                       i <- dplyr::coalesce(i, "")
                       a <- strsplit(i, split_sep, fixed = TRUE)
                       a<-unlist(a, use.names=FALSE)
                       gnl[k] <- ""
                       for (l in a) {
                         l <- stringr::str_trim(l)
                         is_match <- if (is_human_species) {
                           l != "" && l == toupper(l)
                         } else {
                           stringr::str_detect(l, "^[A-Z][a-z]*$")
                         }
                         if (is_match){
                           gnl[k] <- l
                           break
                         }
                       }
                       update_annotation_progress(
                         progress = (k / max(1, gnl_l)) * 100,
                         detail = paste0("Generating gene names (", k, "/", gnl_l, ")...")
                       )
                     }
                     updates <- tibble::tibble(
                       row_id = current_data()$row_id,
                       GeneNames = gnl
                     )
                     update_derived_rows(updates)
                     update_annotation_progress(
                       progress = 100,
                       status = paste0("Gene names generated successfully at ", format(Sys.time(), "%H:%M:%S"), "."),
                       detail = "Completed."
                     )
                   },
                   error = function(e) {
                     update_annotation_progress(
                       progress = 0,
                       status = paste("Gene name generation failed:", e$message),
                       detail = "Stopped."
                     )
                     showNotification(paste("Gene name generation failed:", e$message), type = "error")
                   }
                 )
               }
  )
  #Gene name button observer
  observeEvent(input$add_annotation,
               {
                 update_annotation_progress(
                   progress = 0,
                   status = "Adding GO annotations...",
                   detail = "Preparing GeneNames list..."
                 )
                 tryCatch(
                   {
                     req(current_data())
                     req("GeneNames" %in% names(current_data()))
                     gnl <- current_data()$GeneNames
                     update_annotation_progress(
                       progress = 10,
                       detail = "Fetching GO annotations from mygene..."
                     )
                     res <- queryMany(gnl, scopes='symbol', fields=c('go'), species=input$selected_species, returnall=TRUE, return.as = "records")
                     response_len <- length(res$response)
                     update_annotation_progress(
                       progress = 30,
                       detail = paste0("Fetched ", response_len, " records. Parsing GO fields...")
                     )
                     annotations <- vector(mode="list", length=length(res$response))
                     for (i in seq_along(res$response)){
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
                       if (i %% max(1, floor(max(1, response_len) / 20)) == 0 || i == response_len) {
                         update_annotation_progress(
                           progress = 30 + (i / max(1, response_len)) * 40,
                           detail = paste0("Parsing GO response (", i, "/", response_len, ")...")
                         )
                       }
                     }
                     update_annotation_progress(
                       progress = 70,
                       detail = "Mapping parsed GO annotations back to table rows..."
                     )
                     gnl_l <- length(gnl)
                     cc_add <- c()
                     mf_add <- c()
                     bp_add <- c()
                     for (m in 1:gnl_l){
                       if (gnl[m] %in% names(annotations)){
                         mf <- paste0(annotations[[gnl[m]]]$GO_Molecular_Function, collapse = ', ')
                         bp <- paste0(annotations[[gnl[m]]]$GO_Biological_Process, collapse = ', ')
                         cc <- paste0(annotations[[gnl[m]]]$GO_Cellular_Compartment, collapse = ', ')
                       } else {
                         mf <- ""
                         bp <- ""
                         cc <- ""
                       }
                       cc_add <- c(cc_add,cc)
                       bp_add <- c(bp_add,bp)
                       mf_add <- c(mf_add,mf)
                       if (m %% max(1, floor(max(1, gnl_l) / 20)) == 0 || m == gnl_l) {
                         update_annotation_progress(
                           progress = 70 + (m / max(1, gnl_l)) * 25,
                           detail = paste0("Applying GO annotations to rows (", m, "/", gnl_l, ")...")
                         )
                       }
                     }
                     update_annotation_progress(progress = 96, detail = "Writing GO columns to dataset...")
                     updates <- tibble::tibble(
                       row_id = current_data()$row_id,
                       GO_MF = mf_add,
                       GO_BP = bp_add,
                       GO_CC = cc_add
                     )
                     update_derived_rows(updates)
                     update_annotation_progress(
                       progress = 100,
                       status = paste0("GO annotations added successfully at ", format(Sys.time(), "%H:%M:%S"), "."),
                       detail = "Completed."
                     )
                   },
                   error = function(e) {
                     update_annotation_progress(
                       progress = 0,
                       status = paste("GO annotation failed:", e$message),
                       detail = "Stopped."
                     )
                     showNotification(paste("GO annotation failed:", e$message), type = "error")
                   }
                 )
               }
  )
  
  
}

# Create Shiny app  
shinyApp(ui, server)
