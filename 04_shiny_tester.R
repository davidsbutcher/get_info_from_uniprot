library(tidyverse)
library(shiny)
library(magrittr)

ui <- fluidPage(
  
  titlePanel("Hello Shiny!"),
  
  sidebarLayout(
    
    sidebarPanel(
      textInput("filedir", "File Path:"),
      numericInput("fdrcutoff", "False Discovery Rate Cutoff:", value = 0.01, min = 0.001, max = 1, step = 0.001),
      numericInput("UPtaxon", "UniProt Taxon:", value = 83333),
      sliderInput("futureworkers", "Number of R Sessions for Future:", min = 0, max = 10, value = 10, step = 1),
      actionButton("runbutton", label = "Run")
    ),
    
    mainPanel(
      textOutput("selected_dir"),
      textOutput("fdr"),
      textOutput("UPtaxon")
      )
    
  )
  
)

server <- function(input, output, session) {
  
  output$selected_dir <-
    
    renderText({ 
      paste("Input Directory: ", input$filedir)
    })

  output$fdr <-
    
    renderText({ 
      paste("FDR: ", input$fdrcutoff)
    })
  
  output$UPtaxon <-
    
    renderText({ 
      paste("UniProt Taxon to search: ", input$UPtaxon)
    })
  
  observeEvent(input$runbutton, {
    
    filedir <- input$filedir
    
    fdr <- input$fdrcutoff
    
    UPtaxon <- UniProt.ws::UniProt.ws(input$UPtaxon)
    
    go_locs <- read_tsv("QuickGO_annotations_20190708.tsv") %>%
      .["GO NAME"] %>% unique() %>% pull()
    
    plan(multisession(workers = input$futureworkers))
    
    tic()
    
    source("01_get_protein_info.R")
    source("02_get_proteoform_info.R")
    source("03_output_plots.R")
    
    totaltime <- capture.output(toc()) %>%
      str_extract("[0-9]+") %>%
      as.numeric %>% 
      `/`(60) %>% round(digits = 2)
    
    print(glue("Elapsed time: {totaltime} min"))
    
  })
  
}

shinyApp(ui, server)