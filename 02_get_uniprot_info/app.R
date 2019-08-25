library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Get Uniprot data for TDReport Proteins"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(fileInput(filelist, multiple = TRUE),
                     textInput(uniprot_taxon, "Uniprot Taxon #", placeholder = "83333")
                     ),
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
