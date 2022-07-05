#
# Application to explore differential gene expression in 
# bulk RNAseq data comparing uninjured caudal fin tissue to regenerating tissue
#
# Allows user to select gene of interest (searchable from dropdown list), 
# returns expression information, and plot of normalized counts from DESeq
#

library(shiny)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(data.table)

# Define non-reactive elements
data <- fread("data/Poss_RNAseq_DESeq_output.txt")
data <- data %>% .[,-1] 
colnames(data) <- c("ENSEMBL ID", "Base Mean", "log2 Fold Change", "Shrunken Log Fold Change Estimates", "p Value", "Adjusted p Value", "Threshold", "Gene Name")
#deseq_obj <- load("data/Poss_RNAseq_dds.rds")

# Define UI for application 
ui <- fluidPage(
    
    # Application title
    titlePanel("Regenerating Caudal Fin Tissue: Uninjured vs. 4DPA"),
    
    #User selects gene name from external gene name list 
    selectizeInput("gene", label = "Enter gene symbols of interest:", choices = NULL, multiple = TRUE),
    
    #output controls that tell Shiny where to put rendered output (text and plot)
    textOutput("text"),
     tableOutput("static")
    #plotOutput("plot", width = "400px")
    
    )

# Define server logic for generating summary stats
server <- function(input, output, session) {
    
    gene <- reactive({
        c(input$gene)
    })
    
    updateSelectizeInput(session, "gene", choices = data$`Gene Name`, server = TRUE)
    
    output$text <- renderText({
        "Gene-level summary stats from DESseq2 analysis:"
    })
    
    output$static <- renderTable({
        data[data$`Gene Name`%in%gene(),]
    })
    
   # output$table <- renderTable({
        #dataset()
   # })
    
    

}

# Run the application 
shinyApp(ui = ui, server = server)
