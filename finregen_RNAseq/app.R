#
# Application to explore differential gene expression in 
# bulk RNAseq data comparing uninjured caudal fin tissue to regenerating tissue
#
# Allows user to select gene of interest (searchable from dropdown list), 
# returns expression information and plot of normalized counts between time points
#

library(shiny)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(data.table)

# Define non-reactive elements

#load and wrangle the DEseq output file 
data <- fread("data/Poss_RNAseq_DESeq_output.txt")
data <- data %>% .[,-1] %>% as.data.frame()
colnames(data) <- c("ENSEMBL ID", "Base Mean", "log2 Fold Change", "Shrunken Log Fold Change Estimates", "p Value", "Adjusted p Value", "Threshold", "Gene Name")

#load metadata and expression matrix, run DESeq2 
metadata <- fread("data/Poss_metadata.txt")
metadata <- metadata %>% as.data.frame()
rownames(metadata) <- metadata$V1
metadata <- metadata[,-1] %>% as.matrix()
matrix <- fread("data/Poss_count_matrix.txt")
matrix <- matrix %>% as.data.frame()
rownames(matrix) <- matrix$V1
matrix <- matrix[,-1] %>% as.matrix()

dds <- DESeqDataSetFromMatrix(countData=matrix, 
                              colData=metadata, 
                              design=~dpa)
#specify the order of dpa factor
dds$dpa <- factor(dds$dpa, levels = c("uninjured","4DPA"))
dds <- DESeq(dds)


# Define UI for application 
ui <- fluidPage(
    
    # Application title
    titlePanel("Regenerating Caudal Fin Tissue: Uninjured vs. 4DPA"),
    
    #User selects gene name from external gene name list 
    selectizeInput("gene", label = "Enter gene symbols of interest:", choices = NULL, multiple = FALSE),
    
    #output controls that tell Shiny where to put rendered output (text and plot)
    textOutput("text"),
    tableOutput("static"),
    plotOutput("plot")
    
    )

# Define server logic for generating summary stats and plot
server <- function(input, output, session) {
    
    updateSelectizeInput(session, "gene", choices = data$`Gene Name`, server = TRUE)
    
    gene <- reactive({
        c(input$gene)
    })
    
    ensembl_id <- reactive({
      c(data[data$`Gene Name`%in%gene(), 1])
    })
    
    output$text <- renderText({
        "Gene-level summary stats from DESeq2 analysis:"
    })
    
    output$static <- renderTable({
        data[data$`Gene Name`%in%gene(), 1:7]
    })

    output$plot <- renderPlot({
        ensembl_id <- c(as.character(data[data$`Gene Name`%in%gene(), 1]))
        plotCounts(dds, gene = ensembl_id, intgroup = "dpa", xlab = "Timepoint")
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
