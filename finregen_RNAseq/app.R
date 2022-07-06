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
results <- results(dds)

#specify size of plots
size = 450
#specify max y value for MA plot 
ymax <- 5

# Define UI for application 
ui <- fluidPage(
    
    # Application title
    titlePanel("Regenerating caudal fin tissue: bulk RNAseq"),
    
    sidebarLayout(
      sidebarPanel(
        #User selects gene name from external gene name list 
        selectizeInput("gene", label = "Enter gene symbols of interest:", choices = NULL, multiple = FALSE),
        
      ),
    
    mainPanel(
    #output controls that tell Shiny where to put rendered output (text and plot)
    #desriptive text above stat table
    textOutput("text"),
    #table of stats from DEseq2
    tableOutput("static"),
    #plot of normalized counts for selected gene 
    plotOutput("counts_plot", width = size, height = size, click = "counts_plot_click"),
    #sample-level count information (based on clicked point)
    textOutput("counts_info"),
    #text title for MA plot
    verbatimTextOutput("ma_plot_title"),
    #MA plot - with adjustable slider and clickable points 
    plotOutput("ma_plot", width = size, height = size),
    #MA plot slider 
    sliderInput("alpha", "Adjusted p-value threshold", min = 0, max = 0.5, value = 0.1, step = 0.001, width = size)
    )
    
  )
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

    output$counts_plot <- renderPlot({
        ensembl_id <- c(as.character(data[data$`Gene Name`%in%gene(), 1]))
        p <- plotCounts(dds, gene = ensembl_id, intgroup = "dpa", returnData = TRUE)
        ggplot(p, aes(x=dpa, y=count)) +
          geom_point(position = position_jitter(w=0.1, h=0)) +
          scale_y_log10(breaks=c(25, 100, 400)) +
          xlab("Timepoint") +
          theme_classic(base_size = 20)
    })
    
    output$counts_info <- renderText({
      counts_str <- function(e) {
        if(is.null(e)) return("NULL \n")
        paste0("Normalized counts: ", round(e$y, 1), "\n")
        
      }

    paste0(counts_str(input$counts_plot_click)) 
    })

    output$ma_plot <- renderPlot({
      par( mar = c(5, 5, 3, 2), cex.main = 1.5, cex.lab = 1.35)
      plotMA(results, ylim = c(-ymax, ymax), alpha = input$alpha)
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
