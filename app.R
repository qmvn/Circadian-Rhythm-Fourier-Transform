#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(tidyverse)
library(dplyr)
library(dryR)
library(janitor)
library(data.table)
library(leaflet)


#loading datasets
load("C:/Users/quynh/OneDrive/Desktop/Summer research/research/simple.RData")

load("C:/Users/quynh/OneDrive/Desktop/Summer research/research/dataset.RData")

load("C:/Users/quynh/OneDrive/Desktop/Summer research/research/sim.RData")
dryListsim

source("C:/Users/quynh/OneDrive/Desktop/Summer research/research/dry_plot_harmonic.R")
#loading `gene` information
dry_plot_harmonic(dryListsim, "sin1")

view(dryListsim[["results"]])


parameters <- dryListharmonic[["parameters"]]

model_1 <- parameters[parameters[, "chosen_model"] == 1,]

model_2 <- parameters[parameters[, "chosen_model"] == 2,]

model_3 <- parameters[parameters[, "chosen_model"] == 3,]

model_4 <- parameters[parameters[, "chosen_model"] == 4,]

model_5 <- parameters[parameters[, "chosen_model"] == 5,]

length1 <- as.numeric(nrow(model_1))
length2 <- as.numeric(nrow(model_2))
length3 <- as.numeric(nrow(model_3))
length4 <- as.numeric(nrow(model_4))
length5 <- as.numeric(nrow(model_5))

dataset <- data.frame(
    rn = rownames(dryListharmonic[["parameters"]]),
    chosen_model = dryListharmonic[["parameters"]][, "chosen_model"],
    bicharmonic = dryListharmonic[["parameters"]][, "chosen_model_BICW"],
    bicsimple = simpledryList[["parameters"]][, "chosen_model_BICW"]
    
)

bicyay <- as.numeric(nrow(dataset[dataset$bicharmonic > dataset$bicsimple,]))

bicnay <- as.numeric(nrow(dataset[dataset$bicharmonic <= dataset$bicsimple,]))

#tibfull <- setDT(parameters, keep.rownames = TRUE)[]

#tib <- tibfull[,c("rn","chosen_model")]


# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("DryR models of the sample data given"),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("select_model",
                        "Select a model:",
                        choices = c("Neither conditions cycling" = 1, "Cond 1 (18C) cycling, Cond 2 (25C) not cycling" = 2, "Cond 1 (18C) not cycling, Cond 2 (25C) cycling" = 3, "Both conditions cycling with the same differential model" = 4, "Both conditions cycling with different differential models"=5), selected = 1
            ),
            selectInput("BICW", 
                               strong("BICW scores"), 
                               choices = c("12/24 > 24 model BIC weight" = 1, 
                                              "24 > 12/24 model BIC weight" = 2), 
                                             
                               selected = 2),
            selectInput("select_gene",
                        "Select a gene:",
                        choices = c(rownames(dryListharmonic[["parameters"]])[1])
            ),
            p("BIC Weight of 12 and 24 hour harmonic model: ", htmlOutput("bicw12")),
            p("BIC Weight of 24 hour harmonic model: ", htmlOutput("bicw24")),
            p("Genes with neither gene cycling: ", strong(length1)),
            br(),
            strong("Summary:"),
            p("Genes where 12/24 hr modeling was more abundant than 24 hr modeling: ", strong(bicyay)),
            p("Genes where 24 hr modeling was more abundant than 12/24 hr modeling: ", strong(bicnay)),
            p("Genes with neither condition cycling: ", strong(length1)),
            p("Genes with Cond 1 (18C) cycling, Cond 2 (25C) not cycling: ", strong(length2)),
            p("Genes with Cond 1 (18C) not cycling, Cond 2 (25C) cycling: ", strong(length3)),
            p("Genes with both conditions cycling with the same differential model: ", strong(length4)),
            p("Genes with both conditions cycling with different differential models: ", strong(length5))
            
            
            
            
        ),
        # Show a plot of the generated distribution
        mainPanel(
            h3("Modeling with both 12 and 24 hour harmonics"),
            plotOutput("dryRharmonic"),
            h3("Modeling with a 24 hour harmonic"),
            plotOutput("dryR"),
            h3("Fit Formula: "),
            withMathJax(
                helpText("\\(\\mu*y\\) + \\(ay*cos(2\\pi*\\frac{x}{24})\\)  + \\(by*sin(2\\pi*\\frac{x}{24})\\) + \\(cy*cos(2\\pi*\\frac{x}{12})\\) + \\(dy*sin(2\\pi*\\frac{x}{12})\\)")
            ),
            p("Where \\(x\\) represents the time point of transform (ZT), \\(y\\) represents the specific gene expression amount at a given time (ZT), \\(\\mu\\) represents the average level of expression for a given gene, and \\(a\\), \\(b\\), \\(c\\), and \\(d\\) are various phase [arctan\\(\\frac{a}{b}\\) and arctan\\(\\frac{c}{d}\\)] and amplitude coefficients [log2-fold change peak-to-trough; 2sqrt{\\(a^2+b^2\\)} and $$\\sqrt{2}$$
{\\(c^2+d^2\\)}] for a given gene expression pattern."),
            h3("Statistics"),
            h4("Gene name:", strong(textOutput("result", inline = TRUE))), 
            br(), 
            h4(strong("Parameters with 12 and 24 hr Harmonics: "), br()),
            p(tableOutput("parameters")),
            p(tableOutput("parameters2")),
            p(tableOutput("parameters3")),
            
            h4(strong("Parameters with a 24 hr Harmonic: "), br()),
            p(tableOutput("sparameters")),
            p(tableOutput("sparameters2")),
            
            tabsetPanel(
                tabPanel("Synthetic Data", leafletOutput("synthetic"),
                         p("Modeling with 12 and 24 hour harmonics"),
                         plotOutput("modelharmonic"),
                         p("Modeling with 24 hour harmonic"),
                         plotOutput("modelsimple"),
                         selectInput("syntheticmod", "Select a model", 
                                     choices = c("Neither conditions cycling" = 1, "Cond 1 (18C) cycling, Cond 2 (25C) not cycling" = 2, "Cond 1 (18C) not cycling, Cond 2 (25C) cycling" = 3, "Both conditions cycling with the same differential model" = 4, "Both conditions cycling with different differential models"=5), selected = 1),
                         selectInput("select_modgene",
                                     "Select a gene:",
                                     choices = c(rownames(dryListharmonic[["parameters"]])[1]))
))
            
        )
    )
)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observe({
        data <- data.frame(
            rn = rownames(dryListharmonic[["parameters"]]),
            chosen_model = dryListharmonic[["parameters"]][, "chosen_model"],
            bicharmonic = dryListharmonic[["parameters"]][, "chosen_model_BICW"],
            bicsimple = simpledryList[["parameters"]][, "chosen_model_BICW"]
            
        )
        if(input$BICW == 1){
            selected_genes <- data %>%
                filter(chosen_model == as.numeric(input$select_model) & bicharmonic >= bicsimple) %>%
                pull(rn) %>%
                as.character()
        }else{
            selected_genes <- data %>%
                filter(chosen_model == as.numeric(input$select_model) & bicharmonic < bicsimple) %>%
                pull(rn) %>%
                as.character()
        }
        # Set to blank if there are no gene names
        if (is.null(selected_genes)) {
            selected_genes <- character(0)
        }
        
        
        
        # Set the labels to the new gene names
        updateSelectInput(session, "select_gene",
                          label = paste("Select a gene:"),
                          choices = selected_genes,
                          selected = head(selected_genes, 1)
        )
    })
    output$dryRharmonic <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the plot with the specified gene
        dry_plot_harmonic(dryList = dryListharmonic, gene = input$select_gene)
    })
    
    output$dryR <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the plot with the specified gene
        dry_plot(dryList = simpledryList, gene = input$select_gene)
    })
    
    output$result = renderText(as.character(input$select_gene))
    
    subset <- dryListharmonic[["parameters"]]
    rownames(subset)
    output$parameters = renderTable(subset[input$select_gene,1:10])
    output$parameters2 = renderTable(subset[input$select_gene,11:20])
    output$parameters3 = renderTable(subset[input$select_gene,21:26])
    
    subset2 <- simpledryList[["parameters"]]
    rownames(subset)
    output$sparameters = renderTable(subset2[input$select_gene,1:10])
    output$sparameters2 = renderTable(subset2[input$select_gene,11:16])
    

    output$bicw12 <- renderUI({
        subset[input$select_gene, "chosen_model_BICW"] <- paste("<b>",subset[input$select_gene, "chosen_model_BICW"],"</b>")
        HTML(paste(subset[input$select_gene, "chosen_model_BICW"]))
    })
    
    output$bicw24 <- renderUI({
        subset2[input$select_gene, "chosen_model_BICW"] <- paste("<b>",subset2[input$select_gene, "chosen_model_BICW"],"</b>")
        HTML(paste(subset2[input$select_gene, "chosen_model_BICW"]))
    })

    
    output$synthetic <- renderLeaflet({
        
        observe({
            datasim <- data.frame(
                rnmod = rownames(dryListsim[["parameters"]]),
                chosen_mod = dryListsim[["parameters"]][, "chosen_model"]
                
            )
            selected_modgenes <- datasim %>%
                filter(chosen_mod == as.numeric(input$syntheticmod)) %>%
                pull(rnmod) %>%
                as.character()
            # Set to blank if there are no gene names
            if (is.null(selected_modgenes)) {
                selected_modgenes <- character(0)
            }
            
            
            
            # Set the labels to the new gene names
            updateSelectInput(session, "select_modgene",
                              label = paste("Select a gene:"),
                              choices = selected_modgenes,
                              selected = head(selected_modgenes, 1)
            )
        })
        
        
    })

    output$modelharmonic <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the plot with the specified gene
        dry_plot_harmonic(dryList = dryListsim, gene = input$select_modgene)
    })
    
    output$modelsimple <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the plot with the specified gene
        dryR::dry_plot(dryList = dryListsim, gene = input$select_modgene)
    })

}


# Run the application 
shinyApp(ui = ui, server = server)

