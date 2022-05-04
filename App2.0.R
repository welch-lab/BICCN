
#Shiny App 2.0
#Load Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(rliger)
library(shiny)
library(gplots)
library(shinythemes)

#Import Data
setwd("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/MarkerGeneBrowser/Base_Reference_Files/")

#Please input your own prior to running:
dataset.vector = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region/AUD/Analysis1_AUD/Images/GeneExpression_AUD_Analysis1.RDS")
results.table = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region/AUD/Analysis1_AUD/Analysis1_AUD_Results_Table.RDS")

ui <- fluidPage(
  navbarPage("Marker Gene Browser", theme = shinytheme("lumen"), 
             tabPanel("Marker Genes", fluid = TRUE,
                      sidebarLayout(
                        sidebarPanel (width =2,
                                      titlePanel("Select Datasets & Markers"),
                                      fluidRow(          selectInput(inputId = "Dataset", 
                                                                     label = "Dataset:",
                                                                     choices= names(dataset.vector),
                                                                     width = '400px'),
                                                         selectInput(inputId = "MarkerGenes", 
                                                                     label = "Marker Genes",
                                                                     choices = c("meth", "atac", "rna")),
                                                         selectInput(inputId = "CellSubtype",
                                                                     label = "CellSubtype",
                                                                     choices = c("Exc", "Inh", "NonN")),
                                                         selectInput(inputId = "Resolution", 
                                                                     label = "High Resolution?",
                                                                     choices = c("lowRcluster", "highRcluster"))
                                                         
                                      )
                                      
                        ),
                        mainPanel(
                          fluidRow(column(10, plotOutput("dotplots")))
                          
                        )
                      )),
             tabPanel("Associated Cell Types", fluid = TRUE,
                      sidebarLayout(
                        sidebarPanel (width =2,
                                      titlePanel("Select Datasets & Markers"),
                                      fluidRow(          selectInput(inputId = "Dataset", 
                                                                     label = "Dataset:",
                                                                     choices= names(dataset.vector),
                                                                     width = '400px'),
                                                         selectInput(inputId = "MarkerGenes", 
                                                                     label = "Marker Genes",
                                                                     choices = c("meth", "atac", "rna")),
                                                         selectInput(inputId = "CellSubtype",
                                                                     label = "CellSubtype",
                                                                     choices = c("Exc", "Inh", "NonN")),
                                                         selectInput(inputId = "Resolution", 
                                                                     label = "High Resolution?",
                                                                     choices = c("lowRcluster", "highRcluster"))
                                                         
                                      )
                                      
                        ),
                        mainPanel(
                          fluidRow(column(10, plotOutput("markers"))
                                   
                                   
                          )
                        )
                        
                      )
             ),
             tabPanel("Cell Types", fluid = TRUE,
                      fluidRow(
                        selectInput(inputId = "ResolutionC", 
                                    label = "Resolution?",
                                    choices = c("lowRcluster", "highRcluster"))),
                      numericInput(inputId = "Cluster",
                                   label = "Select Cluster", value = 0, min = 0, max = NA,step = NA),
                      fluidRow(
                        
                        plotOutput("celltypes"))
                      
             )
  )
)







server <- function(input, output) {
  output$markers = renderPlot(
    height = 1000,
    width = 1000,
    {
      meow = dotplot_generator(Dataset = input$Dataset, markerGenes = input$MarkerGenes, CellSubtype = input$CellSubtype, resolution = input$Resolution)
      marks = marker_table(meow,markerGenes = input$MarkerGenes, CellSubtype = input$CellSubtype, resolution = input$Resolution)
      marks
    })
  output$dotplots = renderPlot(
    width = 800,
    height = 1000,
    {
      test = dotplot_generator(Dataset = input$Dataset, markerGenes = input$MarkerGenes, CellSubtype = input$CellSubtype, resolution = input$Resolution)
      DotPlot(test, resolution = input$Resolution)
    })
  
  output$celltypes = renderPlot( 
    width = 1800,
    height = 600,
    {
      types = graph_celltypes(input$ResolutionC, input$Cluster)
      types})
}

shinyApp(ui, server)

