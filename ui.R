# ui.R

library(shiny)

shinyUI(fluidPage(
  titlePanel("QTL Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("traits", label = h5("Plot which trait?"), 
                  choices = list("Trait 1" = 1,
                                 "Trait 2" = 2,
                                 "Trait 3" = 3),
                  selected = 1),
      br(),
      selectInput("chromosome", label = h5("Plot which chromosome(s)?"), 
                  choices = list("All" = 0,
                                 "Chromosome 1" = 1,
                                 "Chromosome 2" = 2,
                                 "Chromosome 3" = 3,
                                 "Chromosome 4" = 4,
                                 "Chromosome 5" = 5),
                  selected = 0),
      conditionalPanel(condition = "input.chromosome != 0", uiOutput("slider"))

    ),
    
    mainPanel(
      plotOutput("graph")
    )
  )
    
))