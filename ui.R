# ui.R

library(shiny)

shinyUI(fluidPage(
  titlePanel("QTL Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("traits", label = h5("Plot which trait?"), 
                  choices = list("Trait 1" = 1,
                                 "Trait 2" = 2,
                                 "Trait 3" = 3,
                                 "Traits 1 and 2" = 4),
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
      br(),
      selectInput("distance", label = h5("Plot genetic or physical distance?"), 
                  choices = list("Genetic Distance" = 1,
                                 "Physical Distance" = 2),
                  selected = 1),
      br(),
      conditionalPanel(condition = "input.chromosome != 0", uiOutput("slider")),
      br(),
      conditionalPanel(condition = "input.chromosome !=0", uiOutput("sort_checkbox"))
    ),
    
    mainPanel(
      uiOutput("plots")
    )
  )
    
))