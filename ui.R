# ui.R

library(shiny)

shinyUI(fluidPage(
  titlePanel("QTL Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("traits", label = h5("Plot which trait?"), 
                  choices = list("Trait 1" = "trait1",
                                 "Trait 2" = "trait2",
                                 "Trait 3" = "trait3",
                                 "Traits 1 and 2" = 4),
                  selected = 1),
      br(),
      selectInput("chromosome", label = h5("Plot which chromosome(s)?"), 
                  choices = list("All" = 0,
                                 "Chromosome 1" = 1,
                                 "Chromosome 2" = 2,
                                 "Chromosome 3" = 3,
                                 "Chromosome 4" = 4,
                                 "Chromosome 5" = 5,
                                 "Chromosome 6" = 6,
                                 "Chromosome 7" = 7,
                                 "Chromosome 8" = 8,
                                 "Chromosome 9" = 9,
                                 "Chromosome 10" = 10),
                  selected = 0),
      br(),
      selectInput("ex_graph", label = h5("Plot t-statistic or fold change"), 
                  choices = list("t-statistic" = 1,
                                 "Fold Change" = 2),
                  selected = 1),
      br(),
      conditionalPanel(condition = "input.chromosome != 0", uiOutput("slider")),
      br(),
      downloadButton('download_table', 'Download Full Gene Table')
    ),
    
    mainPanel(
      plotOutput("qtl_graph"),
      br(),
      br(),
      conditionalPanel(condition = "input.chromosome != 0", plotOutput("expression_graph")),
      br(),
      br(),
      tableOutput("table")
    )
  )
    
))