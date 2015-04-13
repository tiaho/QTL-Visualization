# ui.R

library(shiny)

# creates the data used for the trait select input
phenotype_data <- read.csv("data/simulated_phenotypes_real_map.csv")
num_traits <- length(grep("trait", names(phenotype_data)))
trait_names <- list()
for (i in 1:num_traits){
  trait_names[i] = i
}
for (i in 1:num_traits){
  names(trait_names)[i] = paste("Trait", i, sep = " ")
}

  
shinyUI(fluidPage(
  titlePanel("QTLVizR"),
  
  sidebarLayout(
    sidebarPanel(     
      selectizeInput("traits", label = h5("Choose up to 3 traits"),
                      multiple = TRUE, selected = "1",
                      choices = trait_names, 
      ),
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
                                 "fold change" = 2),
                  selected = 1),
      br(),
      conditionalPanel(condition = "input.chromosome != 0", uiOutput("slider")),
      br(),
      downloadButton('download_table', 'Download Full Gene Table'),
      br(),
      br(),
      br(),
      "QTLvisR is a web application to help visualize quantitative trait loci (QTL).
      This app allows the user to select physiological and developmental traits from a drop down menu.
      The app displays the support for a physiological or developmental (pQTL) for every genotyped marker across the organism's genome.
      The user can sub-select specific chromosomes, or intervals along a chromosome to focus attention on.
      At the same time, QTL's for gene expression (eQTL) for genes falling within the user defined genomic interval are also visualized.
      This allows the user to quickly examine a genomic region for co-occurrence of pQTL and eQTL."
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