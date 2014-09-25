# server.R

library(ggplot2)
library(grid)

qtl <- read.csv("data/traits.csv")

qtl$background.color <- "1"
qtl$background.color[(qtl$chr %% 2) == 0] <- "0"

shinyServer(function(input, output) {
  
  # slider input
  output$slider <- renderUI({
    qtlChr <- subset(qtl, chr == input$chromosome)
    sliderInput("region", label = h5("Display a region of the chromosome?"),
     min = 1, max = max(qtlChr$pos), value = c(1, max(qtlChr$pos)))
  })
                           
  # plots the graph
  output$graph <- renderPlot({
    
    # subsets the data depending on the trait selected
    if (input$traits == 1){
      qtlData <- subset(qtl, select = c(chr, pos, trait1_lod, background.color))
    } else if (input$traits == 2){
      qtlData <- subset(qtl, select = c(chr, pos, trait2_lod, background.color))
    } else {
      qtlData <- subset(qtl, select = c(chr, pos, trait3_lod, background.color))
    }
    
    # subsets the data depending on chromosome selected
    if (input$chromosome == 0){
      # do nothing
    } else if (input$chromosome == 1){
      qtlData <- subset(qtlData, chr == 1)
    } else if (input$chromosome == 2){
      qtlData <- subset(qtlData, chr == 2)
    } else if (input$chromosome == 3){
      qtlData <- subset(qtlData, chr == 3)
    } else if (input$chromosome == 4){
      qtlData <- subset(qtlData, chr == 4)
    } else {
      qtlData <- subset(qtlData, chr == 5)
    }
    
    # changes the column name of traitx_lod to lod
    colnames(qtlData)[3] <- "lod"
    
    # extracts the max lod value for the data set
    peak <- max(qtlData$lod)
    
    # start of the plot
    p <- ggplot(qtlData)
    
    # determines the color of the background; depends if all chromosomes are plotted, or just one
    # plots a certain region only for single chromosomes (changes xlim)
    if (input$chromosome == 0){
      p <- p +
        scale_fill_manual(values = c("1" = "gray90", "0" = "white")) +
        scale_x_continuous(breaks = seq(0, 1000, by = 25))
    } else {
      p <- p +
        scale_fill_manual(values = c("1" = "white", "0" = "white")) +
        scale_x_continuous(breaks = seq(0, 1000, by = 25),
         limits = c(input$region[1], input$region[2]))
    }

    # rest of the plot
    p <- p +
      facet_grid(~ chr, scales = "free_x", space = "free_x") +
      geom_rect(aes(fill = background.color), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_line(aes(x = pos, y = lod), size = 2) +
      geom_hline(yintercept = 0.50, color = "red", size = 1) +
      geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
      scale_y_continuous(expand = c(0, 0), limits = c((peak * -0.06), (peak * 1.02))) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90),
            axis.line=element_line(),
            panel.margin = unit(0, "cm")) +
      ggtitle("LOD Curves for QTLs") +
      xlab("Position in cM") +
      ylab("LOD Score") 

    # prints the plot
    print(p)
  })
  
})
