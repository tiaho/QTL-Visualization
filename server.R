# server.R

library(ggplot2)
library(grid)

qtl <- read.csv("data/traits.csv")

qtl$background.color <- "1"
qtl$background.color[(qtl$chr %% 2) == 0] <- "0"

### makes the test data set for plotting gene expression
gene.name <- vector()
for (i in 1:50){
  gene.name[i] = paste("gene", i, sep = "")
}
expression.values <- rnorm(50, 0 ,1)
genetic.position <- runif(50, 0, 100)
physical.position <- genetic.position + 10000
chrom <- vector()
for (i in 1:5){
  chrom <- c(chrom, rep(i, 10))
}
expression.data <- as.data.frame(cbind(gene.name, expression.values, genetic.position, physical.position, chrom))

# converts the values from factor to numeric
expression.data$expression.values <- as.numeric(as.character(expression.data$expression.values))
expression.data$genetic.position <- as.numeric(as.character(expression.data$genetic.position))
expression.data$physical.position <- as.numeric(as.character(expression.data$physical.position))
expression.data$chrom <- as.numeric(as.character(expression.data$chrom))


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



shinyServer(function(input, output) {
  
  # slider input
  output$slider <- renderUI({
    qtlChr <- subset(qtl, chr == input$chromosome)
    max_pos <- ceiling(max(qtlChr$pos, 1))
    sliderInput("region", label = h5("Display a region of the chromosome?"),
     min = 0, max = max_pos, value = c(0, max_pos), step = 1)
  })
  
  # sort checkbox input
  output$sort_checkbox <- renderUI({
    checkboxInput("sort", "Sorts the expression value from most negative to most positive", value = FALSE)
  })
  
  
  
  # Insert the right number of plot output objects into the web page
  output$plots <- renderUI({
    plot_output_list <- lapply(1:2, function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname)
#       plotOutput(plotname, height = 280, width = 250)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (i in 1:1) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")
      
      output[[plotname]] <- renderPlot({
        # qtl graph
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
          qtl_plot <- ggplot(qtlData)
          
          # determines the color of the background; depends if all chromosomes are plotted, or just one
          # plots a certain region only for single chromosomes (changes xlim)
          if (input$chromosome == 0){
            qtl_plot <- qtl_plot +
                          scale_fill_manual(values = c("1" = "gray90", "0" = "white")) +
                          scale_x_continuous(breaks = seq(0, 1000, by = 25))
          } else {
            qtl_plot <- qtl_plot +
                          scale_fill_manual(values = c("1" = "white", "0" = "white")) +
                          scale_x_continuous(breaks = seq(0, 1000, by = 25),
                            limits = c(input$region[1], input$region[2]))
          }
          
          # rest of the plot
          qtl_plot <- qtl_plot +
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
        
          
          
        # expression graph
        
        # subsets the data depending on chromosome selected
        data4plotting <- subset(expression.data, chrom == input$chromosome)
        
        # plots either the genetic distance or physical distance
        # plots genetic distance
        if (input$distance == 1){
          if (input$sort == TRUE){
            data4plotting$genetic.position <- factor(data4plotting$genetic.position, levels = data4plotting[order(data4plotting$expression.values, decreasing = TRUE), "genetic.position"])
          }
          expression_plot <- ggplot(data4plotting) +
            scale_x_continuous(limits = c(input$region[1], input$region[2])) +
            geom_histogram(aes(genetic.position, expression.values), stat="identity") +
            coord_flip() +
            xlab("Genetic Position") +
            ylab("Relative Gene Expression")
          
        }
        # plots physical distance
        else if(input$distance == 2){
          if (input$sort == TRUE){
            data4plotting$physical.position <- factor(data4plotting$physical.position, levels = data4plotting[order(data4plotting$expression.values, decreasing = TRUE), "physical.position"])
          }
          expression_plot <- ggplot(data4plotting) +
            scale_x_continuous(limits = c(input$region[1], input$region[2])) +
            geom_histogram(aes(physical.position, expression.values), stat="identity") +
            coord_flip() +
            xlab("Physical Position") +
            ylab("Relative Gene Expression")
        }
        
      # plots the plots
      multiplot(qtl_plot, expression_plot, cols=1)
      
      })
    })
  }
})
