# server.R

library(ggplot2)
library(reshape2)
library(grid)

# qtl <- read.csv("data/traits.csv")
qtl <- read.csv("data/allele_specific_test_p_adjusted.csv")

qtl$background.color <- "1"
qtl$background.color[(qtl$chr %% 2) == 0] <- "0"

num_genes = length(qtl$gene_name)

### makes the test data set for plotting gene expression
qtl$genetic.position <- runif(num_genes, 0, 100)
qtl$genetic.position <- round(qtl$genetic.position, 1)
qtl$physical.position <- qtl$genetic.position + 10000
qtl$chr <- ceiling(runif(num_genes, 0, 5))
qtl$trait1 <- qtl[,2]
qtl$trait2 <- qtl[,3]
# qtl$expression.values <- rnorm(num_genes, 0 ,1)
qtl$A <- rnorm(num_genes, 0 ,1)
qtl$B <- rnorm(num_genes, 0 ,1)

# converts the values from factor to numeric
# expression.data$expression.values <- as.numeric(as.character(expression.data$expression.values))
# expression.data$A <- as.numeric(as.character(expression.data$A))
# expression.data$B <- as.numeric(as.character(expression.data$B))
# expression.data$genetic.position <- as.numeric(as.character(expression.data$genetic.position))
# expression.data$physical.position <- as.numeric(as.character(expression.data$physical.position))
# expression.data$chrom <- as.numeric(as.character(expression.data$chrom))

# changes data frame from wide format to long format
qtl <- melt(qtl, id.vars = c("gene_name", "trait1", "trait2", "chr", "genetic.position", "physical.position", "background.color"),
                              measure.vars = c("A", "B"),
                              variable.name = "allele",
                              value.name = "expression.value")

qtl <- melt(qtl, id.vars = c("gene_name", "chr", "genetic.position", "physical.position", "allele", "expression.value", "background.color"),
            measure.vars = c("trait1", "trait2"),
            variable.name = "trait",
            value.name = "lod")

shinyServer(function(input, output) {
  
  # slider input
  output$slider <- renderUI({
    qtlChr <- subset(qtl, chr == input$chromosome)
    max_pos <- ceiling(max(qtlChr$genetic.position, 1))
    sliderInput("region", label = h5("Display a region of the chromosome?"),
     min = 0, max = max_pos, value = c(0, max_pos), step = 1)
  })

  # qtl graph
  output$qtl_graph <- renderPlot({
    # subsets the data depending on the trait selected
    if (input$traits == 1){
      qtlData <- subset(qtl, trait == "trait1", select = c(chr, genetic.position, lod, background.color))
    } else if (input$traits == 2){
      qtlData <- subset(qtl, trait == "trait2", select = c(chr, genetic.position, lod, background.color))
#     } else if (input$traits == 3){
#       qtlData <- subset(qtl, select = c(chr, genetic.position, trait3_lod, background.color))
    } else {
      qtlData <- subset(qtl, select = c(chr, genetic.position, lod, background.color))
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
    
#     # changes the column name of traitx_lod to lod
#     colnames(qtlData)[3] <- "lod"
#     if (input$traits == 4){
#       colnames(qtlData)[4] <- "lod2"
#     }          
#     
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
      geom_line(aes(x = genetic.position, y = lod), size = 2) +
      geom_hline(yintercept = 0.50, color = "red", size = 1) +
      geom_segment(aes(x = genetic.position, xend = genetic.position), y = (peak * -0.02), yend = (peak * -0.05)) +
      scale_y_continuous(expand = c(0, 0), limits = c((peak * -0.06), (peak * 1.02))) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90),
            axis.line=element_line(),
            panel.margin = unit(0, "cm")) +
      ggtitle("LOD Curves for QTLs") +
      xlab("Position in cM") +
      ylab("LOD Score") +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    
    if (input$traits == 4){
      qtl_plot +
        geom_line(aes(x = genetic.position, y = lod2), size = 2, color = "blue")
    } else {
      qtl_plot
    }
  })
  
  # expression graph
  output$expression_graph <- renderPlot({
    # subsets the data depending on chromosome selected
    if (input$traits == 4){
      data4plotting <- subset(qtl, chr == input$chromosome & (trait == "trait1" | trait == "trait2"))
    } else {
      data4plotting <- subset(qtl, chr == input$chromosome & trait == input$traits)
    }
    
    # plots either the genetic distance or physical distance
    # plots genetic distance
    if (input$distance == 1){
      expression_plot <- ggplot(data4plotting) +
                          scale_x_continuous(limits = c(input$region[1], input$region[2]))
      if (input$traits == 4){
        expression_plot <- expression_plot + 
                            geom_point(aes(genetic.position, expression.value, color = factor(allele), shape = factor(trait)), size = 3, stat="identity", position = "jitter") +
                            scale_shape_discrete(name = "Trait")
        } else{
        expression_plot <- expression_plot + 
                            geom_point(aes(genetic.position, expression.value, color = factor(allele)), size = 3, stat="identity", position = "jitter")
      }
      expression_plot +
        xlab("Genetic Position in cM") +
        ylab("Relative Gene Expression") +
        scale_color_discrete(name = "Allele") +
        theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    }
    # plots physical distance
    else if(input$distance == 2){
      data4plotting <- subset(data4plotting, genetic.position >= input$region[1] & genetic.position <= input$region[2])
      if (input$traits == 4){
        expression_plot <- ggplot(data4plotting) +
                            geom_point(aes(physical.position, expression.value, color = factor(allele), shape = factor(trait)), size = 3, stat="identity", position = "jitter") +
                            scale_shape_discrete(name = "Trait")
      } else{
        expression_plot <- ggplot(data4plotting) +
                            geom_point(aes(physical.position, expression.value, color = factor(allele)), size = 3, stat="identity", position = "jitter")
        }
      expression_plot +
        xlab("Physical Position in Base Pairs") +
        ylab("Relative Gene Expression") +
        scale_color_discrete(name = "Allele") +
        theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    }
  })
  
  # generates the dataset for users to download - complete list of genes in the region that they are viewing
  download_data <- reactive({
    if (input$chromosome == 0){
      if (input$traits == 4){
        data <- subset(qtl, trait == "trait1" | trait == "trait2")
      } else{
        data <- subset(qtl, trait == input$traits)
      }
    } else{
      data_in_region <- subset(qtl, genetic.position >= input$region[1] & genetic.position <= input$region[2])
      if (input$traits == 4){
        data <- subset(data_in_region, chr == input$chromosome & (trait == "trait1" | trait == "trait2"))
      } else{
        data <- subset(data_in_region, chr == input$chromosome & trait == input$traits)
      }
    }
  })
  
  # allows user to download the full gene table
  output$download_table <- downloadHandler(
    filename = function() { paste("chr", input$chromosome, "trait", input$traits, ".csv", sep="") },
    content = function(file) {
      write.csv(download_data(), file)
    }
  )
  
  # generates the data for the table. outputs gene name, relative expressions for both alleles, genetic position, physical position, trait, and chromosome
  # or just gene name, position, and expression?
  # only shows top 10 upregulated and top 10 downregulated genes
  table_data <- reactive({
    selected.data <- subset(qtl, select = c(gene_name, allele, expression.value, genetic.position, physical.position, trait, chr))
    if (input$chromosome == 0){
      if (input$traits == 4){
        data <- subset(selected.data, trait == "trait1" | trait == "trait2")
      } else{
        data <- subset(selected.data, trait == input$traits)
      }
    } else{
      data_in_region <- subset(selected.data, genetic.position >= input$region[1] & genetic.position <= input$region[2])
      if (input$traits == 4){
        data <- subset(data_in_region, chr == input$chromosome & (trait == "trait1" | trait == "trait2"))
      } else{
        data <- subset(data_in_region, chr == input$chromosome & trait == input$traits)
      }
    }
    sorted_data <- data[order(data$expression.value),]
    sorted_data$physical.position = as.character(round(sorted_data$physical.position, digits = 0))
    sorted_data$chrom = as.character(sorted_data$chrom)
    top10downregulated <- head(sorted_data, n=10)
    top10upregulated <- tail(sorted_data, n=10)
    final_table <- rbind(top10downregulated, top10upregulated)
    unique(final_table) # ensures there are no duplicates (in case the gene list < 10)
  })
  
  # shows the data previously retreived in a table
#   output$table <- renderTable({
#     table_data()
#     }, include.rownames = FALSE)
    
})
