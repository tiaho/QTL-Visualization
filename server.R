# server.R

library(ggplot2)
library(reshape2)
library(grid)

# reads in data
allele_data <- read.csv("data/allele_specific_test_p_adjusted.csv")
eqtl_data <- read.csv("data/transcripts_eqtl_start_stop_eqtl.csv")
phenotype_data <- read.csv("data/simulated_phenotypes_real_map.csv")

# calculates the fold change for the 2 parent alleles
allele_data$fold_change <- abs(allele_data$IMB211 - allele_data$R500)

# extracts the physical location from the phenotypic data
phenotype_data$phy_pos <- sapply(phenotype_data$marker, function (x){strsplit(as.character(x), "x")[[1]][2] })
phenotype_data$phy_pos <- as.numeric(phenotype_data$phy_pos)

# converts the chromosome vaue from Axx to x
phenotype_data$chr <- sapply(as.character(phenotype_data$chr), function(x){ strsplit(x, "A")[[1]][2] })
phenotype_data$chr <- as.numeric(phenotype_data$chr)

# assigns background colors to the phenotype data, based on chromosome
phenotype_data$background.color <- "1"
phenotype_data$background.color[(phenotype_data$chr %% 2) == 0] <- "0"

# merges some data frames together to get the desired data frame for plotting the expression graph
expression_data <- merge(x = allele_data, y = eqtl_data, by.x = "gene_name", by.y = "tx_name")

# converts the chromosome vaue from Axx to x
expression_data$tx_chrom <- sapply(as.character(expression_data$tx_chrom), function(x){ strsplit(x, "A")[[1]][2] })
expression_data$tx_chrom <- as.numeric(expression_data$tx_chrom)


shinyServer(function(input, output, session) {
  
  # slider input
  output$slider <- renderUI({
    qtlChr <- subset(phenotype_data, chr == input$chromosome)
    max_phy_pos <- ceiling(max(qtlChr$phy_pos, 1))
    sliderInput("region", label = h5("Display a region of the chromosome?"),
     min = 0, max = max_phy_pos, value = c(0, max_phy_pos), step = 1)
  })
  
  
  # qtl graph
  output$qtl_graph <- renderPlot({
    
    # determines the trait that is inputted
    trait_num = input$traits
    trait_pos <- grep(trait_num, names(phenotype_data))
    lod <- phenotype_data[, trait_pos]
    
    # extracts the needed info from phenotype_data
    qtlData <- subset(phenotype_data, select = c(chr, pos, phy_pos, background.color))
    qtlData <- cbind(qtlData, lod)
    
    # subsets the data depending on chromosome selected
    if (input$chromosome == 0){
      # do nothing
    } else {
      qtlData <- subset(qtlData, chr == input$chromosome)
    }     
    
    # extracts the max lod value for the data set
    peak <- max(qtlData$lod)
    
    # start of the plot
    qtl_plot <- ggplot(qtlData)
    
    # determines the color of the background; depends if all chromosomes are plotted, or just one
    # plots a certain region only for single chromosomes (changes xlim)
    if (input$chromosome == 0){
      qtl_plot <- qtl_plot +
        scale_fill_manual(values = c("1" = "gray90", "0" = "white"))
    } else {
      qtl_plot <- qtl_plot +
        scale_fill_manual(values = c("1" = "white", "0" = "white")) +
        scale_x_continuous(limits = c(input$region[1], input$region[2]))
    }
    
    # rest of the plot
    qtl_plot <- qtl_plot +
      facet_grid(~ chr, scales = "free_x", space = "free_x") +
      geom_rect(aes(fill = background.color), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_line(aes(x = phy_pos, y = lod), size = 2) +
      geom_hline(yintercept = 4, color = "red", size = 1) +
      geom_segment(aes(x = phy_pos, xend = phy_pos), y = (peak * -0.02), yend = (peak * -0.05)) +
      scale_y_continuous(expand = c(0, 0), limits = c((peak * -0.06), max (5, peak))) +
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
        geom_line(aes(x = phy_pos, y = lod2), size = 2, color = "blue")
    } else {
      qtl_plot
    }
  })
  
  # expression graph
  output$expression_graph <- renderPlot({
    # subsets the data depending on chromosome selected
    data4plotting <- subset(expression_data, tx_chrom == input$chromosome, select = c("tx_start", "t_stat", "fold_change"))

    
    # plots physical distance
      data4plotting <- subset(data4plotting, tx_start >= input$region[1] & tx_start <= input$region[2])
    if (input$ex_graph == 1){ # t-statistic
      expression_plot <- ggplot(data4plotting) +
                          geom_point(aes(tx_start, t_stat, color = fold_change)) +
                          scale_color_continuous(name = "log2(fold change)") +
                          xlab("Physical Position in Base Pairs") +
                          ylab("t-statistic") +
                          theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    } else { # fold change
      expression_plot <- ggplot(data4plotting) +
                          geom_point(aes(tx_start, fold_change, color = t_stat)) +
                          scale_color_continuous(name = "t-stat") +
                          xlab("Physical Position in Base Pairs") +
                          ylab("log2(fold change)") +
                          theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    }
      expression_plot
  })
  
  # generates the dataset for users to download - complete list of genes in the region that they are viewing
  download_data <- reactive({
    if (input$chromosome == 0){
        data <- expression_data
    } else{
      data_in_region <- subset(expression_data, tx_start >= input$region[1] & tx_start <= input$region[2])
        data <- subset(expression_data, tx_chrom == input$chromosome)
    }
  })
  
  # allows user to download the full gene table
  output$download_table <- downloadHandler(
    filename = function() { paste("chr", input$chromosome, ".csv", sep="") },
    content = function(file) {
      write.csv(download_data(), file)
    }
  )
  
  # generates the data for the table. outputs gene name, relative expressions for both alleles, physical position, chromosome, t-statistic, and fold change
  # or just gene name, position, and expression?
  # only shows top 10 upregulated and top 10 downregulated genes
  table_data <- reactive({
    selected.data <- subset(expression_data, select = c(gene_name, R500, IMB211, t_stat, fold_change, tx_chrom, tx_start))
    if (input$chromosome == 0){
        data <- selected.data
    } else{
      data_in_region <- subset(selected.data, tx_start >= input$region[1] & tx_start <= input$region[2])
        data <- subset(data_in_region, tx_chrom == input$chromosome)
    }
    sorted_data <- data[order(data$fold_change),]
    sorted_data$tx_chrom = as.character(sorted_data$tx_chrom)
    top10downregulated <- head(sorted_data, n=10)
    top10upregulated <- tail(sorted_data, n=10)
    final_table <- rbind(top10downregulated, top10upregulated)
    colnames(final_table)[6] <- "chr"
    colnames(final_table)[7] <- "physical_position"
    unique(final_table) # ensures there are no duplicates (in case the gene list < 10)
  })
  
  # shows the data previously retreived in a table
  output$table <- renderTable({
    table_data()
    }, include.rownames = FALSE)
    
})
