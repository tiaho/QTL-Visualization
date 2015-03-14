# server.R

library(ggplot2)
library(reshape2)
library(grid)

allele_data <- read.csv("data/allele_specific_test_p_adjusted.csv")
eqtl_data <- read.csv("data/transcripts_eqtl_start_stop_eqtl.csv")
phenotype_data <- read.csv("data/simulated_phenotypes_real_map.csv")

# calculates the fold change for the 2 parent alleles
allele_data$fold_change <- allele_data$IMB211/allele_data$R500

# extracts the physical location from the phenotypic data
for (i in 1:length(phenotype_data$marker)){
   physical_position <- strsplit(as.character(phenotype_data$marker[i]), "x")[[1]][2]
   phenotype_data$phy_pos[i] <- as.numeric(physical_position)
}


# converts the chromosome vaue from Axx to x
phenotype_data$chr <- as.character(phenotype_data$chr)
for (i in 1:length(phenotype_data$chr)){
  if (phenotype_data$chr[i] == "A01"){
    phenotype_data$chr[i] = 1
  } else if (phenotype_data$chr[i] == "A02"){
    phenotype_data$chr[i] = 2
  } else if (phenotype_data$chr[i] == "A03"){
    phenotype_data$chr[i] = 3
  } else if (phenotype_data$chr[i] == "A04"){
    phenotype_data$chr[i] = 4
  } else if (phenotype_data$chr[i] == "A05"){
    phenotype_data$chr[i] = 5
  } else if (phenotype_data$chr[i] == "A06"){
    phenotype_data$chr[i] = 6
  } else if (phenotype_data$chr[i] == "A07"){
    phenotype_data$chr[i] = 7
  } else if (phenotype_data$chr[i] == "A08"){
    phenotype_data$chr[i] = 8
  } else if (phenotype_data$chr[i] == "A09"){
    phenotype_data$chr[i] = 9
  } else if (phenotype_data$chr[i] == "A10"){
    phenotype_data$chr[i] = 10
  }
}
phenotype_data$chr <- as.numeric(phenotype_data$chr)

# assigns background colors to the phenotype data, based on chromosome
phenotype_data$background.color <- "1"
phenotype_data$background.color[(phenotype_data$chr %% 2) == 0] <- "0"

# merges some data frames together to get the desired data frame for plotting the expression graph
expression_data <- merge(x = allele_data, y = eqtl_data, by.x = "gene_name", by.y = "tx_name")

# converts the chromosome vaue from Axx to x
expression_data$tx_chrom <- as.character(expression_data$tx_chrom)
for (i in 1:length(expression_data$tx_chrom)){
  if (expression_data$tx_chrom[i] == "A01"){
    expression_data$tx_chrom[i] = 1
  } else if (expression_data$tx_chrom[i] == "A02"){
    expression_data$tx_chrom[i] = 2
  } else if (expression_data$tx_chrom[i] == "A03"){
    expression_data$tx_chrom[i] = 3
  } else if (expression_data$tx_chrom[i] == "A04"){
    expression_data$tx_chrom[i] = 4
  } else if (expression_data$tx_chrom[i] == "A05"){
    expression_data$tx_chrom[i] = 5
  } else if (expression_data$tx_chrom[i] == "A06"){
    expression_data$tx_chrom[i] = 6
  } else if (expression_data$tx_chrom[i] == "A07"){
    expression_data$tx_chrom[i] = 7
  } else if (expression_data$tx_chrom[i] == "A08"){
    expression_data$tx_chrom[i] = 8
  } else if (expression_data$tx_chrom[i] == "A09"){
    expression_data$tx_chrom[i] = 9
  } else if (expression_data$tx_chrom[i] == "A10"){
    expression_data$tx_chrom[i] = 10
  }
}
expression_data$tx_chrom <- as.numeric(expression_data$tx_chrom)

shinyServer(function(input, output) {
  
  # slider input
  output$slider <- renderUI({
    qtlChr <- subset(phenotype_data, chr == input$chromosome)
    max_pos <- ceiling(max(qtlChr$pos, 1))
    sliderInput("region", label = h5("Display a region of the chromosome?"),
     min = 0, max = max_pos, value = c(0, max_pos), step = 1)
  })

  # qtl graph
  output$qtl_graph <- renderPlot({
    # subsets the data depending on the trait selected
    if (input$traits == "trait1"){
      qtlData <- subset(phenotype_data, select = c(chr, pos, phy_pos, trait1, background.color))
    } else if (input$traits == "trait2"){
      qtlData <- subset(phenotype_data, select = c(chr, pos, phy_pos, trait2, background.color))
#     } else {
#       qtlData <- subset(phenotype_data, select = c(chr, pos, lod, background.color))
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
    } else if (input$chromosome == 5){
      qtlData <- subset(qtlData, chr == 5)
    } else if (input$chromosome == 6){
      qtlData <- subset(qtlData, chr == 6)
    } else if (input$chromosome == 7){
      qtlData <- subset(qtlData, chr == 7)
    } else if (input$chromosome == 8){
      qtlData <- subset(qtlData, chr == 8)
    } else if (input$chromosome == 9){
      qtlData <- subset(qtlData, chr == 9)
    } else if (input$chromosome == 10){
      qtlData <- subset(qtlData, chr == 10)
    }
    
    # changes the column name of traitx to lod
    colnames(qtlData)[4] <- "lod"
    if (input$traits == 4){
      colnames(qtlData)[5] <- "lod2"
    }          
    
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
      geom_hline(yintercept = 4, color = "red", size = 1) +
      geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
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
        geom_line(aes(x = pos, y = lod2), size = 2, color = "blue")
    } else {
      qtl_plot
    }
  })
  
  # expression graph
  output$expression_graph <- renderPlot({
    # subsets the data depending on chromosome selected
    data4plotting <- subset(expression_data, tx_chrom == input$chromosome, select = c("tx_start", "t_stat", "fold_change"))

    
    # plots physical distance
#       data4plotting <- subset(data4plotting, genetic.position >= input$region[1] & genetic.position <= input$region[2])
    if (input$ex_graph == 1){ # t-statistic
      expression_plot <- ggplot(data4plotting) +
                          geom_point(aes(tx_start, t_stat)) +
                          xlab("Physical Position in Base Pairs") +
                          ylab("t-statistic") +
                          theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    } else { # fold change
      expression_plot <- ggplot(data4plotting) +
                          geom_point(aes(tx_start, fold_change)) +
                          xlab("Physical Position in Base Pairs") +
                          ylab("Fold Change") +
                          theme(axis.text = element_text(size=12), axis.title = element_text(size=16), title = element_text(size=16))
    }
      expression_plot
  })
  
  # generates the dataset for users to download - complete list of genes in the region that they are viewing
  download_data <- reactive({
    if (input$chromosome == 0){
        data <- expression_data
    } else{
#       data_in_region <- subset(expression_data, genetic.position >= input$region[1] & genetic.position <= input$region[2])
        data <- subset(expression_data, chr == input$chromosome)
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
#       data_in_region <- subset(selected.data, genetic.position >= input$region[1] & genetic.position <= input$region[2])
        data <- subset(selected.data, tx_chrom == input$chromosome)
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
