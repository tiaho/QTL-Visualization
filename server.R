# server.R

library(ggplot2)
library(reshape2)
library(grid)
library(scales)

# reads in data
allele_data <- read.csv("data/allele_specific_test_p_adjusted.csv")
eqtl_data <- read.csv("data/transcripts_eqtl_start_stop_eqtl.csv")
phenotype_data <- read.csv("data/real_traits.csv")

# calculates the fold change for the 2 parent alleles
allele_data$fold_change <- allele_data$R500 - allele_data$IMB211

# extracts the physical location from the phenotypic data
phenotype_data$phy_pos <- sapply(phenotype_data$marker, function (x){strsplit(as.character(x), "x")[[1]][2] })
phenotype_data$phy_pos <- as.numeric(phenotype_data$phy_pos)

# converts the chromosome vaue from Axx to x
phenotype_data$chr <- sapply(as.character(phenotype_data$chr), function(x){ strsplit(x, "A")[[1]][2] })
phenotype_data$chr <- as.numeric(phenotype_data$chr)

# assigns background colors to the phenotype data, based on chromosome
phenotype_data$background.color <- "1"
phenotype_data$background.color[(phenotype_data$chr %% 2) == 0] <- "0"

# converts the physical position from bases to megabases
phenotype_data$phy_pos <- round(phenotype_data$phy_pos / 1000000, digits = 3)

# merges some data frames together to get the desired data frame for plotting the expression graph
expression_data <- merge(x = allele_data, y = eqtl_data, by.x = "gene_name", by.y = "tx_name")

# converts the chromosome vaue from Axx to x
expression_data$tx_chrom <- sapply(as.character(expression_data$tx_chrom), function(x){ strsplit(x, "A")[[1]][2] })
expression_data$tx_chrom <- as.numeric(expression_data$tx_chrom)

# converts the physical position from bases to megabases
expression_data$tx_start <- round(expression_data$tx_start / 1000000, digits = 3)
expression_data$tx_end <- round(expression_data$tx_end / 1000000, digits = 3)


shinyServer(function(input, output, session) {
  
  # slider input
  output$slider <- renderUI({
    qtlChr <- subset(phenotype_data, chr == input$chromosome)
    max_phy_pos <- max(qtlChr$phy_pos, 1)
    sliderInput("region", label = h5("Display a region of the chromosome?"),
     min = 0, max = max_phy_pos, value = c(0, max_phy_pos), step = 0.001)
  })
  
  
  # qtl graph
  output$qtl_graph <- renderPlot({
    
    # determines the trait that is inputted
    trait_num = input$traits[1]
#     trait_pos <- grep(trait_num, names(phenotype_data))
    trait_pos <- as.numeric(trait_num) + 3
    lod <- phenotype_data[, trait_pos]
    trait <- rep(names(phenotype_data)[trait_pos], length(lod))
    
    # extracts the needed info from phenotype_data
    qtlData <- subset(phenotype_data, select = c(chr, pos, phy_pos, background.color))
    qtlData <- cbind(qtlData, lod, trait)
    
    # if there's a second trait selected...
    if (length(input$traits) > 1){
      # determines the trait that is inputted
      trait_num2 = input$traits[2]
      trait_pos2 <- grep(trait_num2, names(phenotype_data))
      lod <- phenotype_data[, trait_pos2]
      trait <- rep(trait_num2, length(lod))

      # extracts the needed info from phenotype_data
      qtlData2 <- subset(phenotype_data, select = c(chr, pos, phy_pos, background.color))
      qtlData2 <- cbind(qtlData2, lod, trait)

      qtlData <- rbind(qtlData, qtlData2)
    }
    
    # if there's a third trait selected...
    if (length(input$traits) > 2){
      # determines the trait that is inputted
      trait_num3 = input$traits[3]
      trait_pos3 <- grep(trait_num3, names(phenotype_data))
      lod <- phenotype_data[, trait_pos3]
      trait <- rep(trait_num3, length(lod))
      
      # extracts the needed info from phenotype_data
      qtlData3 <- subset(phenotype_data, select = c(chr, pos, phy_pos, background.color))
      qtlData3 <- cbind(qtlData3, lod, trait)
      
      qtlData <- rbind(qtlData, qtlData3)
    }
    
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
      geom_line(aes(x = phy_pos, y = lod, color = as.factor(trait)), size = 2) +
      geom_hline(yintercept = 3, color = "red", size = 1) +
      geom_segment(aes(x = phy_pos, xend = phy_pos), y = (peak * -0.02), yend = (peak * -0.05)) +
      scale_y_continuous(expand = c(0, 0), limits = c((peak * -0.06), max (5, peak))) +
      ggtitle("LOD Curves for QTLs") +
      xlab("Position in Megabases") +
      ylab("LOD Score") +
      theme(axis.text = element_text(size=12),
            axis.text.x = element_text(angle = 90),
            axis.title = element_text(size=16),
            axis.title.y = element_text(vjust = 2),
            axis.line=element_line(),
            title = element_text(size=16),
            legend.position = "top",
            plot.margin = unit(c(0, 0, 0, 0.55), "cm"),
            panel.margin = unit(0, "cm"))
    
    if (length(input$traits) == 1){
      qtl_plot <- qtl_plot +
                    scale_colour_manual(name = "Trait", values = "black")
    } else if (length(input$traits) == 2){
      qtl_plot <- qtl_plot +
                    scale_colour_manual(name = "Trait", values =c("black", "blue"))
    } else if (length(input$traits) == 3){
      qtl_plot <- qtl_plot +
                    scale_colour_manual(name = "Trait", values =c("black","blue", "green"))
    }

    qtl_plot <- qtl_plot + guides(fill=FALSE)

    # plots the graph
    qtl_plot
    
  })
  
  # expression graph
  output$expression_graph <- renderPlot({
    # subsets the data depending on chromosome selected
    data4plotting <- subset(expression_data, tx_chrom == input$chromosome & tx_start >= input$region[1] & tx_start <= input$region[2],
                            select = c("tx_start", "t_stat", "fold_change"))
    t_limits <- c(-max(abs(data4plotting$t_stat),na.rm=T),max(abs(data4plotting$t_stat),na.rm=T))
    fc_limits <- c(-max(abs(data4plotting$fold_change),na.rm=T),max(abs(data4plotting$fold_change),na.rm=T))

    if (input$ex_graph == 1){ # t-statistic
      expression_plot <- ggplot(data4plotting) +
                          geom_point(aes(tx_start, t_stat, color = fold_change)) +
                          scale_colour_gradientn(colours = c("red", "red1", "red2", "red3", "black", "blue3", "blue2", "blue1", "blue"), name = "log2(fold change)", limits=fc_limits) +
                          ylab("t-statistic")
    } else { # fold change
      expression_plot <- ggplot(data4plotting) +
                          geom_point(aes(tx_start, fold_change, color = t_stat)) +
                          scale_colour_gradientn(colours = c("red", "red1", "red2", "red3", "red4", "black", "blue4", "blue3", "blue2", "blue1", "blue"), name = "t-statistic",limits=t_limits) +
                          ylab("log2(fold change)")
                          
                          
    }
    expression_plot +
      xlab("Position in Megabases") +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=16),
            title = element_text(size=16),
            legend.position = "top")
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
    filename = function() {
      if (input$chromosome == 0){
        c("all_chr.csv")
      } else {
        paste("chr", input$chromosome, "_pos", input$region[1], ":", input$region[2], ".csv", sep="")
      }
    },
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
