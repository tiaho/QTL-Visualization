library(ggplot2)
library(reshape)

setwd("~/Desktop/Lab/Maloof Lab/QTL_Visualization/")

# reads in the data and formats them as needed
eqtl_data <- read.csv("data/transcripts_eqtl_start_stop_eqtl.csv")

phenotype_data <- read.csv("data/real_traits.csv")
phenotype_data$chr <- sapply(as.character(phenotype_data$chr), function(x){ strsplit(x, "A")[[1]][2] })
phenotype_data$chr <- as.numeric(phenotype_data$chr)
phenotype_data$phy_pos <- sapply(phenotype_data$marker, function (x){strsplit(as.character(x), "x")[[1]][2] })
phenotype_data$phy_pos <- as.numeric(phenotype_data$phy_pos)

bin_data <- read.csv("~/Desktop/Lab/Maloof Lab/QTL_Visualization/data/gene_marker_ranges.csv")
bin_data$chr <- sapply(as.character(bin_data$chr), function(x){ strsplit(x, "A")[[1]][2] })
bin_data$chr <- as.numeric(bin_data$chr)

# xtest <- c(0, 15, 20)
# ytest <- c("trait1", "trait2", "trait3")
# filltest <- c(1, 2, 3)
# dftest <- data.frame(xtest, ytest, filltest)
# 
# qplot(xtest, ytest, data = dftest, geom="tile", fill = filltest)
# 
# ggplot(dftest, aes(x = xtest, y = ytest)) +
#   geom_tile(aes(fill = filltest))


# merges the eqtl dat and bin data
eb_merge <- merge(bin_data, eqtl_data, by.x = "tx_name", by.y = "tx_name")

# test vectors, just a bit of each df
# testp <- head(phenotype_data, 10)
# testb <- head(bin_data, 50)

# sorts the data frames by chr, then position
testp <- phenotype_data[with(phenotype_data, order(chr, phy_pos)), ]
testb <- bin_data[with(bin_data, order(chr, bin.start)), ]

# matches up corresponding lines in the eqtl-bin data and the phenotype data
count <- 1
tmp_vectors
tmp
for (i in 1:nrow(testb)){
  if (count > nrow(testp)){
    break
  } else if (testp$phy_pos[count] >= testb$bin.start[i] & testp$phy_pos[count] <= testb$bin.end[i]){
    tmp <- cbind(testb[i,], testp[count,])
    if (i == 1){
      tmp_vectors <- tmp
    } else {
      tmp_vectors <- rbind(tmp_vectors, tmp)
    }
  } else {
    count = count + 1
    tmp <- cbind(testb[i,], testp[count,])
    tmp_vectors <- rbind(tmp_vectors, tmp)
  }
}

# makes the result a data frame
new_data <- as.data.frame(tmp_vectors)

# removes the duplicate chromosome column at the beginning
one_chr_data <- new_data[, 2:ncol(new_data)]

# melts the data frames
melted_data <- melt(one_chr_data, id = c("chr", "bin.start", "bin.end",
                                      "tx_name", "pos", "marker",
                                      "phy_pos", "background.color"))
# plot!
ggplot(melted_data, aes(x = bin.start, y = variable)) +
  geom_tile(aes(fill = value)) +
  facet_grid(~ chr, scales = "free_x", space = "free_x")
