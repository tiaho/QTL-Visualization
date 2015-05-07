# combine_traits.R
# combines the 5 traits files into 1

library(reshape)

# reads in the data
biomass <- read.csv("data/biomass.csv")
height <- read.csv("data/height.csv")
leaf_length <- read.csv("data/leaf_length.csv")
period <- read.csv("data/period18C.csv")
leaf_area <- read.csv("data/specific_leaf_area.csv")

# combines it all into a data frame
data <- biomass[1:2]
data$marker <- row.names(biomass)
data$biomass <- biomass$lod
data$height <- height$lod
data$leaf_length <- leaf_length$lod
data$period <- period$lod
data$leaf_area <- leaf_area$lod

# writes out to a new file
write.table(data, file = "data/real_traits.csv", sep = ",")
