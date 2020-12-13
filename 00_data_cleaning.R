setwd("~/Desktop/FLAX/LILE_common_garden_proj")

data <- read.csv("raw_data/LILE_yield_data_2013_seed_wt.csv", header=T, na.strings = c("","NA")) #use na.strings to fill blank cells with NA

library(dplyr)

# this code is to fill in empty rows of the flax data from USFS. As of 5/29 only the seed_wt data sheet has been cleaned/filled using this code. Empty rows are being filled with the values in the row up to 3 rows above it (data was only entered for the first out of 5 replicates per plot, this is to fill in the source, block, and field row numbers for each of the other 4 replicates, which share the same values). It will also delete rows in which all elements are NA. It could be easily modified to fill empty rows of other data sets

# basic test
test <- c(14,NA,NA,NA,NA,12,NA,NA,NA,NA)
for(i in 1:length(test)){
  if(is.na(test[i]) & is.na(test[i+1])){
    test[i] <- test[i-1]
  }
}

# test on a single column
head(data$source)
for(i in 1:length(data$source)){
  #print(data$source[i])
  if(is.na(data$source[i]) & is.na(data$source[i+1])){
    data$source[i] <- data$source[i-1]
  }
}

# run on the seed weight data
for(i in 1:(ncol(data)-3)){ #skip the last three columns because they don't need to be filled
  for(j in 1:length(data[,i])){ #iterate over elements in ith column
    if(is.na(data[,i][j]) & is.na(data[,i][j+1])){ #if an element is na and the next element is also na, fill in this element with the preceeding value
      data[,i][j] <- data[,i][j-1]
    }
  }
}

# get rid of spacer rows, which had all NAs in each cell#
data <- data[rowSums(is.na(data)) != ncol(data),] 

# write df as csv
write.csv(data, "~/Desktop/FLAX/LILE_common_garden_proj/filled_LILE3_yield_data_2013_seed_wt.csv", row.names = F)

#### cleaning fruit fill data
data <- read.csv("raw_data/LILE_yield_data_2013_fruit_fill.csv", header=T, na.strings = c("","NA"))

# remove completely empty rows
data <- data %>% filter_all(any_vars(!is.na(.)))

# only need to fill in first two columns, 'source' and 'trt'
for(i in 1:2){ #skip the last three columns because they don't need to be filled
  for(j in 1:length(data[,i])){ #iterate over elements in ith column
    if(is.na(data[,i][j])){ #if an element is na and the next element is also na, fill in this element with the preceeding value
      data[,i][j] <- data[,i][j-1]
    }
  }
}
View(data)
write.csv(data, "~/Desktop/FLAX/LILE_common_garden_proj/cleaned_LILE_yield_data_2013_fruit_fill.csv", row.names = F)

#### stem data
stem_data <- read.csv("raw_data/LILE_yield_data_2013_stem_and_fruit.csv", header=T, na.strings = c("","NA"))

head(stem_data)
# remove empty rows
stem_data <- stem_data %>% filter_all(any_vars(!is.na(.)))

# only need to fill in first two columns, 'source' and 'trt'
for(i in 1:6){ #only need to fill the first 6 columns
  for(j in 1:length(stem_data[,i])){ #iterate over elements in ith column
    if(is.na(stem_data[,i][j])){ #if an element is na and the next element is also na, fill in this element with the preceeding value
      stem_data[,i][j] <- stem_data[,i][j-1]
    }
  }
}
View(stem_data)

write.csv(stem_data, "~/Desktop/FLAX/LILE_common_garden_proj/cleaned_LILE_yield_data_2013_stem_and_fruit.csv", row.names = F)
