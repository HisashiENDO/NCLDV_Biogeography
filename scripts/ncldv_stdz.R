######################################################################
# This script aims to normalize the abundance profiles of NCLDV PolBs 
# 1. Load raw data and cobvert it to the data frame
# 2. Standardize raw data 
#######################################################################

# Import library
library(tidyverse)

################## Prepare data frame ##################
# Load abundance profile
ncldv_raw <- read_tsv("PolB_NCLDV_table.txt", col_names = FALSE) #including sample lavel

# Cobvert fo dataframe for the downstream analysis
ncldv_raw_rev <- as.data.frame(t(ncldv_raw))
write_tsv(ncldv_raw_rev, "ncldv_raw_rev.txt", col_names = FALSE)
ncldv_df <- read_tsv("ncldv_raw_rev.txt", col_names = TRUE) #including sample lavel
class(ncldv_df)
#############################################################
#
#

################## Data standardization ##################
# Remove samples with the coverage lower than 50
ncldv_select_df <- ncldv_df %>% filter(rowSums(ncldv_df[,2:6819]) > 50)

# Remove id column
ncldv_select_freq <- ncldv_select_df[,-1]
rownames(ncldv_select_freq) <- ncldv_select_df$ids
  
# Standardized by the sample with the lowest sum of length-normalized PolB abundance value. 
data_sel_row <- rowSums(ncldv_select_freq)
sum_min <- min(data_sel_row)
row_min <- which.min(data_sel_row)
ncldv_select_freq[ncldv_select_freq[,] == 0] <- NA 
factor_min <- min(ncldv_select_freq[row_min,], na.rm = TRUE) 
ncldv_select_freq[is.na(ncldv_select_freq[,])] <- 0 
ncldv_select_freq2 <- as.data.frame(ncldv_select_freq)
start <- 1
end <- as.integer(length(ncldv_select_freq2[,1])) 
data_stdz <- NULL 

for (i in start:end){
  stdz <- ncldv_select_freq[i,]*(sum_min/rowSums(ncldv_select_freq[i,]))
  stdz[stdz[1,] < factor_min] <- 0
  data_stdz <- rbind(data_stdz, stdz)
}

#C heck if processing completed successfully
data_stdz_row <- rowSums(data_stdz)
min(data_stdz_row) == sum_min  # FALSE should be returned.
data_stdz_col <- colSums(data_stdz)
barplot(data_stdz_row, ylim=c(0,120))
barplot(data_stdz_col)
data_stdz[data_stdz[,] == 0] <- NA
min(data_stdz[,], na.rm=TRUE) == factor_min # TRUE should be returned.
data_stdz[is.na(data_stdz[,])] <- 0

# Standardize rowSums to 100 for each sample
start <- 1 #forループの開始値
data_stdz <- as.data.frame(data_stdz)
end2 <- as.integer(length(data_stdz[,1]))
data_stdz2 <- NULL
for (i in start:end2){
  stdz <- data_stdz[i,]*(100/rowSums(data_stdz[i,])) 
  data_stdz2 <- rbind(data_stdz2, stdz) 
}

# Check if empty colum exists
if (min(colSums(data_stdz2)) != 0){
  print("There is no blank gene")
} else {
  print("There are gene colum(s) which needs to be removed")
}

### Important!!!! ####
# Remove line of "TARA_155_SRF_lt-0.22_G" because it has only 1 gene (TARA_155_SRF_lt-0.22_G_scaffold3310_2_gene25079)
data_stdz3 <- data_stdz2[-262,]

# Save standardized df 
write_tsv(data_stdz3, "ncldv_stdz_df.txt", col_names = TRUE)
#############################################################
#
#

