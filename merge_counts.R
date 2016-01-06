# script to merge a bunch of .counts files, which are white-space separated files with two columns: the first containing the gene names, the second containing the counts for each of those genes
# saves the resulting dataframe as `counts.csv`
# run like `Rscript merge_counts.R /path/to/counts`

# print session info for reference / debugging purposes
sessionInfo()

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
# change the working directory to whatever the directory passed as the arugment is
setwd(args[1])

# find all files that end in *counts:
files <- list.files(path = ".", pattern = "*.counts")

# check to make sure you have more than 1 file:
if(length(files) < 2){
  stop("fewer than two files match the criteria. Make sure the path is correct.")
}

# start by iniating with the first dataframe
x <- read.table(files[1], col.names=c("counts","gene"))

# in each iteration of the loop, add one additional column to the dataframe
for(i in 2:length(files)){
  
  # read in the next data frame
  y <- read.table(files[i], col.names=c("counts","gene"))
  # merge
  x <- merge(x,y,all=T, by="gene")
  
}

# change the column names to reflect what individual the reads came from:
names(x)[2:ncol(x)] <- files %>% strsplit(.,".",fixed=T) %>% lapply(., `[[`, 1) %>% unlist

# when merging with `all = T`, R substitutes NAs instead of 0.
# replace NAs with zeros
x[is.na(x)] <- 0

write.csv(x, file = "counts.csv")

cat("merged dataframe saved as counts.csv.")
