# script to merge a bunch of .counts files, which are white-space separated files with two columns: the first containing the gene names, the second containing the counts for each of those genes
# saves the resulting dataframe as `counts.csv`
# run like `Rscript merge_counts.R /path/to/counts`

if(!"magrittr" %in% installed.packages())install.packages("magrittr",repos="http://cran.us.r-project.org")

require(magrittr)

# print session info for reference / debugging purposes
sessionInfo()

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
# change the working directory to whatever the directory passed as the arugment is
setwd(args[1])

# print working directory
getwd()

# check to make sure `counts.csv` isn't in the directory already:
if(file.exists("counts.csv")){
  stop("there's already a counts.csv in the directory. Either rename it or delete it then try again.")
}

# find all files that end in `.counts`:
files <- list.files(pattern = "*\\.counts$")

print(files)

# check to make sure you have more than 1 file:
if(length(files) < 2){
  stop("fewer than two files match the criteria. Make sure the path is correct.")
}

# start by iniating with the first dataframe
print(files[1])
a <- files[1] %>% strsplit(.,".",fixed=T) %>% lapply(., `[[`, 1) %>% unlist
x <- read.table(files[1], col.names=c(a,"gene"))

print(head(x))

# in each iteration of the loop, add one additional column to the dataframe
for(i in 2:length(files)){
  
  print(files[i])
  print(i)
  
  # read in the next data frame
  y <- read.table(files[i], col.names=c("counts","gene"))
  # merge
  n <- files[(i-1):i] %>% strsplit(.,".",fixed=T) %>% lapply(., `[[`, 1) %>% unlist
  x <- merge(x,y,all=T, by="gene", suffixes=n)
  
}

# get rid of `counts`  prefix on column names
names(x) %<>% gsub("counts","",.)

# check to make sure all colum names are unique
if(names(x) %>% unique %>% length != ncol(x)){
  stop("problem assigned column unique names")
}

# when merging with `all = T`, R substitutes NAs instead of 0.
# replace NAs with zeros
x[is.na(x)] <- 0

write.csv(x, file = "counts.csv", row.names=F)

cat("merged dataframe saved as counts.csv.")
