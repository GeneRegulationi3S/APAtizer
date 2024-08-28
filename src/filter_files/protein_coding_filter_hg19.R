required_packages <- c(
  "data.table", "readr"
)

install_missing_packages <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }
}

install_missing_packages(required_packages)

library("data.table")
library("readr")

# set working directory
setwd("TRIMMED_htseq/")
dir.create("FILTERED")
# create a list of all txt files in the directory
txt_files <- list.files(pattern = "\\.txt",full.names = TRUE)
data <- lapply(txt_files, function(x)read.table(x, header = FALSE, sep = "\t", dec = "."))
counts <- read_delim("../src/filter_files/filter_file_hg19.txt", delim = "\t", 
                     escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
# set path for filtered files
filtered_path <- "FILTERED/"

# loop through each file in data and filter by first column of counts, then save
for (i in seq_along(data)) {
  # extract the original file name (without the extension)
  file_name <- gsub("\\.txt", "", basename(txt_files[i]))
  
  # filter the data frame by first column of counts
  filtered_df <- data[[i]][data[[i]]$V1 %in% counts$X1,]
  
  # write the filtered data frame to a file with the original file name
  write.table(filtered_df, file.path(filtered_path, paste0(file_name, ".txt")), 
              row.names = FALSE,col.names = FALSE, sep = "\t", dec = ".")
}
