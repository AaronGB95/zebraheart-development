##-----------------------------------------------------------------------------
##
## runall.R
##
## script that executes the whole project
##
##-----------------------------------------------------------------------------

dir_docs <- "01-documentation/"
dir_src <- "02-scripts/"
dir_data <- "03-data/"
dir_output <- "04-output/"

annotations <- read.table(paste0(dir_docs, "annotations.txt"),
                          sep = "\t",
                          header = TRUE)