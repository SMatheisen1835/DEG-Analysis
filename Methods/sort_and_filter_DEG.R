# This script is used to filter out and keep only DEGs that have an
# adj.p.val <=0.05 and |logFC| >=0.58

# Versions: R = 4.2.2, dplyr = , openxlsx = 4.2.7.1

###############################################################################

# setting the working directory to make sure the paths work as intended 
library(rstudioapi)
script_dir <- dirname(getActiveDocumentContext()$path)
setwd(paste0(script_dir, "/../"))
getwd()

# this is the path from the wd to the GSE folder 
path_to_table_folder <- "Material/GSE137268/"
# these are the paths from the GSE-folder to the tables
table_paths <- c("Case=Controlled Asthma/GSE137268.CA.tsv", 
                 "Case=Severe Asthma/GSE137268.SA.tsv", 
                 "Case=Uncontrolled Asthma/GSE137268.UA.tsv")

for (table_path in table_paths){
  # Read the table
  data <- read.table(paste0(path_to_table_folder, table_path), header = TRUE, sep = "\t")
  
  # filter and keep only if row has  adj.P.val < 0.05 and |logFC|>0.58 
  library(dplyr)
  filtered_data <- data %>%
    filter(adj.P.Val <= 0.05) %>%
    filter(abs(logFC) >= 0.58)
  
  # write the resulting  table in an excel file
  library(openxlsx)
  
  # extracts the GSE number from the path
  gse = substr(path_to_table_folder, 10, nchar(path_to_table_folder))
  # extracts the orignal file name but cuts off the .tsv
  file_name = substr(table_path, regexpr("/",table_path), nchar(table_path)-3)
  # pastes together the path and the file name and adds "_filtered.xlsx"
  table_path_and_name = paste0("Results/", gse, file_name, "_filtered.xlsx")
  
  print(paste0("Filtered table will be created as: ", table_path_and_name))
  
  write.xlsx(filtered_data, table_path_and_name, rowNames = FALSE)  
}
