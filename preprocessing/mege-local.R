library(readr)
library(dplyr)

# Read in all CSVs in the data/preprocessed directory

files <- list.files("data/preprocessed", pattern = "*.csv", full.names = TRUE)
data <- lapply(files, function(x){
  read_csv(x, col_types = "cdcllcddddllc")
}) %>% 
  bind_rows() %>% 
  dplyr::filter(electrode %in% c("Cz", "Pz"))

write_csv(data, "data/final/eeg.csv")
