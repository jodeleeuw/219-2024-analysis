library(readr)
library(dplyr)

# Read in all CSVs in the data/preprocessed directory

files <- list.files("data/preprocessed", pattern = "*.csv", full.names = TRUE)
data <- lapply(files, function(x){
  read_csv(x, col_types = "cdcllcddlc") %>% 
    dplyr::filter(electrode %in% c("Cz", "Pz")) %>%
    select(-v_range, -good_segment)
}) %>% 
  bind_rows()

write_csv(data, "data/final/eeg.csv")
