library(readr)
library(dplyr)

# Read in all CSVs in the data/preprocessed directory

which_phase <- "phase_2"

files <- list.files(paste0("data/", which_phase, "/preprocessed"), pattern = "*.csv", full.names = TRUE)

cat(files)

data <- lapply(files, function(x){
  read_csv(x, col_types = "cdcllcddlc") %>% 
    dplyr::filter(electrode %in% c("Cz", "Pz")) %>%
    select(-v_range, -good_segment)
}) %>% 
  bind_rows()

write_csv(data, paste0("data/", which_phase, "/final/eeg.csv"))
