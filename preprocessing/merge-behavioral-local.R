library(jsonlite)
library(dplyr)
library(readr)

behavior.files <- list.files('data/raw/beh', pattern="json", full.names = T)

merged.data <- lapply(behavior.files, fromJSON) %>% bind_rows()

write_csv(merged.data, file="data/final/behavioral.csv")
