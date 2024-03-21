library(jsonlite)
library(dplyr)
library(readr)
library(tidyr)

which_phase <- "phase_2"

behavior.files <- list.files(paste0('data/', which_phase, '/raw/beh'), pattern="json", full.names = T)

merged.data <- lapply(behavior.files, fromJSON) %>% bind_rows()

response_column <- merged.data$response

null_positions <- which(sapply(response_column, is.null))

new_values <- rep(NA, length(response_column))

new_values[-null_positions] <- unlist(response_column, use.names = F)

merged.data$response <- new_values

write_csv(merged.data, file=paste0("data/", which_phase, "/final/behavioral.csv"))
