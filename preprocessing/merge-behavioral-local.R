library(jsonlite)
library(dplyr)
library(readr)
library(tidyr)

behavior.files <- list.files('data/raw/beh', pattern="json", full.names = T)

merged.data <- lapply(behavior.files, fromJSON) %>% bind_rows()

response_column <- merged.data$response

null_positions <- which(sapply(response_column, is.null))

new_values <- rep(NA, length(response_column))

new_values[-null_positions] <- unlist(response_column, use.names = F)

merged.data$response <- new_values

write_csv(merged.data, file="data/final/behavioral.csv")
