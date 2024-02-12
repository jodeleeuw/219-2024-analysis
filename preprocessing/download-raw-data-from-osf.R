library(osfr)
library(dplyr)

osf_retrieve_node("v7mbj") %>%
  osf_ls_files(path="Behavioral", n_max=Inf) %>%
  osf_download(path="data/raw/beh", conflicts="skip")

osf_retrieve_node("v7mbj") %>%
  osf_ls_files(path="EEG", n_max=Inf) %>%
  osf_download(path="data/raw/eeg", conflicts="skip")
