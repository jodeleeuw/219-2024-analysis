library(osfr)
library(dplyr)

get_phase_1 <- FALSE
get_phase_2 <- TRUE

if(get_phase_1){
osf_retrieve_node("v7mbj") %>%
  osf_ls_files(path="Behavioral", n_max=Inf) %>%
  osf_download(path="data/phase_1/raw/beh", conflicts="skip")

osf_retrieve_node("v7mbj") %>%
  osf_ls_files(path="EEG", n_max=Inf) %>%
  osf_download(path="data/phase_1/raw/eeg", conflicts="skip")
}

# phase 2
if(get_phase_2){
osf_retrieve_node("gd24h") %>%
  osf_ls_files(path="Behavioral", n_max=Inf) %>%
  osf_download(path="data/phase_2/raw/beh", conflicts="skip")

osf_retrieve_node("gd24h") %>%
  osf_ls_files(path="EEG", n_max=Inf) %>%
  osf_download(path="data/phase_2/raw/eeg", conflicts="skip")
}
