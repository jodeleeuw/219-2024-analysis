library(readr)
library(stringr)
source("preprocessing/preprocessing-eeg-functions.R")

subjects <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14")

for(subject in subjects){
  
  final.path <- paste0("data/preprocessed/subject-", subject,"-epochs.csv")
  print(paste0("Working on ",subject))
  if(file.exists(final.path)){
    print(paste0("epochs.csv already exists - skipping"))
    next
  }
  tryCatch({
    eeg.file <- list.files('data/raw/eeg', pattern=paste0("subject-", subject), full.names = T)
    beh.file <- paste0('data/raw/beh/219_2024_behavioral_', subject,'.json')
    
    preprocessed.data <- preprocess_eeg(
      file = eeg.file, 
      beh.file = beh.file,
      subject_id = subject,
      sampling_rate = 500,
      filter_low = 0.1,
      filter_high = 70,
      notch_filter_low = 59,
      notch_filter_high = 61,
      segment_begin = -200,
      segment_end = 1000,
      segment_offset = 0,
      bad_segment_range = 200,
      which_electrodes = c("Fp1", "Fp2", "Cz", "Pz"))
    
    write_csv(preprocessed.data, file=paste0("data/preprocessed/subject-", subject,"-epochs.csv"))
  },
  error = function(cond){
    message(cond)
  },
  warning = function(cond){
    message(cond)
  })
}
