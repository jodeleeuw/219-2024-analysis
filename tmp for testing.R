source("preprocessing/preprocessing-eeg-functions.R")

subject <- "01"

eeg.file <- list.files('data/raw/eeg', pattern=paste0("subject-", subject), full.names = T)
beh.file <- paste0('data/raw/beh/219_2024_behavioral_', subject,'.json')

sampling_rate = 500
filter_low = 0.1
filter_high = 30
segment_begin = -200
segment_end = 1000
segment_offset = 0
bad_segment_range = 500
eye_threshold = 70
which_electrodes = c("Fp1", "Fp2", "Cz", "Pz")
eye_channels=c("Fp1", "Fp2")
target_channels=c("Cz", "Pz")

epochs <- read_eeg_tidy(eeg.file, beh.file, which_electrodes) %>%
  #linked_ears_rereference() %>%
  bandpass(low=filter_low, high=filter_high, sampling_rate=sampling_rate) %>%
  #notch(low=notch_filter_low, high=notch_filter_high, sampling_rate=sampling_rate) %>%
  segment(start = segment_begin, end = segment_end, offset=segment_offset, sampling_rate=sampling_rate)


