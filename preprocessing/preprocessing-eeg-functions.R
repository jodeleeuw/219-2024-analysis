library(dplyr)
library(tidyr)
library(purrr)
library(ravetools)
library(eegkit)
source('preprocessing/extract-events.R')

preprocess_eeg <- function(file, beh.file, subject_id, sampling_rate=500, filter_low=0.1, filter_high=70, notch_filter_low=59, notch_filter_high=61, segment_begin=-100, segment_end=1000, segment_offset=0, bad_segment_range=200, which_electrodes=c('Cz', 'Fz', 'Fp1', 'Fp2')){
  data <- read_eeg_tidy(eeg.file, beh.file, which_electrodes) %>%
    #linked_ears_rereference() %>%
    bandpass(low=filter_low, high=filter_high, sampling_rate=sampling_rate) %>%
    notch(low=notch_filter_low, high=notch_filter_high, sampling_rate=sampling_rate) %>%
    segment(start = segment_begin, end = segment_end, offset=segment_offset, sampling_rate=sampling_rate) %>%
    artifact_rejection(max_range=bad_segment_range) %>%
    baseline_correct() %>%
    mutate(subject = subject_id)
  
  return(data)
}

read_eeg_tidy <- function(eeg.file, beh.file, which_electrodes=NA) {
  head <- edfReader::readEdfHeader(eeg.file)
  signals <- edfReader::readEdfSignals(head, signals="Ordinary")
  
  data <- map_df(signals, "signal")
  data$sample_id <- 1:nrow(data)
  
  data$`Packet Counter` <- NULL
  data$ACC21 <- NULL
  data$ACC22 <- NULL
  data$ACC23 <- NULL
  data$`ExG 1` <- NULL
  data$TRIGGER <- NULL
  
  tidy_data <- pivot_longer(data, 1:19, names_to = "electrode", values_to = "v") %>% 
    select(sample_id, electrode, v, A2) 
  
  if(!any(is.na(which_electrodes))){
    tidy_data <- tidy_data %>% dplyr::filter(electrode %in% which_electrodes)
  }
  
  events <- extract_events(eeg.file, beh.file)
  
  return(list(signals=tidy_data, events=events))
}

linked_ears_rereference <- function(data){
  df <- data$signals
  df <- df %>% 
    mutate(v = v - (A2/2)) %>%
    select(-A2)
  data$signals <- df
  return(data)
}

bandpass <- function(data, low, high, sampling_rate){
  df <- data$signals
  df <- df %>% 
    group_by(electrode) %>% 
    mutate(v = eegfilter(v, sampling_rate, low, high)) %>%
    ungroup()
  data$signals <- df
  return(data)
}


notch <- function(data, low, high, sampling_rate){
  df <- data$signals
  df <- df %>% 
    group_by(electrode) %>% 
    mutate(v = notch_filter(v, sampling_rate, low, high)) %>%
    ungroup()
  data$signals <- df
  return(data)
}

segment <- function(data, start, end, offset, sampling_rate){
  events <- data$events
  signals <- data$signals
  
  ms_per_sample <- 1000/sampling_rate
  min_sample <- start/ms_per_sample + round(offset/ms_per_sample)
  max_sample <- end/ms_per_sample + round(offset/ms_per_sample)
  
  epochs <- events %>% 
    group_by(event_id) %>%
    reframe(
      t=seq(start, end, ms_per_sample), 
      sample_id=seq(sample_id+min_sample, sample_id+max_sample, 1),
      word_type = word_type,
      is_word = is_word,
      correct = correct
    ) %>%
    left_join(signals, by="sample_id", multiple="all") %>%
    select(-sample_id)
  
  return(epochs)    
}

eyeblink_detector_step <- function(v, window_size=200, sampling_rate=500){
  window <- window_size * (sampling_rate/1000)
  max.diff <- -Inf
  for(i in 1:(length(v)-window)){
    diff <- -mean(v[i:(i+(window/2))]) + mean(v[(i+(window/2)+1):(i+window)])
    diff <- abs(diff)
    if(diff > max.diff){
      max.diff <- diff
    }
  }
  return(max.diff)
}

artifact_rejection <- function(epochs, max_range=200, eye_threshold=70){
  d <- epochs %>% 
    group_by(event_id, electrode) %>%
    summarize(v_range = max(v) - min(v), eye_range = eyeblink_detector_step(v))
  
  eyeblinks <- d %>% 
    group_by(event_id) %>%
    summarize(eyeblink = any(eye_range > eye_threshold))
  
  good.segments <- d %>%
    left_join(eyeblinks, by=c("event_id")) %>%
    group_by(event_id, electrode) %>%
    mutate(good_segment = v_range <= max_range & !eyeblink) %>%
    ungroup()
  
  
  epochs.ar <- epochs %>%
    left_join(good.segments, by=c("event_id", "electrode"))
  
  return(epochs.ar)
}

baseline_correct <- function(epochs){
  baseline.means <- epochs %>%
    group_by(electrode, event_id) %>%
    dplyr::filter(t <= 0) %>%
    summarize(baseline.mean = mean(v))
  
  epoch.bc <- epochs %>%
    left_join(baseline.means, by=c("electrode", "event_id")) %>%
    mutate(v = v - baseline.mean) %>%
    select(-baseline.mean)
  
  return(epoch.bc)
}