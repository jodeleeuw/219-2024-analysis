library(dplyr)
library(tidyr)
library(purrr)
library(ravetools)
library(eegkit)
source('preprocessing/extract-events.R')

preprocess_eeg <- function(file, beh.file, subject_id, sampling_rate=500, filter_low=0.1, filter_high=70, notch_filter_low=59, notch_filter_high=61, segment_begin=-100, segment_end=1000, segment_offset=0, bad_segment_range=200, eye_threshold=70, which_electrodes=c('Cz', 'Fz', 'Fp1', 'Fp2')){
  data <- read_eeg_tidy(eeg.file, beh.file, which_electrodes) %>%
    #linked_ears_rereference() %>%
    bandpass(low=filter_low, high=filter_high, sampling_rate=sampling_rate) %>%
    #notch(low=notch_filter_low, high=notch_filter_high, sampling_rate=sampling_rate) %>%
    segment(start = segment_begin, end = segment_end, offset=segment_offset, sampling_rate=sampling_rate) %>%
    artifact_rejection(max_range=bad_segment_range, eye_threshold = eye_threshold) %>%
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
    select(sample_id, electrode, v, A2) %>%
    select(-A2) #because we are not doing linked ear re-referencing
  
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
    mutate(v = as.numeric(eegfilter(v, sampling_rate, low, high))) %>% #as numeric because eegfilter returns a list column
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

eyeblink_points <- function(v, blink_criteria=15, sampling_rate=500, window_ms=200){
  window_length <- window_ms * (sampling_rate/1000)
  half_window <- window_length / 2
  blink_cov <- c(rep(+1, half_window), rep(-1, half_window))
  cv <- stats::filter(v, blink_cov, method="convolution", sides=2) / length(blink_cov)
  blink_pts <- abs(cv) >= blink_criteria
  blink_group <- cumsum(c(0, abs(diff(as.integer(blink_pts)))))
  blink_mean <- mean(v[blink_pts])
  
  blink.df <- data.frame(v=v, is_blinking=blink_pts, blink_group=blink_group)
  blink.df <- blink.df %>% 
    group_by(blink_group) %>%
    mutate(blink_group_avg=mean(v)) %>%
    ungroup() %>%
    mutate(is_blinking = is_blinking | blink_group_avg > blink_mean)
  
  return(blink.df %>% pull(is_blinking))
}

eyeblink_removal_with_regression <- function(epochs, eye_channels=c("Fp1", "Fp2"), target_channels=c("Cz", "Pz")){
  
  # compute a single average for the eye channels for each condition
  eyes <- epochs %>% 
    dplyr::filter(electrode %in% eye_channels) %>%
    group_by(event_id, word_type, is_word, t) %>%
    summarize(v = mean(v))
  
  # identify all time points that include an eyeblink
  # we will use this to remove the eyeblink from the target channels
  
  eyeblinks <- eyes %>%
    group_by(event_id) %>%
    mutate(is_blinking = eyeblink_points(v, blink_criteria=15, sampling_rate=500, window_ms=200)) %>%
  
  eyeblinks <- eyes %>% 
    group_by(event_id) %>%
    summarize(eyeblink = any(is_blinking))
  
  eyes <- eyes %>%
    left_join(eyeblinks, by="event_id")
  
  #ggplot(eyes %>% dplyr::filter(event_id==3), aes(x=t, y=v)) +
   # geom_line()
  
  # compute the average for the eye channels and for the target channels for each condition
  # this generate the event-related signal
  
  eye <- epochs %>% 
    dplyr::filter(electrode %in% eye_channels) %>%
    group_by(word_type, is_word, t) %>%
    summarize(v = mean(v))
  
  target <- epochs %>% 
    dplyr::filter(electrode %in% target_channels) %>%
    group_by(electrode, word_type, is_word, t) %>%
    summarize(v = mean(v))
 
  # ggplot(target, aes(x=t, y=v, color=word_type, linetype=is_word)) +
  #   geom_line() +
  #   facet_grid(word_type~electrode)
  
  # subtract the average of the eye channels from each individual
  # trial in the eye channels. trials are defined by event_id.
  # merge to a single "virtual" EOG electrode
  
  eye_single_trials <- epochs %>%
    dplyr::filter(electrode %in% eye_channels) %>%
    group_by(event_id, word_type, is_word, t) %>%
    summarize(v = mean(v)) %>%
    left_join(eye, by=c("word_type", "is_word", "t")) %>%
    mutate(v = v.x - v.y) %>%
    select(-v.x, -v.y)
  
  # do the same for target channels
  
  target_single_trials <- epochs %>% 
    dplyr::filter(electrode %in% target_channels) %>%
    left_join(target, by=c("electrode", "word_type", "is_word", "t")) %>%
    mutate(v = v.x - v.y) %>%
    select(-v.x, -v.y)
  
  # now we can compute the regression coefficients for each trial
  # in the eye channels and the target channels.
  library(broom)
  regression_data <- target_single_trials %>%
    left_join(eye_single_trials, by=c("event_id", "t", "word_type", "is_word")) %>%
    group_by(electrode, word_type, is_word) %>%
    nest() %>%
    mutate(regression = map(data, function(x){
      m <- lm(v.x ~ v.y, data=x)
      return(tidy(m, conf.int=FALSE))
    })) %>%
    unnest(regression) %>%
    select(-data, -std.error, -statistic, -p.value) %>%
    dplyr::filter(term == "v.y") %>%
    select(-term)
  
  # now we can use the coefficients to remove the eyeblink from the target channels
  # for each trial
  
  epochs.corrected <- epochs %>% 
    dplyr::filter(electrode %in% target_channels) %>%
    left_join(regression_data, by=c("electrode", "word_type", "is_word")) %>%
    left_join(eyes, by=c("event_id", "t", "word_type", "is_word")) %>%
    mutate(v = v.x - (v.y * estimate))
  
  target.corrected <- epochs.corrected %>% 
    group_by(electrode, word_type, is_word, t) %>%
    summarize(v = mean(v))

  ggplot(target.corrected, aes(x=t, y=v, color=word_type, linetype=is_word)) +
    geom_line() +
    coord_cartesian(ylim=c(-20, 10)) +
    facet_grid(word_type~electrode)
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