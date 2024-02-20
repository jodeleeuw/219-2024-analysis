library(edfReader)
library(purrr)
library(dplyr)
library(jsonlite)
library(fuzzyjoin)

eeg.file <- "data/raw/eeg/subject-25-eeg_2024-02-14.bdf"
beh.file <- "data/raw/beh/219_2024_behavioral_25.json"

extract_events <- function(eeg.file, beh.file){
  
  # read in EEG data
  head <- readEdfHeader(eeg.file)
  signals <- readEdfSignals(head, signals="Ordinary")
  
  eeg.data <- map_df(signals, "signal")
  eeg.data$sample_id <- 1:nrow(eeg.data)
  
  ## read in behavioral data
  
  behavioral.data <- fromJSON(beh.file)
  trials <- behavioral.data %>% 
    dplyr::filter(task=="response") %>% 
    dplyr::filter(phase!="practice") %>%
    select(word_type, is_word, correct, time_elapsed) %>%
    mutate(time_since_first_trial = (time_elapsed - min(time_elapsed))/1000) %>%
    select(-time_elapsed) %>%
    mutate(row_id = 1:n())
  
  ## extract EEG events
  
  # 256 = moral / fashion word
  # 512 = moral / fashion non-word
  # 128 = non-moral / non-fashion word
  # 768 = non-moral / non-fashion non-word
  
  events <- eeg.data %>% 
    select(sample_id, TRIGGER) %>% 
    dplyr::filter(TRIGGER != lag(TRIGGER)) %>%
    dplyr::filter(TRIGGER %in% c(128, 256, 512, 768)) %>%
    mutate(time = sample_id/500)
  
  # find all rows where the time between neighboring rows is less than 0.5 seconds
  duplicate.events <- events %>% 
    dplyr::filter((time - lag(time) < 0.5) | (lead(time) - time < 0.5)) %>%
    dplyr::filter(TRIGGER != 768) %>% 
    pull(sample_id)
  
  filtered.events <- events %>% 
    dplyr::filter(!sample_id %in% duplicate.events) %>%
    mutate(event_type = case_when(
      TRIGGER == 256 ~ "moral_word",
      TRIGGER == 512 ~ "moral_nonword",
      TRIGGER == 128 ~ "nonmoral_word",
      TRIGGER == 768 ~ "nonmoral_nonword"
    )) %>%
    mutate(row_id = 1:n()) 
  
  # check if the first event is way before the others,
  # indicating that the experimenter dragged the window 
  # under the sensor after starting recording
  
  while(filtered.events$time[2] - filtered.events$time[1] > 10){
    # remove the first row
    filtered.events <- filtered.events %>%
      slice(-1) %>%
      mutate(row_id = 1:n())
  }
  
  removed_one <- TRUE
  while(removed_one){
    # remove the first row where diff is less than 1.9
    fake.event <- filtered.events %>% 
     mutate(diff = time - lag(time)) %>%
      dplyr::filter(diff < 1.9) %>%
      slice(1) %>%
      pull(sample_id)
    
    if(length(fake.event) == 0){
      removed_one <- FALSE
    }
    
    filtered.events <- filtered.events %>%
      dplyr::filter(!sample_id %in% fake.event) %>%
      mutate(row_id = 1:n())
  }
  
  # add a check that compares time of the event to the time in the behavioral
  # data file to align the two datasets when there are fewer events recorded
  # than trials in the behavioral data.
  
  filtered.events <- filtered.events %>% 
    mutate(time_since_first_trial = time - min(time))
  
  final_events <- trials %>%
    difference_left_join(filtered.events, by="time_since_first_trial", max_dist = 2) %>%
    mutate(event_id = 1:n()) %>%
    dplyr::filter(!is.na(sample_id)) %>%
    select(event_id, sample_id, time, event_type, word_type, is_word, correct)
  
 
  # final_events <- filtered.events %>% 
  #   left_join(trials, by="row_id") %>%
  #   mutate(event_id = 1:n()) %>%
  #   select(event_id, sample_id, time, event_type, word_type, is_word, correct)
  
  return(final_events)
}
