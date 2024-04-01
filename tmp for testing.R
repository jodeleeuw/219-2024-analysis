source("preprocessing/preprocessing-eeg-functions.R")

subject <- "24"

eeg.file <- list.files('data/raw/eeg', pattern=paste0("subject-", subject), full.names = T)
beh.file <- paste0('data/raw/beh/219_2024_behavioral_', subject,'.json')

sampling_rate = 500
filter_low = 0.1
filter_high = 30
segment_begin = -400
segment_end = 1200
segment_offset = 0
bad_segment_range = 500
eye_threshold = 70
which_electrodes = c("Fp1", "Fp2", "Cz", "Pz")
eye_channels=c("Fp1", "Fp2")
target_channels=c("Cz", "Pz")

# epochs <- read_eeg_tidy(eeg.file, beh.file, which_electrodes) %>%
#   #linked_ears_rereference() %>%
#   bandpass(low=filter_low, high=filter_high, sampling_rate=sampling_rate) %>%
#   #notch(low=notch_filter_low, high=notch_filter_high, sampling_rate=sampling_rate) %>%
#   segment(start = segment_begin, end = segment_end, offset=segment_offset, sampling_rate=sampling_rate) %>%
#   baseline_correct()


data <- read_eeg_tidy(eeg.file, beh.file, which_electrodes)

sample_window <- 100000:105000

unfiltered_data_cz <- data$signals %>%
  dplyr::filter(electrode == "Cz") %>%
  dplyr::filter(sample_id %in% sample_window)

library(ggplot2)
ggplot(unfiltered_data_cz, aes(x=sample_id, y=v)) +
  geom_line() +
  labs(title="Unfiltered Cz data", x="Sample ID", y="Voltage (uV)")+
  theme_minimal()

data_filtered <- data %>%
  bandpass(low=filter_low, high=filter_high, sampling_rate=sampling_rate)

filtered_data_cz <- data_filtered$signals %>%
  dplyr::filter(electrode == "Cz") %>%
  dplyr::filter(sample_id %in% sample_window)

ggplot(filtered_data_cz, aes(x=sample_id, y=v)) +
  geom_line() +
  labs(title="Filtered Cz data", x="Sample ID", y="Voltage (uV)")+
  theme_minimal()

epochs <- data_filtered %>%
  segment(start = segment_begin, end = segment_end, offset=segment_offset, sampling_rate=sampling_rate)

example_epoch_cz <- epochs %>%
  dplyr::filter(electrode == "Cz") %>%
  dplyr::filter(event_id == 111)

ggplot(example_epoch_cz, aes(x=t, y=v)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title="Example epoch, Cz", x="Time (ms)", y="Voltage (uV)")+
  theme_minimal()+
  theme(panel.grid=element_blank())

epochs_first_removal <- epochs %>%
  large_anomaly_removal(1000) %>%
  filter_epochs_with_bad_channels(expected_channels = 4)

# find removed event_id
removed_event_id <- epochs %>% 
  anti_join(epochs_first_removal, by=c("event_id", "electrode")) %>%
  pull(event_id) %>%
  unique()

# plot removed epoch example

example_epoch_removed <- epochs %>%
  dplyr::filter(electrode == "Fp2") %>%
  dplyr::filter(event_id == 302)

ggplot(example_epoch_removed, aes(x=t, y=v)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title="Example epoch, Fp2, removed", x="Time (ms)", y="Voltage (uV)")+
  theme_minimal()+
  theme(panel.grid=element_blank())

# baseline correction example

epochs_baseline_corrected <- epochs_first_removal %>%
  baseline_correct(start_baseline=-200)

example_epoch_baseline_corrected <- epochs_baseline_corrected %>%
  dplyr::filter(electrode == "Cz") %>%
  dplyr::filter(event_id == 111)

ggplot(example_epoch_cz, aes(x=t, y=v)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title="Example epoch, Cz, baseline corrected", x="Time (ms)", y="Voltage (uV)")+
  annotate("rect", xmin = -200, xmax = 0, ymin = -50, ymax = 50, alpha = .2)+
  theme_minimal()+
  theme(panel.grid=element_blank())

ggplot(example_epoch_baseline_corrected, aes(x=t, y=v)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title="Example epoch, Cz, baseline corrected", x="Time (ms)", y="Voltage (uV)")+
  annotate("rect", xmin = -200, xmax = 0, ymin = -50, ymax = 50, alpha = .2)+
  theme_minimal()+
  theme(panel.grid=element_blank())

# eye blink removal 


example_epoch_fp1 <- epochs %>%
  dplyr::filter(electrode %in% c("Fp1", "Fp2")) %>%
  dplyr::filter(event_id == 111) %>%
  group_by(t) %>%
  summarize(v = mean(v))

ggplot(example_epoch_fp1, aes(x=t, y=v)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title="Example epoch, eye channel", x="Time (ms)", y="Voltage (uV)")+
  theme_minimal()+
  theme(panel.grid=element_blank())


epochs_eye_removed <- epochs_baseline_corrected %>%
  eyeblink_removal_with_regression(eye_channels=c("Fp1", "Fp2"), target_channels=c("Cz", "Pz"), blink_criteria=10, sampling_rate=sampling_rate, windows_ms=200)

# find where the eye blinks were removed

events_with_blinks <- epochs_baseline_corrected %>%
  left_join(epochs_eye_removed, by=c("event_id", "electrode", "t")) %>%
  dplyr::filter(v.x != v.y) %>%
  pull(event_id) %>%
  unique()

example_epoch_eye_removed <- epochs_eye_removed %>%
  dplyr::filter(electrode == "Cz") %>%
  dplyr::filter(event_id == 111)

ggplot(example_epoch_eye_removed, aes(x=t, y=v)) +
  geom_line(data=example_epoch_baseline_corrected, aes(x=t, y=v), color="red") +
  geom_line() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title="Example epoch, Cz, eye blink removed", x="Time (ms)", y="Voltage (uV)")+
  theme_minimal()+
  theme(panel.grid=element_blank())




