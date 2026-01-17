# code to prepare the simulated EEG data
eeg_data <- read.csv(file = here::here("./data-raw/raw_df.csv") )

# exporting to .rda
usethis::use_data(eeg_data, overwrite = TRUE, compress = "xz")

# code to prepare the cross-temporal data
timegen_data <- read.csv(file = here::here("./data-raw/timegen_data.csv") )

# exporting to .rda
usethis::use_data(timegen_data, overwrite = TRUE, compress = "xz")
