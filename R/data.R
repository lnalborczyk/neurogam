#' Simulated EEG data
#'
#' Simulated EEG data.
#'
#' @format A data frame with 502000 rows and 5 variables:
#' \describe{
#'   \item{participant}{Character variable indicating the participant's ID.}
#'   \item{condition}{Integer variable indicating the condition.}
#'   \item{trial}{Integer variable indicating the trial number.}
#'   \item{time}{Numeric variable indicating the time step (in seconds).}
#'   \item{eeg}{Numeric variable indicating the EEG signal.}
#' }
#'
"eeg_data"

#' Time-generalisation decoding example dataset
#'
#' A tidy dataset containing time-generalisation decoding results (AUC values)
#' for multiple participants, indexed by training time and testing time.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{train_time}{Numeric. Training time (in seconds) relative to event onset.}
#'   \item{test_time}{Numeric. Testing time (in seconds) relative to event onset.}
#'   \item{auc}{Numeric. Area under the ROC curve (AUC) for the decoding model.}
#'   \item{participant}{Integer. Participant identifier.}
#' }
#'
#' @details
#' Each row corresponds to one participant's AUC for a given (train_time, test_time)
#' pair. This dataset can be used to demonstrate aggregation, visualisation, and
#' statistical modelling of time-generalization matrices.
#'
"timegen_data"
