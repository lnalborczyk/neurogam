# Time-generalisation decoding example dataset

A tidy dataset containing time-generalisation decoding results (AUC
values) for multiple participants, indexed by training time and testing
time.

## Usage

``` r
timegen_data
```

## Format

A data frame with the following columns:

- train_time:

  Numeric. Training time (in seconds) relative to event onset.

- test_time:

  Numeric. Testing time (in seconds) relative to event onset.

- auc:

  Numeric. Area under the ROC curve (AUC) for the decoding model.

- participant:

  Integer. Participant identifier.

## Details

Each row corresponds to one participant's AUC for a given (train_time,
test_time) pair. This dataset can be used to demonstrate aggregation,
visualisation, and statistical modelling of time-generalization
matrices.
