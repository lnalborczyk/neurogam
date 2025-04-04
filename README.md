
# Precise temporal localisation of M/EEG effects with Bayesian generalised additive multilevel models

The goal of `neurogam` is to provide utilities for estimating onset and
offset of effects in M/EEG data.

## Installation

You can install the development version of `neurogam` from GitHub with:

``` r
install.packages("remotes")

remotes::install_github(
    repo = "https://github.com/lnalborczyk/neurogam",
    dependencies = TRUE
    )
```

## Usage

Below we fit a Bayesian generalised additive multilevel model (BGAMM) to
estimate the onset and offset of effect from simulated EEG data.

``` r
# loading the neurogam package
library(neurogam)

# importing some simulated EEG data
data(eeg_data)
head(eeg_data)

# fitting the BGAMM to identify clusters
results <- testing_through_time(data = eeg_data, threshold = 3)
```

``` r
# displaying the identified clusters
print(results$clusters)
#>   cluster_id cluster_onset cluster_offset
#> 1          1         0.172          0.348
```

``` r
# plotting the identified clusters superimposed with the data
plot(results)
```

<img src="man/figures/README-fig-clusters-1.png" width="100%" />
