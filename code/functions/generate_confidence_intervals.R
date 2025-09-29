# Calculate Confidence Intervals for Simulated Data
#
# This function calculates the 95% confidence intervals and median values
# from a matrix of simulated values.
#
# Inputs:
#   simulated_df:  A matrix where rows represent different simulations and
#                  columns represent time periods
#
# Details:
#   The function calculates three summary statistics for each column (time period):
#   - Upper bound (97.5th percentile) of the 95% confidence interval
#   - Lower bound (2.5th percentile) of the 95% confidence interval
#   - Median value (50th percentile)
#
#   These are useful for comparing observed archaeological data
#   with simulated null models and for plotting.
#
# Returns:
#   A list containing:
#     - hicount: Vector of upper bounds (97.5th percentile) for each time period
#     - lowcount: Vector of lower bounds (2.5th percentile) for each time period
#     - medcount: Vector of median values for each time period
get_conf_intervals <- function(simulated_df) {
  get_high_conf <- function(y) {
    quantile(y, probs = 0.975)
  }

  # function defining lower limit of 95% confidence envelop
  get_low_conf <- function(y) {
    quantile(y, probs = 0.025)
  }

  hicount <- apply(simulated_df, 2, get_high_conf) # upper limit of 95% confidence envelope of randomised start date of sites
  lowcount <- apply(simulated_df, 2, get_low_conf) # lower limit of 95% confidence envelope of randomised start date of sites
  medcount <- apply(simulated_df, 2, median) # median of 95% confidence envelope of randomised start date of sites

  return(list(
    hicount = hicount,
    lowcount = lowcount,
    medcount = medcount
  ))
}
