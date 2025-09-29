
# Functions for Montecarlo Mimulation on Archaeological Sites

# Simulate Archaeological Site Durations
#
# Generates random duration values for archaeological sites based on a normal distribution
#
# Inputs:
#   sites: A data frame containing site information
#   nsim: Number of simulations to run
#   mean_duration: Mean duration in years for the normal distribution
#   sd_duration: Standard deviation in years for the normal distribution
#   min_duration: Minimum allowed duration for sites
#
# Details:
#   Generates random durations from a normal distribution.
#
# Returns:
#   A vector of simulated site durations in years
simulate_site_durations <- function(sites, nsim, mean_duration, sd_duration, min_duration) {
  cat("\nGenerating random durations:\n")

  pb <- txtProgressBar(min = 0, max = nsim, style = 3, width = 100)

  for (b in 1:nsim) {
    setTxtProgressBar(pb, b)
    # cat(paste(b, "; ", sep = ""))
    site_dur <- rnorm(n = nrow(sites), mean = mean_duration, sd = sd_duration)
    while (any(site_dur < min_duration)) {
      site_dur[site_dur < min_duration] <- rnorm(n = length(site_dur[site_dur < min_duration]), mean = mean_duration, sd = sd_duration) # enforce truncation at min_duration
    }
  }
  close(pb)
  cat("\n")

  return(site_dur)
}

# Count Sites Across Time Blocks Based on Simulated Durations
#
# Distributes sites across time blocks according to simulated durations
#
# Inputs:
#   sites: A data frame with site information including StartBP, EndBP, and Duration
#   nsim: Number of simulations to run
#   site_dur: Vector of simulated site durations in years (generated in the simulate_site_durations function)
#   time_blocks: Vector defining the time block boundaries
#   block_width: Width of each time block in years
#   new_columns: Vector of column names for the time blocks
#
# Details:
#   For each site, if the simulated duration is >= the actual duration, uses the
#   original start/end dates. Otherwise, randomly places the site within its
#   original date range. Counts are distributed proportionally across time blocks.
#
# Returns:
#   A data frame where each row represents a simulation and columns represent
#   time blocks with summed site counts
# Notes: find another way instead of site_counts or sites to generate the empty df
count_simulated_sites <- function(sites, nsim, site_dur, time_blocks, block_width, new_columns) {
  sim_df <- as.data.frame(matrix(ncol = length(new_columns), nrow = nsim))
  cat("\nRunning Montecarlo Simulation:\n")
  pb <- txtProgressBar(min = 0, max = nsim, style = 3, width = 100)
  colnames(sim_df) <- new_columns

  for (b in 1:nsim) {
    tmp_site_counts <- sites
    tmp_site_counts[new_columns] <- NA
    tmp_site_counts <- tmp_site_counts[, new_columns]
    # cat(paste(b, "; ", sep = ""))
    setTxtProgressBar(pb, b)

    for (a in seq_len(nrow(sites))) {
      if (site_dur[a] >= sites$Duration[a]) {
        site_start <- sites$StartBP[a]
        site_end <- sites$EndBP[a]
      } else {
        site_start <- round(runif(1, max = sites$StartBP[a], min = sites$EndBP[a] + site_dur[a]), 0)
        site_end <- round(site_start - site_dur[a])
      }

      site_years <- seq(site_end, site_start, by = 1)
      site_years_hist <- hist(site_years, breaks = time_blocks, plot = FALSE)
      tmp_site_counts[a, ] <- site_years_hist$counts / block_width
    }

    sim_df[b, ] <- colSums(tmp_site_counts, na.rm = TRUE)
  }
  close(pb)
  cat("\n")

  return(sim_df)
}

# Run Monte Carlo Simulation for Archaeological Site Durations
#
# Wrapper function to perform a Monte Carlo simulation for archaeological site durations
#
# Inputs:
#   sites: A data frame with site information
#   nsim: Number of simulations to run
#   mean_duration: Mean duration in years for the normal distribution
#   sd_duration: Standard deviation in years for the normal distribution
#   min_duration: Minimum allowed duration for sites
#   time_blocks: Vector defining the time block boundaries
#   block_width: Width of each time block in years
#   new_columns: Vector of column names for the time blocks
#
# Details:
#   This function combines the site duration simulation and site counting
#   steps into a complete Monte Carlo workflow. It simulates random site
#   durations and distributes sites across time blocks accordingly.
#
# Returns:
#   A list containing:
#     - site_dur: Vector of simulated site durations in years (generated from the simulate_site_durations function)
#     - sim_df: Data frame of simulation results across time blocks (generated from the count_simulated_sites function)
montecarlo_sim <- function(sites, nsim, mean_duration, sd_duration, min_duration, time_blocks, block_width, new_columns) {
  site_dur <- simulate_site_durations(sites, nsim, mean_duration, sd_duration, min_duration)

  sim_df <- count_simulated_sites(sites, nsim, site_dur, time_blocks, block_width, new_columns)

  return(list(site_dur = site_dur, sim_df = sim_df))
}
