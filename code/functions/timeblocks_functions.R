# Functions for Creating and Managing Time Blocks
# NOTE: Create a more general wrapper function that combines all timeblock functions

# Generate Time Blocks
#
# Creates time blocks based on a specified time window and block width
#
# Inputs:
#   time_window: A named list with elements "start" and "end" defining the time range
#   block_width: Width of each time block in years
#
# Details:
#   This function divides the specified time window into time blocks.
#   Each block is named with its starting year followed by "start".
#
# Returns:
#   A list containing:
#     - time_blocks: Numeric vector of time blocks start
#     - new_columns: Character vector of column names for the blocks
#
gen_time_blocks <- function(time_window, block_width) {
  time_blocks <- seq(time_window[["start"]], time_window[["end"]], by = block_width)
  new_columns <- paste(time_blocks[1:(length(time_blocks) - 1)], "start", sep = "")

  return(list(time_blocks = time_blocks, new_columns = new_columns))
}

# Generate Empty Tables for Time Block Analysis
#
# Creates empty data frames for storing site data across time blocks
#
# Inputs:
#   sites: A data frame containing site information
#   new_columns: A character vector of column names for time blocks (from gen_time_blocks)
#
# Details:
#   This function prepares three empty data frames to be filled with:
#   1. Raw site counts per time block
#   2. Aoristic weights per time block
#   3. Site area sums per time block
#   Each data frame will have rows equal to the number of sites and columns
#   equal to the number of time blocks.
#
# Returns:
#   A list containing three empty data frames:
#     - site_counts: For storing raw counts
#     - aorist_weights: For storing aoristic weights
#     - site_area_sums: For storing area calculations
#
# NOTE: Consider if it's worth turning into a vectorized operation using apply functions
gen_new_tables <- function(sites, new_columns) {
  aorist_weights <- site_area_sums <- site_counts <- sites

  site_counts[new_columns] <- NA
  site_counts <- site_counts[, new_columns]

  aorist_weights[new_columns] <- NA
  aorist_weights <- aorist_weights[, new_columns]

  site_area_sums[new_columns] <- NA
  site_area_sums <- site_area_sums[, new_columns]

  list(
    site_counts = site_counts,
    aorist_weights = aorist_weights,
    site_area_sums = site_area_sums
  )
}

# Update Time Block Tables with Archaeological Data
#
# Fills prepared data frame from gen_new_tables() with site data distributed across time blocks
#
# Inputs:
#   site_start: (numeric) vector of site start dates (BP)
#   site_end: (numeric) vector of site end dates (BP)
#   time_blocks: Numeric vector of time block break points
#   block_width: (integer) Width of each time block in years
#   site_counts: Empty data frame for storing site counts
#   aorist_weights: Empty data frame for storing aoristic weights
#   site_area_sums: Empty data frame for storing site area calculations
#   area_field: (numeric) vector of site areas
#
# Details:
#   For each site, this function:
#   1. Generates all years between start and end dates
#   2. Distributes these years across time blocks
#   3. Calculates raw counts
#   4. Calculates aoristic weights
#   5. Calculates area values
#
# Returns:
#   A list containing the three updated data frames:
#     - site_counts: Raw site counts per time block
#     - aorist_weights: Aoristic weights per time block
#     - site_area_sums: Site area calculations per time block
update_time_blocks <- function(site_start, site_end, time_blocks, block_width, site_counts, aorist_weights, site_area_sums, area_field) {

  cat("\nUpdating time blocks with archaeological proxies:\n")
  pb <- txtProgressBar(min = 0, max = nrow(site_counts), style = 3, width = 100)
  for (i in seq_along(site_start)) {
    setTxtProgressBar(pb, i)
    # cat(paste(i, "; ", sep = ""))

    site_years <- seq(site_end[i], site_start[i], by = 1)
    site_years_hist <- hist(site_years, breaks = time_blocks, plot = FALSE)
    time_wts <- site_years_hist$counts / block_width
    site_counts[i, ] <- time_wts
    aorist_weights[i, ] <- time_wts / sum(time_wts)
    site_area_sums[i, ] <- time_wts * area_field[i]
  }
  close(pb)
  cat("\n")

  list(
    site_counts = site_counts,
    aorist_weights = aorist_weights,
    site_area_sums = site_area_sums
  )
}
