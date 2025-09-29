# Calculate Running Correlation Between Paleoclimate and Archaeological Proxies
#
# This function calculates running (moving window) correlations between paleoclimate data
# and various archaeological proxies over time.
#
# Inputs:
#   paleoclim: A numeric vector of paleoclimate values with time period names
#   spd: A numeric vector of Summed Probability Distribution values
#   aorist: (Optional) A numeric vector of aoristic sum values
#   simulated: (Optional) A numeric vector of simulated radiocarbon data values
#   area: (Optional) A numeric vector of area values
#   counts: (Optional) A numeric vector of count values
#   blockwidth: The width of each time block in years
#   start: The start date for analysis (BP)
#   end: The end date for analysis (BP)
#   method: The correlation method to use (e.g., "pearson", "spearman")
#          note: The function has been tested with pearson correlation only
#   cor_width: The number of time blocks to include in each correlation window
#
# Details:
#   The function computes running correlations between paleoclimate data and each provided
#   archaeodemographic proxy. For each correlation window, it calculates both the correlation
#   coefficient and its p-value. The function handles missing values by limiting analysis
#   to the last available complete data point. All input vectors must have the same length.
#   For background, see:
#   Roberts, N. (2021). Boon or Curse? The Role of Climate Change in the Rise and
#   Demise of Anatolian Civilizations, in WINDS OF CHANGE Environment and Society in Anatolia,
#   edited by Christopher H Roosevelt and John Haldon. pp.5–35.
#
# Returns:
#   A list containing:
#     - spd_all_corr: Data frame of all SPD correlations with start/end dates, coefficients and p-values
#     - spd_significant_corr: Data frame of significant SPD correlations (p < 0.05)
#     - Similar output for aorist, simulated, area, and counts (if provided)
#     - moving_window: The total time span in years covered by each correlation window
running_cor <- function(paleoclim, spd, aorist = NULL, simulated = NULL, area = NULL, counts = NULL, blockwidth, start, end, method, cor_width) {
  # Ensure numeric values
  start <- as.numeric(start)
  end <- as.numeric(end)

  # Check for matching lengths
  if (
    length(paleoclim) != length(spd) ||
      (!is.null(aorist) && length(paleoclim) != length(aorist)) ||
      (!is.null(simulated) && length(paleoclim) != length(simulated)) ||
      (!is.null(area) && length(paleoclim) != length(area)) ||
      (!is.null(counts) && length(paleoclim) != length(counts))
  ) {
    stop("All provided input vectors must have the same length")
  }

  if (cor_width >= length(paleoclim)) {
    stop("Correlation width must be less than the length of the input vectors")
  }

  # Check for missing values in the paleoclim database
  if (any(is.na(as.numeric(paleoclim)))) {
    # if there are, limit the analysis on the last useful point
    last_index <- min(which(is.na(paleoclim))) - 1
    #  extract the name of the last useful interval
    last_range <- names(paleoclim[last_index])
    # Extract the string after the "-" in the col vector and convert to numeric
    end_new <- as.numeric(strsplit(last_range, "-")[[1]][2])

    warning(paste0("\nNAs found \nAnalysis has been limited up until ", end_new, " BP \n"))
  } else {
    end_new <- end
    last_index <- length(paleoclim)
  }

  # Extract time ranges from names of paleoclim
  time_ranges <- names(paleoclim)
  range_breaks <- numeric(length(time_ranges) * 2)

  for (i in seq_along(time_ranges)) {
    parts <- as.numeric(strsplit(time_ranges[i], "-")[[1]])
    range_breaks[i * 2 - 1] <- parts[1] # start of range
    range_breaks[i * 2] <- parts[2] # end of range
  }

  range_breaks <- sort(unique(range_breaks))

  # Helper function to calculate running correlation and p-values
  get_running_cor <- function(proxy, proxy_name) {
    # Show basic info about what we're calculating
    cat(paste0(
      "\nCalculating ",
      proxy_name,
      " correlation (",
      length(proxy),
      " time periods, window size: ",
      cor_width,
      " blocks)...\n"
    ))

    # Calculate running correlation
    run_cor <- gtools::running(
      proxy[1:last_index],
      paleoclim[1:last_index],
      fun = cor,
      method = method,
      width = cor_width
    )

    # Generate matrix for p-values
    pvalue <- numeric(length(run_cor))

    # Create and use progress bar
    pb <- txtProgressBar(min = 0, max = length(run_cor), style = 3)

    # Loop for calculating p-values
    block_start <- cor_width - 1
    for (i in seq_along(run_cor)) {
      # Calculate correlation test
      test <- cor.test(
        proxy[i:c(block_start + i)],
        paleoclim[i:c(block_start + i)],
        method = method
      )

      pvalue[i] <- test$p.value

      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)

    # Create breaks for plotting
    n_cors <- length(run_cor)
    window_size <- cor_width

    # Create central point of each correlation window
    midpoints <- numeric(n_cors)
    for (i in 1:n_cors) {
      # Get range names for the first and last points in this window
      first_range <- names(paleoclim)[i]
      last_range <- names(paleoclim)[i + window_size - 1]

      # Extract start of first range and end of last range
      first_start <- as.numeric(strsplit(first_range, "-")[[1]][1])
      last_end <- as.numeric(strsplit(last_range, "-")[[1]][2])

      # Use midpoint of entire window
      midpoints[i] <- (first_start + last_end) / 2
    }

    # Calculate start and end points from midpoints
    half_width <- (blockwidth * cor_width) / 2
    brks_start <- midpoints - half_width
    brks_end <- midpoints + half_width

    # Store results in data frame
    all_corr <- data.frame(
      brks_start = floor(brks_start),
      brks_end = floor(brks_end),
      run_cor = run_cor,
      pvalue = pvalue
    )

    # Get significant correlations
    significant_corr <- subset(all_corr, pvalue < 0.05)

    # Show completion info
    cat(paste0(
      "Completed ",
      proxy_name,
      " correlation (",
      nrow(significant_corr),
      " significant results)\n"
    ))

    return(list(all_corr = all_corr, significant_corr = significant_corr))
  }

  # Calculate correlations for each proxy
  spd_results <- get_running_cor(spd, "SPD")

  # Calculate correlations for aoristic data if providedd
  if (!is.null(aorist)) {
    aorist_results <- get_running_cor(aorist, "Aorist")
  } else {
    aorist_results <- list(all_corr = NULL, significant_corr = NULL)
  }

  # Calculate correlations for simulated data if provided
  if (!is.null(simulated)) {
    sim_results <- get_running_cor(simulated, "Simulated")
  } else {
    sim_results <- list(all_corr = NULL, significant_corr = NULL)
  }

  # Calculate correlations for area if provided
  if (!is.null(area)) {
    area_results <- get_running_cor(area, "Area")
  } else {
    area_results <- list(all_corr = NULL, significant_corr = NULL)
  }

  # Calculate correlations for counts if provided
  if (!is.null(counts)) {
    counts_results <- get_running_cor(counts, "Counts")
  } else {
    counts_results <- list(all_corr = NULL, significant_corr = NULL)
  }

  # Add the new results to the return list
  return(list(
    spd_all_corr = spd_results$all_corr,
    spd_significant_corr = spd_results$significant_corr,
    aorist_all_corr = aorist_results$all_corr,
    aorist_significant_corr = aorist_results$significant_corr,
    sim_all_corr = sim_results$all_corr,
    sim_significant_corr = sim_results$significant_corr,
    area_all_corr = area_results$all_corr,
    area_significant_corr = area_results$significant_corr,
    counts_all_corr = counts_results$all_corr,
    counts_significant_corr = counts_results$significant_corr,
    moving_window = blockwidth * cor_width
  ))
}
