# Generate Time Slices for Paleoclimate and Archaeological Proxies
#
# This function creates time slices from paleoclimate data and various archaeological
# proxies for temporal analysis.
#
# Inputs:
#   paleoclim: A dataframe containing paleoclimate data
#   cname_paleoclim: (character) The column name in paleoclim containing the climate values
#   cname_bp: (character) The column name containing BP (Before Present) dates
#   spd: A Summed Probability Distribution object from rcarbon::spd()
#   aorist: A named numeric vector of aoristic weights with time period names
#   simulated: A named numeric vector of simulated radiocarbon data with time period names
#   area: A named numeric vector of estimated area size with time period names
#   counts: A named numeric vector of site counts with time period names
#   start: (numeric) The start date for analysis (BP)
#   end: (numeric) The end date for analysis (BP)
#   range: (numeric) The number of years to include in each time slice
#
# Details:
#   The function divides the time period between start and end into equal slices.
#   For each proxy, it calculates summary statistics for each time slice:
#   - For paleoclimate data: the mean value (inverted)
#   - For SPD and other proxies: the sum within each time slice
#   If the time range doesn't divide evenly by the slice range, the function
#   adjusts the end date and issues a warning, suggesting another end date.
#
# Returns:
#   A list containing:
#     - spd: Summed SPD values for each time slice
#     - paleoclim: Mean (inverted) paleoclimate values for each time slice
#     - spd_ts: Time-slice matrix for SPD data
#     - paleoclim_ts: Time-slice matrix for paleoclimate data
#     - aorist_sums: Summed aoristic weights for each time slice
#     - simulated_sums: Summed simulated data for each time slice
#     - area_sums: Summed area values for each time slice
#     - counts_sums: Summed count values for each time slice
gen_time_slices <- function(paleoclim, cname_paleoclim, cname_bp, spd, aorist, simulated, area, counts, start, end, range) {
  # ensure numeric values for start and end
  start <- as.numeric(start)
  end <- as.numeric(end)

  dummy_bp <- data.frame(bp = start:end)

  # Perform a left outer join: Return all rows from the left table,
  # and any rows with matching keys from the right table
  paleoclim_merged <- merge(dummy_bp, paleoclim, by = cname_bp, all.x = TRUE)

  # remove duplicated rows if any
  if(any(duplicated(paleoclim_merged$bp))) {
    paleoclim_merged <- paleoclim_merged[!duplicated(paleoclim_merged$bp), ]
    cat("After removing duplicates:", nrow(paleoclim_merged), "\n")
  }

  # Check if the data length is a sub-multiple or multiple of the number of rows (range) for the matrix creation
  if (nrow(paleoclim_merged) %% range != 0) {
  # Find how many complete blocks we can fit
  complete_blocks <- floor(nrow(paleoclim_merged) / range)

  # Calculate the row index for the last complete block
  last_row <- complete_blocks * range

  # Get the BP year at that position
  suggested_bp <- paleoclim_merged[last_row, "bp"]

  warning(paste0(
    "\nComputed data length is not a submultiple or multiple of the number of rows \n",
    "This means there is a mismatch between the length of the input vector and the size of the desired matrix \n",
    "Please consider replacing the end argument with: ", suggested_bp, " BP\n"
  ))

  # Adjust data frame to have correct number of rows and continue
  paleoclim_merged <- paleoclim_merged[1:last_row, ]
}

  # Create time slice labels
  bp_matrix <- matrix(paleoclim_merged$bp, nrow = range)

  # Retrieve the max and min of each matrix column
  max_bp <- apply(bp_matrix, 2, max)
  min_bp <- apply(bp_matrix, 2, min)
  # Generate a new character vector, each element represent the time-slice in the format min_bp - max_bp
  time_slice_labels <- sapply(seq_along(max_bp), function(x) {
    paste0(min_bp[x], "-", max_bp[x])
  })

  # generate time slice of z-scores
  paleoclim_time_slices <- as.data.frame(matrix(paleoclim_merged[[cname_paleoclim]], nrow = range))
  colnames(paleoclim_time_slices) <- time_slice_labels


  # Average the climate proxy value (z-score) within a N years time-slice
  paleoclim_mean <- colMeans(paleoclim_time_slices, na.rm = TRUE) * -1

  # Sum the SPD within a N years time-slice
  # select period of interest
  spd_df <- as.data.frame(spd$grid)
  spd_sel <- c(spd_df[which(spd_df$calBP >= start & spd_df$calBP <= end), ])
  # Generate time-slice of spd
  # revert is necessary so that we start from the most recent years
  spd_time_slices <- as.data.frame(matrix(rev(spd_sel$PrDens), nrow = range))
  colnames(spd_time_slices) <- time_slice_labels

  spd_sums <- colSums(spd_time_slices)

  # Generate time slices of the aoristic weights
  # Use the original named vector names to generate the BP values column
  aorist_df <- as.data.frame(aorist)
  aorist_df$bp <- rownames(aorist_df)

  # Rearrange and rename columns
  aorist_df <- aorist_df[, c(2, 1)]
  colnames(aorist_df) <- c(cname_bp, "aorist")

  # Remove the "start" from the column names
  aorist_df[, cname_bp] <- as.numeric(gsub("start$", "", aorist_df[, cname_bp]))

  # Perform a left outer join: Return all rows from the left table, and any rows with matching keys from the right table
  aorist_df_merged <- merge(dummy_bp, aorist_df, by = cname_bp, all.x = TRUE)

  # select only values within range
  aoristic_sel <- c(aorist_df_merged[which(aorist_df_merged$bp >= start & aorist_df_merged$bp <= end), ])

  # Generate time-slice of spd
  # revert is necessary so that we start from the most recent years
  aorist_time_slices <- as.data.frame(matrix(aoristic_sel$aorist, nrow = range))
  colnames(aorist_time_slices) <- time_slice_labels

  aorist_sums <- colSums(aorist_time_slices, na.rm = TRUE)

  # Generate time slices of simulated data
  # Use the original named vector names to generate the BP values column
  simulated_df <- as.data.frame(simulated)
  simulated_df$bp <- rownames(simulated_df)

  # Rearrange and rename columns
  simulated_df <- simulated_df[, c(2, 1)]
  colnames(simulated_df) <- c(cname_bp, "aorist")

  # Remove the "start" from the column names
  simulated_df[, cname_bp] <- as.numeric(gsub("start$", "", simulated_df[, cname_bp]))

  # Perform a left outer join: Return all rows from the left table, and any rows with matching keys from the right table
  simulated_df_merged <- merge(dummy_bp, simulated_df, by = cname_bp, all.x = TRUE)

  # select only values within range
  simulated_sel <- c(simulated_df_merged[which(simulated_df_merged$bp >= start & simulated_df_merged$bp <= end), ])

  # Generate time-slice of spd
  # revert is necessary so that we start from the most recent years
  simulated_time_slices <- as.data.frame(matrix(simulated_sel$aorist, nrow = range))
  colnames(simulated_time_slices) <- time_slice_labels

  simulated_sums <- colSums(simulated_time_slices, na.rm = TRUE)

  # Generate time slices of the area data
  area_df <- as.data.frame(area)
  area_df$bp <- rownames(area_df)

  # Rearrange and rename columns
  area_df <- area_df[, c(2, 1)]
  colnames(area_df) <- c(cname_bp, "area")

  # Remove the "start" from the column names
  area_df[, cname_bp] <- as.numeric(gsub("start$", "", area_df[, cname_bp]))

  # Perform a left outer join
  area_df_merged <- merge(dummy_bp, area_df, by = cname_bp, all.x = TRUE)

  # Select only values within range
  area_sel <- c(area_df_merged[which(area_df_merged$bp >= start & area_df_merged$bp <= end), ])

  # Generate time-slice of area
  area_time_slices <- as.data.frame(matrix(area_sel$area, nrow = range))
  colnames(area_time_slices) <- time_slice_labels

  area_sums <- colSums(area_time_slices, na.rm = TRUE)

  # Generate time slices of the counts data
  counts_df <- as.data.frame(counts)
  counts_df$bp <- rownames(counts_df)

  # Rearrange and rename columns
  counts_df <- counts_df[, c(2, 1)]
  colnames(counts_df) <- c(cname_bp, "counts")

  # Remove the "start" from the column names
  counts_df[, cname_bp] <- as.numeric(gsub("start$", "", counts_df[, cname_bp]))

  # Perform a left outer join
  counts_df_merged <- merge(dummy_bp, counts_df, by = cname_bp, all.x = TRUE)

  # Select only values within range
  counts_sel <- c(counts_df_merged[which(counts_df_merged$bp >= start & counts_df_merged$bp <= end), ])

  # Generate time-slice of counts
  counts_time_slices <- as.data.frame(matrix(counts_sel$counts, nrow = range))
  colnames(counts_time_slices) <- time_slice_labels

  counts_sums <- colSums(counts_time_slices, na.rm = TRUE)

  # Add the new elements to the return list
  return(list(
    spd = spd_sums,
    paleoclim = paleoclim_mean,
    spd_ts = spd_time_slices,
    paleoclim_ts = paleoclim_time_slices,
    aorist_sums = aorist_sums,
    simulated_sums = simulated_sums,
    area_sums = area_sums,
    counts_sums = counts_sums
  ))
}
