# Summarize Archaeological Proxies
#
# This function calculates summary statistics for different archaeological proxies
# by summing values across time blocks.
#
# Inputs:
#   site_counts: A matrix of site counts with rows for sites and columns for time blocks
#   aorist_weights: A matrix of aoristic weights with rows for sites and columns for time blocks
#   site_area_sums: A matrix of site areas with rows for sites and columns for time blocks
#   simulated_df: A matrix of simulated radiocarbon dates with rows for simulations and columns for time blocks
#
# Details:
#   The function calculates the column sums for each proxy matrix, this is a necessary step before normalizing the proxies.
#   For site areas, NA values are removed before summing.
#
# Returns:
#   A list containing:
#     - site_blocksum_counts: Vector of summed site counts for each time block
#     - site_blocksum_weights: Vector of summed aoristic weights for each time block
#     - site_blocksum_area: Vector of summed site areas for each time block
#     - site_simulated: Vector of summed simulated radiocarbon values for each time block
#
summarise_proxies <- function(site_counts, aorist_weights, site_area_sums, simulated_df) {
  site_blocksum_counts <- colSums(site_counts)
  site_blocksum_weights <- colSums(aorist_weights)
  site_blocksum_area <- colSums(site_area_sums, na.rm = TRUE)
  site_simulated <- colSums(simulated_df)

  return(list(
    site_blocksum_counts = site_blocksum_counts,
    site_blocksum_weights = site_blocksum_weights,
    site_blocksum_area = site_blocksum_area,
    site_simulated = site_simulated
  ))
}
