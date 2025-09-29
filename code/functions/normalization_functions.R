# Collection of Functions for Normalizing Data

# Min-Max Normalization
#
# Scales a numeric vector to range from 0 to 1
#
# Inputs:
#   x:  A numeric vector to normalize
# Returns:
#   Normalized vector with values between 0 and 1
normalize_vector <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Normalize Multiple Archaeological Proxies
#
# Applies normalize_vector to three archaeological proxies
#
# Inputs:
#   count_proxy: A numeric vector of site counts
#   area_proxy: A numeric vector of site areas
#   weight_proxy: A numeric vector of aoristic weights
# Returns:
#   A data frame with three normalized columns: count, area, and weight
normalize_archaeo_proxies <- function(count_proxy, area_proxy, weight_proxy) {
  count_proxy_norm <- normalize_vector(count_proxy)
  area_proxy_norm <- normalize_vector(area_proxy)
  weight_proxy_norm <- normalize_vector(weight_proxy)
  return(data.frame(count = count_proxy_norm, area = area_proxy_norm, weight = weight_proxy_norm))
}

# Normalize SPD Data
#
# Normalizes a Summed Probability Distribution data frame
#
# Inputs:
#   spd: A data frame with calBP dates and probability values
# Returns:
#   A data frame with original dates (x) and normalized probability values (y)
# Note:
#   Column names could be changed to more descriptive ones in future versions
normalize_spd <- function(spd) {
  spdY <- normalize_vector(as.numeric(spd[, 2]))
  spdX <- as.numeric(spd[, "calBP"])
  return(data.frame(x = spdX, y = spdY))
}

# Normalize Simulation Data Frame
#
# Normalizes each row of a simulated proxy data frame
#
# Inputs:
#   sim_df: A data frame of simulation results with rows indicating the iterations
# Returns:
#   A normalized version of the input with values between 0 and 1
#   NA values are replaced with 0 to avoid computation errors in other operations
normalize_sim_df <- function(sim_df) {
  sim_df_norm <- sim_df
  for (a in 1:nrow(sim_df_norm)) {
    x <- as.numeric(sim_df[a, ])
    tmp_norm <- normalize_vector(x)
    sim_df_norm[a, ] <- tmp_norm
  }
  sim_df_norm[is.na(sim_df_norm)] <- 0
  return(sim_df_norm)
}
