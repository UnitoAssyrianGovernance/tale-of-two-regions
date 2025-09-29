# Calculate Z-Score
#
# This function calculates the z-score of data in a dataframe column.
# Inputs:
#   df: A dataframe containing the data to be standardized
#   cname: (character) The name of the column in df containing the data to standardize
#   cname_new: (character) The name for the new column that will contain the z-scores
#
# Details:
#   The function calculates z-scores using the formula z = (x - μ) / σ,
#   where x is the value, μ is the mean, and σ is the standard deviation.
#   This standardizes the data to have a mean of 0 and a standard deviation of 1.
#   See Finné, M., Woodbridge, J., Labuhn, I., Roberts, C.N., 2019.
#   Holocene hydro-climatic variability in the Mediterranean: A synthetic multi-proxy reconstruction.
#   The Holocene 29, 847–863. https://doi.org/10.1177/0959683619826634
#
# Returns:
#   The original dataframe with an additional column containing the z-scores
get_z_score <- function(df, cname, cname_new) {

  # Remove d180 NA from df (if any)
  df <- df[!is.na(df[[cname]]),]

  z <- (df[[cname]] - mean(df[[cname]])) / sd(df[[cname]])

  df[[cname_new]] <- z
  return(df)
}
