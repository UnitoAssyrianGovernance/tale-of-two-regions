# Calculate Before Present (BP) Dates
#
# This function converts calendar dates to Before Present (BP) dates.
#
# Inputs:
#   df: A dataframe containing the date columns to convert
#   start_date: The start date column in calendar years (negative for BCE, positive for CE)
#   end_date: The end date column in calendar years
#   duration: Logical; if TRUE, calculate the duration between start and end dates (default: TRUE)
#
# Details:
#   The function converts calendar dates to BP dates using the formula BP = (calendar year * -1) + 1950.
#   The function adds new columns to the dataframe: StartBP, EndBP, and optionally Duration.
#
# Returns:
#   The original dataframe with additional columns for StartBP, EndBP, and optionally Duration
get_bp_dates <- function(df, start_date, end_date, duration = TRUE) {
  StartBP <- (start_date * -1) + 1950
  EndBP <- (end_date * -1) + 1950

  # add new columns for storing results
  df[, "StartBP"] <- StartBP
  df[, "EndBP"] <- EndBP

  if (duration == TRUE && "StartBP" %in% colnames(df)) {
    df[, "Duration"] <- (df$StartBP - df$EndBP)
  }

  print(str(df))

  return(df)
}
