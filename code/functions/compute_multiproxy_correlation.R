# Calculate Correlation Between Archaeological Proxies and SPD Data
#
# This function calculates Pearson correlation coefficients between
# multiple archaeological proxies and Summed Probability Distribution (SPD) data.
#
# Inputs:
#   archaeo_proxies: A list of archaeological proxies organized in time blocks
#   spd: A data frame containing SPD data with columns for dates (calBP/x)
#         and probability densities (PrDens/y). This should be the standard output
#         of rcarbon::spd()
#   max_bp: (Optional) A numeric value specifying the maximum "Before Present" date
#            to include in the analysis. If provided, only data from this date
#            and older will be used in the correlation.
#   row_names: (Optional) A character vector of custom names for the rows and
#               columns of the correlation matrix
#
# Details:
#   The function processes SPD data to align with the structure of archaeological proxies,
#   organizes the data into a combined dataframe, and calculates Pearson correlations.
#   It handles different column naming conventions in the SPD data and can filter by age.
#
# Returns:
#   A numeric matrix of correlation coefficients between all proxies and SPD data,
#   rounded to 2 decimal places
calculate_correlation <- function(archaeo_proxies, spd, max_bp = NULL, row_names = NULL) {

  # Determine which column names are used in the spd data frame
  cal_col <- if("calBP" %in% names(spd)) "calBP" else "x"
  pr_col <- if("PrDens" %in% names(spd)) "PrDens" else "y"

  # revert SPD to match archaeo_proxies order of data
  spd_rev <- rev(spd[[pr_col]][-length(spd[[pr_col]])])

  # nrow is equal to the block_width, ncol is equal to the length of a list item inside archaeo_proxies
  # the nrows variable is a small hackish way to get the block width without needing another function parameter
  # might need to improve this in the future
  nrows <- diff(as.numeric(gsub("start", "", names(archaeo_proxies[[1]][1:2]))))
  ncols <- length(archaeo_proxies[[1]])
  spd <- as.vector(colSums(matrix(spd_rev, nrow = nrows, ncol = ncols)))

  # Create a data frame with proxies
  # bind elements of a list with do.call
  proxies <- as.data.frame(cbind(do.call(cbind, archaeo_proxies), spd))

  if (is.null(max_bp)) {
    cor_result <- cor(proxies, method = "pearson")[, ]
    cor_result <- round(cor_result, digits = 2)
  } else {
    # If max_bp is provided
    # Filter data to remove dates later than X BP and exclude selected time-blocks
    row_index <- which(rownames(proxies) == paste0(max_bp, "start"))
    cor_result <- cor(proxies[row_index:nrow(proxies), ], method = "pearson")[, ]
    cor_result <- round(cor_result, digits = 2)
  }

  # Rename rows and columns if custom names are provided
  if (!is.null(row_names)) {
    rownames(cor_result) <- row_names
    colnames(cor_result) <- row_names
  }

  return(cor_result)
}
