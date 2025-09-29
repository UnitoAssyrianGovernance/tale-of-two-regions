# Filter Archaeological Sites by Time Period
#
# Filters a dataset of archaeological sites based on time period patterns,
# with options to exclude specific patterns and handle NA values.
#
# Inputs:
#   sites: An sf object containing archaeological site data. Must include columns: Period (character) and SizeHa (numeric)
#   period_pattern: A regular expression pattern to define periods to use
#   exclude_pattern: (Optional) A regular expression pattern to exclude specific periods
#   remove_na_size: Logical; if TRUE, removes sites with NA SizeHa values (default: TRUE)
#
# Details:
#   This function performs several operations:
#   1. Filters sites based on the period_pattern using regular expressions
#   2. Optionally excludes sites matching the exclude_pattern
#   3. Handles NA values in the SizeHa column according to the remove_na_size parameter
#   4. Orders the remaining sites by SizeHa in descending order
#   5. Removes duplicate sites based on their spatial geometry, keeping the largest site
#      (first occurrence since it was sorted by size before)
#
# Returns:
#   A filtered and ordered data frame or sf object containing only the sites
#   that match the specified criteria, with duplicates removed
filter_sites_by_period <- function(sites, period_pattern, exclude_pattern = NULL, remove_na_size = TRUE) {
  filtered <- sites[grepl(period_pattern, sites$Period), ]

  if (!is.null(exclude_pattern)) {
    filtered <- filtered[!grepl(exclude_pattern, filtered$Period), ]
  }

  if (!remove_na_size && any(is.na(filtered$SizeHa))) {
    warning("Some SizeHa values are NA and will be ordered last")
  }

  # Remove NAs if requested
  if (remove_na_size) {
    filtered <- filtered[!is.na(filtered$SizeHa), ]
  }

  # Order by size (with or without NAs)
  filtered <- filtered[order(filtered$SizeHa, decreasing = TRUE, na.last = TRUE), ]

  # Remove duplicates, keeping the first occurrence (which will be the largest due to ordering)
  filtered[!duplicated(sf::st_geometry(filtered)), ]
}
