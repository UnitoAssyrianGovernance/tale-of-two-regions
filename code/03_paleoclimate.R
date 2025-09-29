## Process paleoclimate data, calculate z-scores and plot with radiocarbon data
## Last edited 2025-09-22
## Authors: Andrea Titolo (refactored code), Alessio Palmisano (original code)
## Packages version at the time of analysis
## here: 1.0.1

# Set-up -----------------------------------------------------------------------

# Source the code for plotting coloured rainbows
source(here::here("code/functions/calculate_z_score.R"))
source(here::here("code/functions/plot_bb_blocks.R"))

# Load required library
required_packages <- c("here")

# Check if renv is active
if ("renv" %in% loadedNamespaces()) {
  cat("renv is active\n")
  cat("All packages should be present already\n")
  cat("If you run into troubles, please install them manually\n")
} else {
  # If not, see if packages are installed on the system
  cat("renv not active\n")
  cat("Using system packages\n")
  new_packages <- required_packages[
    !(required_packages %in% installed.packages()[, "Package"])
  ]
  if (length(new_packages)) {
    cat(
      "Installing missing packages:",
      paste(new_packages, collapse = ", "),
      "\n"
    )
    install.packages(new_packages, repos = "https://cran.rstudio.com/")
  } else {
    cat("All required packages already installed\n")
  }
}

# Load all required packages
lapply(required_packages, library, character.only = TRUE)

# Load data -----------------------------------------------------------------------

# List files in the directory that contain the word "palaeoclim_"
files <- list.files(
  here::here("data/raw/csv"),
  pattern = "palaeoclim_",
  full.names = TRUE
)

# Load the files listed above inside a list
paleoclimate <- lapply(files, read.csv, header = TRUE, sep = ",")

# Extract proxy name from the original file and use it to name list elements
names(paleoclimate) <- gsub(
  ".csv",
  "",
  gsub("palaeoclim_", "", basename(files))
)

# Explore data -----------------------------------------------------------------------

# Calculate mean sampling interval for paleoclimate records in use
# This correspoond to the differences between consecutive time points
time_diffs <- lapply(paleoclimate, function(df) {
  # Sort by bp to ensure chronological order
  df <- df[order(df$bp), ]

  # Calculate differences between consecutive time points
  diff(df$bp)
})

# Calculate the mean intervals (useful for Table 2 in the paper)
mean_intervals <- data.frame(
  "Sampling Interval" = round(sapply(time_diffs, mean))
)

write.csv(
  mean_intervals,
  here::here("output/csv", "paleoclimate_mean_intervals.csv"),
  row.names = TRUE
)

# Z-score -----------------------------------------------------------------------

# Get the z-score for each dataframe in the list and add it as a new column
# Apply the function on all the d18o column of each dataframe in the list
paleoclim_z <- lapply(
  seq_along(paleoclimate),
  function(x) get_z_score(df = paleoclimate[[x]], cname = "d18o", cname_new = "z_score")
)

# Make sure list elements names are the same as the original
names(paleoclim_z) <- names(paleoclimate)

# Data manipulation -----------------------------------------------------------------------

# Identify the element in the list with more rows with nrow
row_max_n <- which.max(sapply(paleoclim_z, nrow))

# Make sure it is positioned at the beginning of the list
paleoclim <- append(paleoclim_z[-row_max_n], paleoclim_z[row_max_n], after = 0)

# Rename each d18o column of each dataframe in the paleoclim list to the name of each list element
# This will help with the merge later
for (i in seq_along(paleoclim)) {
  colnames(paleoclim[[i]])[colnames(paleoclim[[i]]) == "d18o"] <- names(paleoclim)[i]
  colnames(paleoclim[[i]])[colnames(paleoclim[[i]]) == "z_score"] <- paste0("z_score_", names(paleoclim)[i])
}

# Merge all paleoclimate dataframes into a single unified dataframe
# Reduce() applies merge() sequentially: merge(df1, df2), then merge(result, df3), etc.
# Because we made sure the longest df was first, the rows are aligned correctly
# all = TRUE keeps all rows from both dataframes (full outer join)
paleoclim_long <- Reduce(function(x, y) merge(x, y, by = "bp", all = TRUE), paleoclim)

# Save processed csv
write.csv(x = paleoclim_long, file = here::here("data/processed/csv", "paleoclimate_all.csv"), row.names = FALSE)

# Plotting -----------------------------------------------------------------------

# Plot with spd -----------------------------------------------------------------------

# Load the results of the radiocarbon analyses
load(here::here("output/rda/radiocarbon_data.RData"))

# Plot palaeoclimate data
png(file = here::here("output/png", "fig-07.png"), width = 1600, height = 1800, res = 150)
layout(matrix(c(1, 2, 3, 4, 5), 5, 1, byrow = TRUE), widths = 6, heights = c(1, 1, 1, 1, 2.5))
par(mar = c(0, 4, 1, 3)) # c(bottom, left, top, right)
# layout.show(n=5)

# Make sure we are only plotting existing values
bp <- paleoclim_long$bp[!is.na(paleoclim_long$z_score_jeita_cave)]
z_score <- paleoclim_long$z_score_jeita_cave[!is.na(paleoclim_long$z_score_jeita_cave)]
ylim <- c(max(z_score, na.rm = TRUE), min(z_score, na.rm = TRUE))

par(mar = c(0, 4, 1, 3)) # c(bottom, left, top, right)
plot(bp, z_score, type = "l", axes = FALSE, xlim = c(6500, 2200), ylim = ylim, lwd = 1.5, col = "darkgreen", xaxt = "n", ylab = "", xaxs = "i", yaxs = "i")
axis(2, at = seq(2.5, -2.5, -0.5), labels = seq(2.5, -2.5, -0.5), lwd = 1, line = 0.2, cex.axis = 1, mgp = c(0, 0.4, 0), col = "black", col.axis = "black")
mtext(2, text = expression(paste(sigma^18, "", "(z-score)")), line = 2.1, cex = 1, col = "black")
abline(v = seq(6500, 2200, -250), lty = "dotted", col = "gray", lwd = 2)
text(x = 6300, y = -1.5, labels = "Jeita Cave (7)", font = 2, col = "darkgreen", cex = 2, adj = c(0, 0))
plot_bb_blocks(logistic_model)

# make sure we are only plotting existing values
bp <- paleoclim_long$bp[!is.na(paleoclim_long$z_score_lake_hula)]
z_score <- paleoclim_long$z_score_lake_hula[!is.na(paleoclim_long$z_score_lake_hula)]
ylim <- c(max(z_score, na.rm = TRUE), min(z_score, na.rm = TRUE))

par(mar = c(0, 4, 0, 3)) # c(bottom, left, top, right)
plot(bp, z_score, type = "l", axes = FALSE, xlim = c(6500, 2200), ylim = ylim, lwd = 1.5, col = "aquamarine4", xaxt = "n", ylab = "", xaxs = "i", yaxs = "i")
axis(4, at = seq(1.5, -1, -0.5), labels = seq(1.5, -1, -0.5), lwd = 1, line = 0.2, cex.axis = 1, mgp = c(0, 0.4, 0), col = "black", col.axis = "black")
mtext(4, text = expression(paste(sigma^18, "", "(z-score)")), line = 2.1, cex = 1, col = "black")
abline(v = seq(6500, 2200, -250), lty = "dotted", col = "gray", lwd = 2)
text(x = 4800, y = -1, labels = "Lake Hula (8)", font = 2, col = "aquamarine4", cex = 2, adj = c(0, 0))
plot_bb_blocks(logistic_model)

# make sure we are only plotting existing values
bp <- paleoclim_long$bp[!is.na(paleoclim_long$z_score_soreq_cave)]
z_score <- paleoclim_long$z_score_soreq_cave[!is.na(paleoclim_long$z_score_soreq_cave)]
ylim <- c(max(z_score, na.rm = TRUE), min(z_score, na.rm = TRUE))

par(mar = c(0, 4, 0, 3)) # c(bottom, left, top, right)
plot(bp, z_score, type = "l", axes = FALSE, xlim = c(6500, 2200), ylim = ylim, lwd = 1.5, col = "skyblue2", xaxt = "n", ylab = "", xaxs = "i", yaxs = "i", xlab = "")
axis(2, at = seq(3, -2, -0.5), labels = seq(3, -2, -0.5), lwd = 1, line = 0.2, cex.axis = 1, mgp = c(0, 0.5, 0), col = "black", col.axis = "black")
mtext(2, text = expression(paste(sigma^18, "", "(z-score)")), line = 2.1, cex = 1, col = "black")
abline(v = seq(6500, 2200, -250), lty = "dotted", col = "gray", lwd = 2)
text(x = 6300, y = -0.8, labels = "Soreq cave (10)", col = "skyblue2", font = 2, cex = 2, adj = c(0, 0))
plot_bb_blocks(logistic_model)

# make sure we are only plotting existing values
bp <- paleoclim_long$bp[!is.na(paleoclim_long$z_score_jerusalem_west_cave)]
z_score <- paleoclim_long$z_score_jerusalem_west_cave[!is.na(paleoclim_long$z_score_jerusalem_west_cave)]
ylim <- c(max(z_score, na.rm = TRUE), min(z_score, na.rm = TRUE))

par(mar = c(0, 4, 0, 3)) # c(bottom, left, top, right)
plot(bp, z_score, type = "l", axes = FALSE, xlim = c(6500, 2200), ylim = ylim, lwd = 1.5, col = "brown", xaxt = "n", ylab = "", xaxs = "i", yaxs = "i", xlab = "")
axis(4, at = seq(2, -1, -0.5), labels = seq(2, -1, -0.5), lwd = 1, line = 0.2, cex.axis = 1, mgp = c(0, 0.5, 0), col = "black", col.axis = "black")
mtext(4, text = expression(paste(sigma^18, "", "(z-score)")), line = 2.1, cex = 1, col = "black")
abline(v = seq(6500, 2200, -250), lty = "dotted", col = "gray", lwd = 2)
text(x = 6000, y = 0.7, labels = "Jerusalem W. cave (9)", col = "brown", font = 2, cex = 2, adj = c(0, 0))
plot_bb_blocks(logistic_model)

xticks <- seq(6500, 2200, -250)

par(mar = c(8, 4, 0, 3)) # c(bottom, left, top, right)
ymax <- max(logistic_model$result$PrDens) * 1.1
plot(logistic_model, ylim = c(0, ymax), xlim = c(6500, 2200), drawaxes = FALSE)
lines(logistic_model$fit$calBP, logistic_model$fit$PrDens, col = "black", lty = "dashed", lwd = 0.5)
abline(v = seq(6500, 2200, -250), lty = "dotted", col = "gray", lwd = 2)
legend(
  x = 6300, ymax * 0.94,
  legend = c("SPD (dates not normalised)", "Logistic Model", "95% MC envelope", "positive deviation", "negative deviation"),
  col = c(1, "black", "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)),
  lty = c(1, 2, 1, 1, 1), lwd = c(0.5, 0.5, 5, 5, 5), cex = 1.2, bg = "white", title = ""
)
text(x = 6200, ymax * 0.92, labels = "Logistic Fit", font = 2, cex = 1.4, adj = c(0, 0.7))
text(x = 6100, ymax * 0.52, cex = 1.4, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(logistic_model$pval, 4))))

text(x = 10400, y= ymax * 0.52, labels = "10.3k \nevent", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
text(x = 9400, y= ymax * 0.52, labels = "9.3k \nevent", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
text(x = 8400, y= ymax * 0.52, labels = "8.2k \nevent", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
text(x = 5250, y= ymax * 0.52, labels = "5.2k \nevent", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
text(x = 4250, y= ymax * 0.52, labels = "4.2k \nevent", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
text(x = 3250, y= ymax * 0.52, labels = "3.2k \nevent", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
text(x = 12800, y= ymax * 0.52, labels = "Younger \nDryas", font = 2, col = "black", cex = 1.5, adj = c(0, 0))
xticks <- seq(6500, 2200, -250)
axis(side = 1, at = xticks, labels = xticks, las = 2, cex.axis = 1.2) # add BP axis
mtext("cal BP", 1, 2.3, at = 6375, adj = 0.4, font = 2, cex = 0.8, las = 2)
axis(side = 1, at = xticks - 50, labels = xticks - 2000, las = 2, cex.axis = 1.2, pos = -0.00009) # add BC/AD axis
mtext("BC", 1, 5.3, at = 6375, adj = 0.4, font = 2, cex = 0.8, las = 2)
dev.off()
