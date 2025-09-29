## Perform aoristic analyses and monte carlo simulation on multiple archaeological proxy data
## Last edited 2025-09-22
## Authors: Andrea Titolo (refactored code), Alessio Palmisano (original code)
## Packages version at the time of analysis
## here: 1.0.1
## sf: 1.0-21

# Set-up -----------------------------------------------------------------------

# Load required libraries
required_packages <- c("sf", "here")

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

# Source functions scripts, to be expanded to simplify the script below
source(here::here("code/functions/calculate_bp_dates.R"))
source(here::here("code/functions/timeblocks_functions.R"))
source(here::here("code/functions/simulation_functions.R"))
source(here::here("code/functions/summarise_proxies.R"))
source(here::here("code/functions/generate_confidence_intervals.R"))
source(here::here("code/functions/normalization_functions.R"))
source(here::here("code/functions/compute_multiproxy_correlation.R"))

# Load data
# Download the dataset from github
# Check if file exists before downloading
if (!file.exists(here::here("data/dataset.gpkg"))) {
  download.file(
    url = "https://github.com/UnitoAssyrianGovernance/villages-to-empire-dataset/raw/refs/heads/main/dataset/dataset.gpkg",
    destfile = here::here("data/dataset.gpkg")
  )
} else {
  print("Dataset file already exists, skipping download...")
}

sites <- terra::vect(x = here("data/dataset.gpkg"), layer = "archaeo_sites")
sites <- terra::as.data.frame(sites)

# Retrieve the oldest and more recent date
min(sites$StartDate)
max(sites$EndDate)

# Keep only a handful of colums
# Ideally, in the future these data will be downloaded directly from the public github repo
sites <- sites[, c(
  "Name",
  "SizeHa",
  "Region",
  "StartDate",
  "EndDate",
  "SiteID"
)]

sites <- get_bp_dates(
  df = sites,
  start_date = sites$StartDate,
  end_date = sites$EndDate,
  duration = TRUE
)

# Retrieve the oldest and more recent date (BP)
max(sites$StartBP)
min(sites$EndBP)

# Create an histogram of sites duration and save it as a png
png(
  file = here::here("output/png/site_duration_histogram.png"),
  width = 1200,
  height = 800,
  res = 150
)

hist(
  sites$Duration, cex.axis = 0.75, xlab = "Site Date Range in Years",
  ylab = "Site Count", main = "", breaks = 30, col = "black"
)

dev.off()

# Variables -----------------------------------------------------------------------

# time blocks length
block_width <- 100
# block_width <- 200

# Simulation-related variables, all of these are arbitrary and need to be adjusted to the observed data
# min_duration is the minimum allowed duration for sites
# sd_duration is the standard deviation in years for the normal distribution
n_simulation <- 100
min_duration <- 10
sd_duration <- 50

# window of analyses
# you might want to cut recent periods, but highest date should always be relative to the source data
# 150 is to add a bit of buffer to the highest date and compensate for the 50 years added by the BP calculation
# the buffer is mostly useful for the plot
time_window <- c(start = 2000, end = max(sites$StartBP) + 150)

# Set filename based on block width
if (block_width == 200) {
  filename <- "200"
} else if (block_width == 100) {
  filename <- "100"
}

# Use normalized spd or not
use_normalized_spd <- TRUE

# if true, set up namings for output data, if not, empty string
if (use_normalized_spd) {
  norm_string <- "_norm_"
} else {
  norm_string <- ""
}

# Data preparation -----------------------------------------------------------------------

# Decide wether to load the dataset already produced by the analyses in this script
# This is useful to avoid repeating the analyses when reproducing only the plots
# NOTE: this does not influence the radiocarbon data, which are necessary for the plotting and are always loaded
load_data <- TRUE

# filter sites that do not span the time window generated above the StartBP and EndBP column
sites <- sites[
  sites$StartBP <= time_window[["end"]] & sites$EndBP >= time_window[["start"]],
]

# Get some quantitative info
# Site phases
nrow(sites)

# Unique sites
length(unique(sites$SiteID))

# hist(sites$Duration, cex.axis = 0.75, xlab = "Site Date Range in Years", ylab = "Site Count", main = "", breaks = 30, col = "black")

# Isolate sites for later sections of the script
sites_samaria <- split(sites, sites$Region)$Samaria
sites_judah <- split(sites, sites$Region)$Judah

# Load the results of the radiocarbon analyses
load(here::here("output/rda/radiocarbon_data.RData")) # Load the results of the Radiocarbon analyses

# Select appropriate spd data from the load above based on the variable use_normalized_spd
if (use_normalized_spd) {
  spd_data <- "all_spd_norm"
  spd_data_sam <- "all_spd_norm_samaria"
  spd_data_jud <- "all_spd_norm_judah"
} else {
  spd_data <- "allspd"
  spd_data_sam <- "all_spd_samaria"
  spd_data_jud <- "all_spd_judah"
}

# Set filename based on block width and load data if required
filename <- ifelse(block_width == 200, 200, ifelse(block_width == 100, 100, NA))

if (load_data && !use_normalized_spd) {
  data_files <- c("samaria_judah_", "samaria_", "judah_")
  for (x in data_files) {
    load(here::here(paste0("output/rda/", x, filename, ".RData")))
  }
  cat("Data files found and loaded")
} else if (load_data && use_normalized_spd) {
  data_files <- c("samaria_judah_norm_", "samaria_norm_", "judah_norm_")
  for (x in data_files) {
    load(here::here(paste0("output/rda/", x, filename, ".RData")))
  }
cat("Data files found and loaded")
} else if (!load_data) {
  cat("Set to not load any data, skipping...")
}

# SPDs are likely larger than the clipped dataset
# we can clip them here using again the time_window created above
spd_range <- time_window[["start"]]:time_window[["end"]]

# filter the all_spd_norm column values for cells containing values from the spd_range above
# get() is used to retrieve the global environment data from the spd_data variable
spd_filtered <- get(spd_data)$grid[get(spd_data)$grid$calBP %in% spd_range, ]

# Aoristic (whole area) -----------------------------------------------------------------------

## Create time blocks (100 yrs) and new columns
time_blocks <- gen_time_blocks(time_window, block_width)$time_blocks
new_columns <- gen_time_blocks(time_window, block_width)$new_columns

## Create new tables for storing results
tables_list <- gen_new_tables(sites, new_columns)

## Loop through and update timeblocks for each site
archaeo_proxies_l <- update_time_blocks(
  site_start = sites$StartBP,
  site_end = sites$EndBP,
  time_blocks = time_blocks,
  block_width = block_width,
  site_counts = tables_list[["site_counts"]],
  site_area_sums = tables_list[["site_area_sums"]],
  aorist_weights = tables_list[["aorist_weights"]],
  area_field = sites$SizeHa
)

# Randomised Start Dates and Duration -----------------------------------------------------------------------

## Randomised start date of site occupation phase (uniform). A sites duration randomly generated from a normal distribution with a mean of 200 years and a standard deviations of 50 years is added to the drawn date

## Create time blocks (100 yrs) and new columns
sim_list <- montecarlo_sim(
  sites = sites,
  nsim = n_simulation, # number of simulations for generating a 95% confidence envelope
  mean_duration = block_width, # mean of randomised site duration == block_width
  min_duration = min_duration, # arbitrary value
  sd_duration = sd_duration, # standard deviation, arbitrary value
  time_blocks = time_blocks,
  block_width = block_width,
  new_columns = new_columns
)

site_dur <- sim_list$site_dur
sim_df <- sim_list$sim_df

# Summarise data -----------------------------------------------------------------------

archaeo_proxies_sum_l <- summarise_proxies(
  site_counts = archaeo_proxies_l[["site_counts"]],
  aorist_weights = archaeo_proxies_l[["aorist_weights"]],
  site_area_sums = archaeo_proxies_l[["site_area_sums"]],
  simulated_df = sim_df
)

conf_intervals_l <- get_conf_intervals(sim_df)

# Normalization -----------------------------------------------------------------------

## Normalise all proxies between 6400-2200 Bp (ca. 4450 - 250 BC)
archaeo_proxies_norm <- normalize_archaeo_proxies(
  count_proxy = archaeo_proxies_sum_l[["site_blocksum_counts"]],
  area_proxy = archaeo_proxies_sum_l[["site_blocksum_area"]],
  weight_proxy = archaeo_proxies_sum_l[["site_blocksum_weights"]]
)

## SPD of radiocarbon dates
spd_norm <- normalize_spd(spd = spd_filtered)

simulated_df_norm <- normalize_sim_df(sim_df = sim_df)

# Get confidence intervals of normalized simulated_df and normalize them
conf_intervals_norm <- as.data.frame(lapply(
  get_conf_intervals(simulated_df_norm),
  normalize_vector
))

# Correlation -----------------------------------------------------------------------

# Calculate Pearson's correlation values between archaeo-demographic proxies

proxies_cor <- calculate_correlation(
  archaeo_proxies = archaeo_proxies_sum_l,
  spd = spd_norm,
  max_bp = 2800,
  row_names = c("Counts", "Aoristic Weight", "Area", "Simulated", "SPD")
)

proxies_cor

# Save correlation to csv
write.csv(
  proxies_cor,
  here::here(paste0(
    "output/csv/proxies_correlation",
    norm_string,
    filename,
    ".csv"
  ))
)

# Samaria -----------------------------------------------------------------------

# filter the all_spd_norm column values for cells containing values from the spd_range above
spd_filtered_sam <- get(spd_data_sam)$grid[
  get(spd_data_sam)$grid$calBP %in% spd_range,
]

# Aoristic (Samaria) -----------------------------------------------------------------------

## Create new tables for storing results
sam_tables_list <- gen_new_tables(sites_samaria, new_columns)

## Loop through and update timeblocks for each site
sam_archaeo_proxies_l <- update_time_blocks(
  site_start = sites_samaria$StartBP,
  site_end = sites_samaria$EndBP,
  time_blocks = time_blocks,
  block_width = block_width,
  site_counts = sam_tables_list[["site_counts"]],
  site_area_sums = sam_tables_list[["site_area_sums"]],
  aorist_weights = sam_tables_list[["aorist_weights"]],
  area_field = sites_samaria$SizeHa
)

# Randomised Start Dates and Duration -----------------------------------------------------------------------

## Randomised start date of site occupation phase (uniform). A sites duration randomly generated from a normal distribution with a mean of 200 years and a standard deviations of 50 years is added to the drawn date

## Create time blocks (100 yrs) and new columns
sam_sim_list <- montecarlo_sim(
  sites = sites_samaria,
  nsim = n_simulation, # number of simulations for generating a 95% confidence envelope
  mean_duration = block_width, # mean of randomised site duration == block_width
  min_duration = min_duration, # arbitrary value
  sd_duration = sd_duration, # standard deviation, arbitrary value
  time_blocks = time_blocks,
  block_width = block_width,
  new_columns = new_columns
)

sam_site_dur <- sam_sim_list$site_dur
sam_sim_df <- sam_sim_list$sim_df

# Summarise data -----------------------------------------------------------------------

sam_archaeo_proxies_sum_l <- summarise_proxies(
  site_counts = sam_archaeo_proxies_l[["site_counts"]],
  aorist_weights = sam_archaeo_proxies_l[["aorist_weights"]],
  site_area_sums = sam_archaeo_proxies_l[["site_area_sums"]],
  simulated_df = sam_sim_df
)

sam_conf_intervals_l <- get_conf_intervals(sam_sim_df)


# Normalization -----------------------------------------------------------------------
## Normalise all proxies between 6400-2200 Bp (ca. 4450 - 250 BC)
sam_archaeo_proxies_norm <- normalize_archaeo_proxies(
  count_proxy = sam_archaeo_proxies_sum_l[["site_blocksum_counts"]],
  area_proxy = sam_archaeo_proxies_sum_l[["site_blocksum_area"]],
  weight_proxy = sam_archaeo_proxies_sum_l[["site_blocksum_weights"]]
)

## SPD of radiocarbon dates
sam_spd_norm <- normalize_spd(spd = spd_filtered_sam)

sam_simulated_df_norm <- normalize_sim_df(sim_df = sam_sim_df)

# Get confidence intervals of normalized simulated_df and normalize them
sam_conf_intervals_norm <- as.data.frame(lapply(
  get_conf_intervals(sam_simulated_df_norm),
  normalize_vector
))

# Correlation -----------------------------------------------------------------------

# Calculate Pearson's correlation values between archaeo-demographic proxies

sam_proxies_cor <- calculate_correlation(
  archaeo_proxies = sam_archaeo_proxies_sum_l,
  spd = sam_spd_norm,
  max_bp = 2800,
  row_names = c("Counts", "Aoristic Weight", "Area", "Simulated", "SPD")
)

sam_proxies_cor

# Save correlation to csv
write.csv(
  sam_proxies_cor,
  here::here(paste0(
    "output/csv/sam_proxies_correlation",
    norm_string,
    filename,
    ".csv"
  ))
)

# Judah -----------------------------------------------------------------------

# filter the all_spd_norm column values for cells containing values from the spd_range above
spd_filtered_jud <- get(spd_data_jud)$grid[get(spd_data_jud)$grid$calBP %in% spd_range, ]

# Aoristic (Judah) -----------------------------------------------------------------------

## Create new tables for storing results
jud_tables_list <- gen_new_tables(sites_judah, new_columns)

## Loop through and update timeblocks for each site
jud_archaeo_proxies_l <- update_time_blocks(
  site_start = sites_judah$StartBP, site_end = sites_judah$EndBP,
  time_blocks = time_blocks, block_width = block_width, site_counts = jud_tables_list[["site_counts"]],
  site_area_sums = jud_tables_list[["site_area_sums"]], aorist_weights = jud_tables_list[["aorist_weights"]],
  area_field = sites_judah$SizeHa
)

# Randomised Start Dates and Duration -----------------------------------------------------------------------

## Randomised start date of site occupation phase (uniform). A sites duration randomly generated from a normal distribution with a mean of 200 years and a standard deviations of 50 years is added to the drawn date

## Create time blocks (100 yrs) and new columns
jud_sim_list <- montecarlo_sim(
  sites = sites_judah,
  nsim = n_simulation, # number of simulations for generating a 95% confidence envelope
  mean_duration = block_width, # mean of randomised site duration == block_width
  min_duration = min_duration, # arbitrary value
  sd_duration = sd_duration, # standard deviation, arbitrary value
  time_blocks = time_blocks,
  block_width = block_width,
  new_columns = new_columns
)

jud_site_dur <- jud_sim_list$site_dur
jud_sim_df <- jud_sim_list$sim_df

# Summarise data -----------------------------------------------------------------------

jud_archaeo_proxies_sum_l <- summarise_proxies(
  site_counts = jud_archaeo_proxies_l[["site_counts"]],
  aorist_weights = jud_archaeo_proxies_l[["aorist_weights"]],
  site_area_sums = jud_archaeo_proxies_l[["site_area_sums"]],
  simulated_df = jud_sim_df
)

jud_conf_intervals_l <- get_conf_intervals(jud_sim_df)


# Normalization -----------------------------------------------------------------------
## Normalise all proxies between 6400-2200 Bp (ca. 4450 - 250 BC)
jud_archaeo_proxies_norm <- normalize_archaeo_proxies(
  count_proxy = jud_archaeo_proxies_sum_l[["site_blocksum_counts"]],
  area_proxy = jud_archaeo_proxies_sum_l[["site_blocksum_area"]],
  weight_proxy = jud_archaeo_proxies_sum_l[["site_blocksum_weights"]]
)

## SPD of radiocarbon dates
jud_spd_norm <- normalize_spd(spd = spd_filtered_jud)

jud_simulated_df_norm <- normalize_sim_df(sim_df = jud_sim_df)

# Get confidence intervals of normalized simulated_df and normalize them
jud_conf_intervals_norm <- as.data.frame(lapply(get_conf_intervals(jud_simulated_df_norm), normalize_vector))


# Correlation -----------------------------------------------------------------------

# Calculate Pearson's correlation values between archaeo-demographic proxies

jud_proxies_cor <- calculate_correlation(
  archaeo_proxies = jud_archaeo_proxies_sum_l,
  spd = jud_spd_norm,
  max_bp = 2800,
  row_names = c("Counts", "Aoristic Weight", "Area", "Simulated", "SPD")
)

jud_proxies_cor

# Save correlation to csv
write.csv(
  jud_proxies_cor,
  here::here(paste0(
    "output/csv/jud_proxies_correlation",
    norm_string,
    filename,
    ".csv"
  ))
)

# Save data -----------------------------------------------------------------------

# Save data (whole area)
save(spd_range, spd_filtered, new_columns,
  time_blocks, archaeo_proxies_l, sim_list, site_dur,
  sim_df, archaeo_proxies_sum_l, conf_intervals_l,
  archaeo_proxies_norm, spd_norm, simulated_df_norm,
  conf_intervals_norm, proxies_cor,
  file = here::here(paste0("output/rda/samaria_judah", norm_string, filename, ".RData"))
)

# Save data (samaria)
save(spd_range, spd_filtered_sam, new_columns,
  time_blocks, sam_archaeo_proxies_l, sam_sim_list, sam_site_dur,
  sam_sim_df, sam_archaeo_proxies_sum_l, sam_conf_intervals_l,
  sam_archaeo_proxies_norm, sam_spd_norm, sam_simulated_df_norm,
  sam_conf_intervals_norm, sam_proxies_cor,
  file = here::here(paste0("output/rda/samaria", norm_string, filename, ".RData"))
)

# Save data (judah)
save(spd_range, spd_filtered_jud, new_columns,
  time_blocks, jud_archaeo_proxies_l, jud_sim_list, jud_site_dur,
  jud_sim_df, jud_archaeo_proxies_sum_l, jud_conf_intervals_l,
  jud_archaeo_proxies_norm, jud_spd_norm, jud_simulated_df_norm,
  jud_conf_intervals_norm, jud_proxies_cor,
  file = here::here(paste0("output/rda/judah", norm_string, filename, ".RData"))
)

# Plotting ---------------------------------------------------------------------

# Variables for the plots
# Calculate site numbers and phases for the time-window of interest (analysis window is larger to properly capture late IA trends)
sites_clip <- sites[sites$StartBP <= time_window[[2]] & sites$EndBP >= time_window[[1]] + 200, ]

# Site phases
site_phases_sam <- nrow(sites_clip[which(sites_clip$Region=="Samaria"),])
site_phases_jud <- nrow(sites_clip[which(sites_clip$Region=="Judah"),])
site_phases <- site_phases_sam + site_phases_jud

# Unique sites
n_sites_sam <- length(unique(sites_clip[which(sites_clip$Region=="Samaria"),]$SiteID))
n_sites_jud <- length(unique(sites_clip[which(sites_clip$Region=="Judah"),]$SiteID))
n_sites <- n_sites_sam + n_sites_jud

mids <- seq(min(time_blocks) + (block_width / 2), max(time_blocks) - (block_width / 2), by = block_width)
x_ticks <- seq(time_window[[2]], time_window[[1]], -200)

png(file = here::here("output/png/fig-06.png"), width = 1000, height = 1400, res = 150)
layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE), widths = 8, heights = c(2.2, 2.2, 2.5))
par(mar = c(1, 1, 0, 1))
par(yaxs = "i")
par(xaxs = "i")
xticks <- x_ticks
par(mar = c(0, 1, 0.5, 1)) # c(bottom, left, top, right)
plot(archaeo_proxies_norm$count ~ mids, lty = "solid", col = "black", cex.axis = 0.5, xaxt = "n", yaxt = "n", type = "l", xlab = "", ylab = "", xlim = c(time_window[[2]], time_window[[1]]+200))
lines(mids, archaeo_proxies_norm$weight, lty = "dashed", col = "black")
lines(mids, archaeo_proxies_norm$area, col = "red")
lines(spd_norm$x, spd_norm$y, col = "darkorange")
polygon(x = c(mids, rev(mids)), y = c(conf_intervals_norm$hicount, rev(conf_intervals_norm$lowcount)), col = rgb(128, 128, 128, alpha = 120, maxColorValue = 255), border = NA)
abline(v = seq(9000, 800, -100), lwd = 1.2, lty = "dotted", col = "grey80")
legend(x = 6600, y = 0.99, legend = c("raw count", "total area", "aoristic weight", "randomised start date", "SPD Samaria + Judah"), lty = c("solid", "solid", "dashed", "solid", "solid", "dashed"), lwd = c(1, 1, 1, 3), col = c("black", "red", "black", "grey75", "darkorange"), bty = "n", cex = 1, bg = "transparent", title = "")
text(x = 6550, y = 0.97, labels = "Samaria + Judah", font = 2, cex = 1.3, adj = c(0, 0.7))
text(x = 5500, y = 0.96, labels = paste("Site-phases = ", site_phases, " Sites = ", n_sites, sep = ""), font = 1, cex = 1, adj = c(0, 0.7))
text(x = 5500, y = 0.91, labels = paste("Radiocarbon data, ", "n = ", nrow(rc_dates), ", Sites = ", length(unique(rc_dates$SiteID)), ", bins = ", length(unique(bins)), sep = ""), font = 1, cex = 1, col = "chocolate", adj = c(0, 0.7))
text(x = 6000, y = 0.57, labels = "6000", font = 1, cex = 1, adj = c(0, 0.5), srt = 90)
text(x = 5000, y = 0.57, labels = "5000", font = 1, cex = 1, adj = c(0, 0.5), srt = 90)
text(x = 4000, y = 0.57, labels = "4000", font = 1, cex = 1, adj = c(0, 0.5), srt = 90)
text(x = 3000, y = 0.64, labels = "3000", font = 1, cex = 1, adj = c(0, 0.5), srt = 90)

par(mar = c(0.5, 1, 0.5, 1)) # c(bottom, left, top, right)
plot(sam_archaeo_proxies_norm$count ~ mids, lty = "solid", col = "black", cex.axis = 0.5, xaxt = "n", yaxt = "n", type = "l", xlab = "", ylab = "", xlim = c(time_window[[2]], time_window[[1]]+200))
lines(mids, sam_archaeo_proxies_norm$weight, lty = "dashed", col = "black")
lines(mids, sam_archaeo_proxies_norm$area, col = "red")
lines(sam_spd_norm$x, sam_spd_norm$y, col = "darkorange")
polygon(x = c(mids, rev(mids)), y = c(sam_conf_intervals_norm$hicount, rev(sam_conf_intervals_norm$lowcount)), col = rgb(128, 128, 128, alpha = 120, maxColorValue = 255), border = NA)
abline(v = seq(9000, 800, -100), lwd = 1.2, lty = "dotted", col = "grey80")
legend(x = 6600, y = 0.99, legend = c("raw count", "total area", "aoristic weight", "randomised start date", "SPD Samaria"), lty = c("solid", "solid", "dashed", "solid", "solid", "dashed"), lwd = c(1, 1, 1, 3), col = c("black", "red", "black", "grey75", "darkorange"), bty = "n", cex = 1, bg = "transparent", title = "")
text(x = 6550, y = 0.97, labels = "Samaria", font = 2, cex = 1.3, adj = c(0, 0.7))
text(x = 5500, y = 0.96, labels = paste("Site-phases = ", site_phases_sam, " Sites = ", n_sites_sam, sep = ""), font = 1, cex = 1, adj = c(0, 0.7))
text(x = 5500, y = 0.91, labels = paste("Radiocarbon data, ", "n = ", nrow(rc_dates_samaria), ", Sites = ", length(unique(rc_dates_samaria$SiteID)), ", bins = ", length(unique(bins_samaria)), sep = ""), font = 1, cex = 1, col = "chocolate", adj = c(0, 0.7))

par(mar = c(7, 1, 0, 1)) # c(bottom, left, top, right)
plot(jud_archaeo_proxies_norm$count ~ mids, lty = "solid", col = "black", cex.axis = 0.5, xaxt = "n", yaxt = "n", type = "l", xlab = "", ylab = "", xlim = c(time_window[[2]], time_window[[1]]+200))
lines(mids, jud_archaeo_proxies_norm$weight, lty = "dashed", col = "black")
lines(mids, jud_archaeo_proxies_norm$area, col = "red")
lines(jud_spd_norm$x, jud_spd_norm$y, col = "darkorange")
polygon(x = c(mids, rev(mids)), y = c(jud_conf_intervals_norm$hicount, rev(jud_conf_intervals_norm$lowcount)), col = rgb(128, 128, 128, alpha = 120, maxColorValue = 255), border = NA)
abline(v = seq(9000, 800, -100), lwd = 1.2, lty = "dotted", col = "grey80")
legend(x = 6600, y = 0.99, legend = c("raw count", "total area", "aoristic weight", "randomised start date", "SPD Judah"), lty = c("solid", "solid", "dashed", "solid", "solid", "dashed"), lwd = c(1, 1, 1, 3), col = c("black", "red", "black", "grey75", "darkorange"), bty = "n", cex = 1, bg = "transparent", title = "")
text(x = 6550, y = 0.97, labels = "Judah", font = 2, cex = 1.3, adj = c(0, 0.7))
text(x = 5500, y = 0.96, labels = paste("Site-phases = ", site_phases_jud, " Sites = ", n_sites_jud, sep = ""), font = 1, cex = 1, adj = c(0, 0.7))
text(x = 5500, y = 0.91, labels = paste("Radiocarbon data, ", "n = ", nrow(rc_dates_judah), ", Sites = ", length(unique(rc_dates_judah$SiteID)), ", bins = ", length(unique(bins_judah)), sep = ""), font = 1, cex = 1, col = "chocolate", adj = c(0, 0.7))
axis(side = 1, at = xticks, labels = xticks, las = 2, cex.axis = 1)
mtext("cal BP", 1, 3, at = 6500, adj = -0.05, srt = 90, font = 2, cex = 0.7, las = 2)
axis(side = 1, at = xticks - 50, labels = xticks - 2000, las = 2, cex.axis = 1, pos = -0.20) # add BC/AD axis
mtext("BC/AD", 1, 6.7, at = 6450, adj = 0, font = 2, cex = 0.7, las = 2)
dev.off()
