## Calibrate radiocarbon dates and generate SPDs for the whole project area and the two regions separately
## Last edited 2025-09-22
## Authors: Andrea Titolo (refactored code), Alessio Palmisano (original code)
## Packages version at the time of analysis
## doParallel: 1.0.17
## rcarbon: 1.5.2
## here: 1.0.1
## sf: 1.0-21

# Set-up -----------------------------------------------------------------------

# Decide weather to load the RData file with the outputs of this analyses or not.
load_data <- TRUE

# Load required libraries
required_packages <- c("here", "rcarbon", "doParallel", "terra")

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

# Load radiocarbon data
rc_dates <- terra::vect(
  x = here::here("data/processed/nerd_xronos.gpkg"),
  layer = "nerd_xronos"
)

# Turn SpatVector to a data.frame
rc_dates <- terra::as.data.frame(rc_dates)

# The generation of some SPDs and null models could be time consuming.
# So, you can load the RData file storing all the outputs generated via the R script here provided.
if (file.exists(here::here("output/rda/radiocarbon_data.RData")) && load_data) {
  load(here::here("output/rda/radiocarbon_data.RData"))
  print("Radiocarbon RData file found and loaded")
} else if (!file.exists(here::here("output/rda/radiocarbon_data.RData" && load_data))) {
  browseURL(
    url = "https://ndownloader.figshare.com/files/58302040?private_link=61729c5d2acfb1f2d6a8",
  )
  cat(
    "==========\nBrowser download page opened\n",
    "You should place this file in the 'output/rda' directory\n",
    "Please re-run the script to load the file\n=========="
  )
} else {
  cat(
    "RData file not found \nPlease run the rest of the script and save the data\n",
    "Or if you want to download from Figshare, change the variable 'load_data' to TRUE"
  )
}

# Data cleaning -----------------------------------------------------------------------

## Drop dates with NA values in the column CRA, with time range wider then needed, and with the value "shell" in the column Material
rc_dates <- rc_dates[!is.na(rc_dates$CRA) & rc_dates$CRA < 15730, ]
rc_dates <- rc_dates[rc_dates$Material != "shell" &
  rc_dates$Material != "shell (marine)" | is.na(rc_dates$Material),
]
rc_dates <- rc_dates[
  rc_dates$CRA != 3955 & rc_dates$CRA != 3800 & rc_dates$CRA != 3740 &
    rc_dates$CRA != 3740 & rc_dates$CRA != 3704 & rc_dates$CRA != 3700 &
    rc_dates$CRA != 3635 | !is.na(rc_dates$LabID) & rc_dates$LabID != "OZC-173",
]

# Query the number of sites
length(unique(rc_dates$SiteName))

# Isolate again samaria and judah dates after filtering above
rc_dates_samaria <- rc_dates[rc_dates$Region == "Samaria", ]
rc_dates_judah <- rc_dates[rc_dates$Region == "Judah", ]

# General SPD parameters -----------------------------------------------------------------------

# Calculate numbers of cores available on the machine
cores_tot <- parallel::detectCores()

n_simulation <- 1000 # number of actual simulations
n_cores <- cores_tot / 2 # multi-core processing (use half the number of total cores)
smoothing_window <- 50 # smoothing of SPDs
bin_cluster <- 100 # bin clustering
real_start_bp <- 6450 # Actual starting date in years Before Present
real_end_bp <- 1950 # Actual ending date in years Before Present
buffer_years <- 1000 # Buffer time in years to extend analysis range
working_start_bp <- real_start_bp + buffer_years # Extended starting date for analysis
working_end_bp <- real_end_bp - buffer_years # Extended ending date for analysis
if (working_end_bp < 0) {
  working_end_bp <- 0
}

# Processing -----------------------------------------------------------------------

## Calibration -----------------------------------------------------------------------

# Calibrate dates (unnormalised for now)
all_dates <- rcarbon::calibrate(
  x = rc_dates$CRA,
  errors = rc_dates$Error,
  calCurves = "intcal20",
  method = "standard",
  normalised = FALSE,
  ncores = n_cores,
  calMatrix = TRUE
)

all_dates_norm <- rcarbon::calibrate(
  x = rc_dates$CRA,
  errors = rc_dates$Error,
  calCurves = "intcal20",
  method = "standard",
  normalised = TRUE,
  ncores = n_cores,
  calMatrix = TRUE
)

samaria_dates <- rcarbon::calibrate(
  x = rc_dates_samaria$CRA,
  errors = rc_dates_samaria$Error,
  calCurves = "intcal20",
  method = "standard",
  normalised = FALSE,
  ncores = n_cores,
  calMatrix = TRUE
)

samaria_dates_norm <- rcarbon::calibrate(
  x = rc_dates_samaria$CRA,
  errors = rc_dates_samaria$Error,
  calCurves = "intcal20",
  method = "standard",
  normalised = TRUE,
  ncores = n_cores,
  calMatrix = TRUE
)

judah_dates <- rcarbon::calibrate(
  x = rc_dates_judah$CRA,
  errors = rc_dates_judah$Error,
  calCurves = "intcal20",
  method = "standard",
  normalised = FALSE,
  ncores = n_cores,
  calMatrix = TRUE
)

judah_dates_norm <- rcarbon::calibrate(
  x = rc_dates_judah$CRA,
  errors = rc_dates_judah$Error,
  calCurves = "intcal20",
  method = "standard",
  normalised = TRUE,
  ncores = n_cores,
  calMatrix = TRUE
)

# Inspect the distribution calibrated dates if necessary
plot(all_dates, 2, calendar = "BCAD")

## Summed Probability Distributions (SPDs) - whole region -----------------------------------------------------------------------

# Generate bins
bins <- rcarbon::binPrep(sites = rc_dates$SiteID, ages = rc_dates$CRA, h = bin_cluster)

all_spd_norm <- rcarbon::spd(
  x = all_dates_norm,
  bins = bins,
  timeRange = c(working_start_bp, working_end_bp),
  datenormalised = TRUE, spdnormalised = TRUE,
  runm = smoothing_window,
  edgeSize = 400)

all_spd <- rcarbon::spd(
  x = all_dates,
  bins = bins,
  timeRange = c(working_start_bp, working_end_bp),
  datenormalised = FALSE,
  runm = smoothing_window)

# Samaria
bins_samaria <- rcarbon::binPrep(
  sites = rc_dates_samaria$SiteID,
  ages = rc_dates_samaria$CRA,
  h = bin_cluster
)

all_spd_norm_samaria <- rcarbon::spd(
  x = samaria_dates_norm,
  bins = bins_samaria,
  timeRange = c(working_start_bp, working_end_bp),
  datenormalised = TRUE,
  spdnormalised = TRUE,
  runm = smoothing_window,
  edgeSize = 400
)

all_spd_samaria <- rcarbon::spd(
  x = samaria_dates,
  bins = bins_samaria,
  timeRange = c(working_start_bp, working_end_bp),
  datenormalised = FALSE,
  runm = smoothing_window
)

# Judah
bins_judah <- rcarbon::binPrep(
  sites = rc_dates_judah$SiteID,
  ages = rc_dates_judah$CRA,
  h = bin_cluster
)

all_spd_norm_judah <- rcarbon::spd(
  x = judah_dates_norm,
  bins = bins_judah,
  timeRange = c(working_start_bp, working_end_bp),
  datenormalised = TRUE,
  spdnormalised = TRUE,
  runm = smoothing_window,
  edgeSize = 400
)
all_spd_judah <- rcarbon::spd(
  x = judah_dates,
  bins = bins_judah,
  timeRange = c(working_start_bp, working_end_bp),
  datenormalised = FALSE,
  runm = smoothing_window
)

# Calculate the Pearson correlation coefficient between the unnormalised and normalised SPDS obtained calibrating all radiocarbon samples
cor.test(all_spd$grid$PrDens, all_spd_norm$grid$PrDens, method = "pearson")

# Do the same for samaria
cor.test(
  all_spd_samaria$grid$PrDens,
  all_spd_norm_samaria$grid$PrDens,
  method = "pearson"
)

# And Judah
cor.test(
  all_spd_judah$grid$PrDens,
  all_spd_norm_judah$grid$PrDens,
  method = "pearson"
)

# calculate the median of bins
bins_median <- rcarbon::binMed(x = all_dates_norm, bins = bins)

# calculate the median of bins for Samaria
bins_median_samaria <- rcarbon::binMed(x = samaria_dates_norm, bins = bins_samaria)

# calculate the median of bins for Judah
bins_median_judah <- rcarbon::binMed(x = judah_dates_norm, bins = bins_judah)

# Exponential model -----------------------------------------------------------------------

set.seed(123)
exp_null_model <- rcarbon::modelTest(
  all_dates_norm,
  errors = rc_dates$Error,
  bins = bins,
  nsim = n_simulation,
  runm = smoothing_window,
  timeRange = c(working_start_bp, working_end_bp),
  model = "exponential",
  method = "calsample",
  ncores = n_cores,
  datenormalised = TRUE,
  spdnormalised = TRUE,
  edgeSize = 400
)

## Logistic model -----------------------------------------------------------------------

set.seed(123)

# Generate a smoothed SPD
spd_smoothed <- rcarbon::spd(
  all_dates_norm,
  timeRange = c(working_start_bp, working_end_bp),
  bins = bins,
  runm = smoothing_window,
  datenormalised = TRUE,
  spdnormalised = TRUE,
  edgeSize = 400
)

# Calculate starting values programmatically
spd_data <- spd_smoothed$grid

# Find asymptote (maximum PrDens) with a 10% buffer, and round it
asym_start <- as.numeric(
  format(signif(max(spd_data$PrDens) * 1.1, digits = 1),
  scientific = FALSE)
)

# Find midpoint (where growth rate changes)
growth_rate <- diff(spd_data$PrDens)/diff(spd_data$calBP)

# Point of maximum decline round to the nearest 100 years
x_mid_start <- round(spd_data$calBP[which.min(growth_rate)], digits = -2)

# Start values should be adjusted depending on the observed SPD
logistic_model_fit <- nls(
  PrDens ~ SSlogis(calBP, Asym, xmid, scale),
  data = spd_smoothed$grid,
  control = nls.control(maxiter = 200),
  start = list(Asym = asym_start, xmid = x_mid_start, scale = -100)
)

# Generate a data frame containing the fitted values
logistic_model_dens <- data.frame(
  calBP = spd_smoothed$grid$calBP,
  PrDens = SSlogis(
    input = spd_smoothed$grid$calBP,
    Asym = coefficients(logistic_model_fit)[1],
    xmid = coefficients(logistic_model_fit)[2],
    scal = coefficients(logistic_model_fit)[3]
  )
)

# Use the modelTest function (returning the raw simulation output - see below)
logistic_model <- rcarbon::modelTest(
  all_dates_norm,
  errors = rc_dates$Error,
  bins = bins,
  nsim = n_simulation,
  timeRange = c(working_start_bp, working_end_bp),
  model = "custom",
  predgrid = logistic_model_dens,
  runm = smoothing_window,
  raw = TRUE, method = "calsample",
  normalised = TRUE,
  spdnormalised = TRUE
)

## Composite Kernel Density Estimates (CKDE) -----------------------------------------------------------------------

# CKDE weighted to emulate an SPD with non-normalised dates
random_dates <- rcarbon::sampleDates(
  all_dates_norm,
  bins = bins,
  nsim = n_simulation,
  boot = TRUE
) # sample random calendar dates

ckde <- rcarbon::ckde(
  random_dates,
  timeRange = c(working_start_bp, working_end_bp),
  bw = 100,
  normalised = TRUE
)

## Regional SDPs -----------------------------------------------------------------------

### Regional SPD - Logistic Model -----------------------------------------------------------------------

#### Samaria -----------------------------------------------------------------------

# Model SPD for each region and compare versus a logistic model
# Fit a logistic model and create a predicted SPD for Samaria
spd_smoothed_sam <- rcarbon::spd(
  samaria_dates_norm,
  timeRange = c(working_start_bp, working_end_bp),
  bins = bins_samaria,
  runm = smoothing_window,
  datenormalised = TRUE,
  spdnormalised = TRUE,
  edgeSize = 400
)

# Calculate starting values programmatically
spd_data_sam <- spd_smoothed_sam$grid

# Find asymptote (maximum PrDens) with a 10% buffer, and round it
asym_start_sam <- format(
  signif(max(spd_data_sam$PrDens) * 1.1, digits = 1),
  scientific = FALSE
)

# Find midpoint (where growth rate changes)
growth_rate_sam <- diff(spd_data_sam$PrDens) / diff(spd_data_sam$calBP)

# Point of maximum decline round to the nearest 100 years
x_mid_start_sam <- round(
  spd_data_sam$calBP[which.min(growth_rate_sam)],
  digits = -2
)

# Start values should be adjusted depending on the observed SPD
logistic_model_fit_sam <- nls(
  PrDens ~ SSlogis(calBP, Asym, xmid, scale),
  data = spd_smoothed_sam$grid,
  control = nls.control(maxiter = 200),
  start = list(Asym = asym_start_sam, xmid = x_mid_start_sam, scale = -100)
)

# Generate a data frame containing the fitted values
logistic_model_dens_sam <- data.frame(
  calBP = spd_smoothed_sam$grid$calBP,
  PrDens = SSlogis(
    input = spd_smoothed_sam$grid$calBP,
    Asym = coefficients(logistic_model_fit_sam)[1],
    xmid = coefficients(logistic_model_fit_sam)[2],
    scal = coefficients(logistic_model_fit_sam)[3]
  )
)

# Use the modelTest function (returning the raw simulation output - see below)
logistic_model_samaria <- rcarbon::modelTest(
  samaria_dates_norm,
  errors = rc_dates_samaria$Error,
  bins = bins_samaria,
  nsim = n_simulation,
  timeRange = c(working_start_bp, working_end_bp),
  model = "custom",
  predgrid = logistic_model_dens_sam,
  runm = smoothing_window,
  raw = TRUE,
  method = "calsample",
  normalised = TRUE,
  spdnormalised = TRUE
)

#### Judah -----------------------------------------------------------------------

## Fit a logistic model and create a predicted SPD for Judah
spd_smoothed_jud <- rcarbon::spd(
  judah_dates_norm,
  timeRange = c(working_start_bp, working_end_bp),
  bins = bins_judah,
  runm = smoothing_window,
  datenormalised = TRUE,
  spdnormalised = TRUE,
  edgeSize = 400
)

# Calculate starting values programmatically
spd_data_jud <- spd_smoothed_jud$grid

# Find asymptote (maximum PrDens) with a 10% buffer, and round it
asym_start_jud <- format(
  signif(max(spd_data_jud$PrDens) * 1.1, digits = 1),
  scientific = FALSE
)

# Find midpoint (where growth rate changes)
growth_rate_jud <- diff(spd_data_jud$PrDens) / diff(spd_data_jud$calBP)

# Point of maximum decline round to the nearest 100 years
x_mid_start_jud <- round(
  spd_data_jud$calBP[which.min(growth_rate_jud)],
  digits = -2
)

# Start values should be adjusted depending on the observed SPD
logistic_model_fit_jud <- nls(
  PrDens ~ SSlogis(calBP, Asym, xmid, scale),
  data = spd_smoothed_jud$grid,
  control = nls.control(maxiter = 200),
  start = list(Asym = asym_start_jud, xmid = x_mid_start_jud, scale = -100)
)

# Generate a data frame containing the fitted values
logistic_model_dens <- data.frame(
  calBP = spd_smoothed_jud$grid$calBP,
  PrDens = SSlogis(
    input = spd_smoothed_jud$grid$calBP,
    Asym = coefficients(logistic_model_fit_jud)[1],
    xmid = coefficients(logistic_model_fit_jud)[2],
    scal = coefficients(logistic_model_fit_jud)[3]
  )
)

# Use the modelTest function (returning the raw simulation output - see below)
logistic_model_judah <- rcarbon::modelTest(
  judah_dates_norm,
  errors = rc_dates_judah$Error,
  bins = bins_judah,
  nsim = n_simulation,
  timeRange = c(working_start_bp, working_end_bp),
  model = "custom",
  predgrid = logistic_model_dens,
  runm = smoothing_window,
  raw = TRUE,
  method = "calsample",
  normalised = TRUE,
  spdnormalised = TRUE
)

### Regional SPD - CKDE -----------------------------------------------------------------------

#### Samaria
random_dates_samaria <- rcarbon::sampleDates(
  samaria_dates_norm,
  bins = bins_samaria,
  nsim = n_simulation,
  boot = TRUE
) # sample random calendar dates

ckde_samaria <- rcarbon::ckde(
  random_dates_samaria,
  timeRange = c(working_start_bp, working_end_bp),
  bw = 100,
  normalised = TRUE
)

#### Judah
random_dates_judah <- rcarbon::sampleDates(
  judah_dates_norm,
  bins = bins_judah,
  nsim = n_simulation,
  boot = TRUE
) # sample random calendar dates

ckde_judah <- rcarbon::ckde(
  random_dates_judah,
  timeRange = c(working_start_bp, working_end_bp),
  bw = 100,
  normalised = TRUE
)

### Regional SPD - Permutation test -----------------------------------------------------------------------

# Test for Regional Departure from Global Pattern
set.seed(123)
perm_test <- rcarbon::permTest(
  x = all_dates_norm,
  bins = bins,
  marks = rc_dates$Region,
  timeRange = c(working_start_bp, working_end_bp),
  runm = smoothing_window,
  nsim = n_simulation,
  datenormalised = TRUE,
  spdnormalised = TRUE
)

# Save data -----------------------------------------------------------------------

# Create directory if it doesn't exist for any reason
if (!dir.exists(here::here("output/rda"))) {
  dir.create(here::here("output/rda"))
} else {
  print("Directory exists, skipping creation...")
}

save(all_dates, all_dates_norm, all_spd, all_spd_norm,
  ckde, ckde, ckde_judah, ckde_judah, ckde_samaria, ckde_samaria,
  all_spd_judah, all_spd_samaria, all_spd_norm_samaria, all_spd_norm_judah,
  bins, bins_judah, bins_samaria, bins_median, bins_median_judah, bins_median_samaria,
  rc_dates, rc_dates_samaria, rc_dates_judah,
  exp_null_model, perm_test, logistic_model, logistic_model_judah, logistic_model_samaria,
  file = here::here("output/rda/radiocarbon_data.RData")
)

# Plotting -----------------------------------------------------------------------

## Plot SPDs (whole area 6500-2200 BP) -----------------------------------------------------------------------

png(file = here::here("output/png/fig-03.png"), width = 1200, height = 1000, res = 150)
layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE), widths = 6, heights = c(2.2, 2.2, 3))
par(mar = c(0, 1, 1, 1)) # c(bottom, left, top, right)

# Plot SPD and cKDE
ymax <- max(all_spd$grid$PrDens)
plot(all_spd_norm, xlim = c(6500, 2200), xaxt = "n", yaxt = "n", rescale = TRUE)
rcarbon::barCodes(bins_median)
abline(v = seq(6500, 2500, -250), lty = "dotted", col = "white")
text(x = 6400, ymax * 5.6, labels = "a. All Dates", font = 2, cex = 0.9, adj = c(0, 0.7))
text(x = 6400, ymax * 4.9, labels = paste("n=", nrow(rc_dates), ", sites=", length(unique(rc_dates$SiteID)), ", bins=", length(unique(bins)), sep = ""), font = 1, cex = 0.9, adj = c(0, 0.7))
legend(x = 6400, ymax * 4.85, legend = c("SPD normalised", "cKDE"), lty = c("solid", "solid"), lwd = c(3, 0.5), col = c("grey90", "lightslateblue"), bty = "n", cex = 0.9)
segments(2750, ymax * 5.9, 2350, ymax * 5.9, lty = "solid", col = "grey23", lwd = (1))
text(x = 2550, ymax * 5.91, labels = "Hallstat Plataeu", lty = c("solid"), col = "grey23", font = 2, cex = 0.7, adj = c(0.5, -0.5))
box()

# Plot CKDE on top
par(new = TRUE)
plot(ckde, xlim = c(6500, 2200), type = "envelope", line.col = "black", fill.col = (rgb(0, 0, 255, alpha = 60, maxColorValue = 255)), xaxt = "n", yaxt = "n")

# Plot Exponential Model
par(mar = c(0, 1, 0, 1)) # c(bottom, left, top, right)
ymax <- max(exp_null_model$result$PrDens) * 1.1
plot(exp_null_model, ylim = c(0, ymax), xlim = c(6500, 2200), drawaxes = FALSE)
lines(exp_null_model$fit$calBP, exp_null_model$fit$PrDens, col = "black", lty = "dashed", lwd = 0.5)
abline(v = seq(6500, 2500, -250), lty = "dotted", col = "white")
legend(
  x = 6400, ymax * 0.94,
  legend = c("SPD (dates normalised)", "Exponential Model", "95% MC envelope", "positive deviation", "negative deviation"),
  col = c(1, "black", "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)),
  lty = c(1, 2, 1, 1, 1), lwd = c(0.5, 0.5, 5, 5, 5), cex = 0.6, bg = "white", title = ""
)
text(x = 6350, ymax * 0.92, labels = "b. Exponential Fit", font = 2, cex = 0.8, adj = c(0, 0.7))
text(x = 6200, ymax * 0.52, cex = 0.7, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(exp_null_model$pval, 4))))
segments(2750, ymax * 0.9, 2350, ymax * 0.9, lty = "solid", col = "grey23", lwd = (1))
text(x = 2550, ymax * 0.91, labels = "Hallstat Plataeu", lty = c("solid"), col = "grey23", font = 2, cex = 0.7, adj = c(0.5, -0.5))
box()

# Plot Logistic Model
par(mar = c(7, 1, 0, 1)) # c(bottom, left, top, right)
ymax <- max(logistic_model$result$PrDens) * 1.1
plot(logistic_model, ylim = c(0, ymax), xlim = c(6500, 2200), drawaxes = FALSE)
lines(logistic_model$fit$calBP, logistic_model$fit$PrDens, col = "black", lty = "dashed", lwd = 0.5)
abline(v = seq(6500, 2500, -250), lty = "dotted", col = "white")
legend(
  x = 6400, ymax * 0.94,
  legend = c("SPD (dates normalised)", "Logistic Model", "95% MC envelope", "positive deviation", "negative deviation"),
  col = c(1, "black", "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)),
  lty = c(1, 2, 1, 1, 1), lwd = c(0.5, 0.5, 5, 5, 5), cex = 0.6, bg = "white", title = ""
)
text(x = 6350, ymax * 0.92, labels = "c. Logistic Fit", font = 2, cex = 0.8, adj = c(0, 0.7))
text(x = 6200, ymax * 0.52, cex = 0.7, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(logistic_model$pval, 4))))
box()

# Add text and segments for periods
segments(6450, ymax * 0.08, 6450, ymax * 0.01, lty = "solid", col = "purple", lwd = (1.2))
text(x = 5950, ymax * 0.01, labels = "Chalcolithic", lty = c("solid"), col = "purple", font = 2, cex = 0.7, adj = c(0.5, -0.5))
segments(6450, ymax * 0.08, 5500, ymax * 0.08, lty = "solid", col = "purple", lwd = (1.2))
segments(5500, ymax * 0.08, 5500, ymax * 0.01, lty = "solid", col = "red", lwd = (1.2))
text(x = 4900, ymax * 0.01, labels = "EBA", lty = c("solid"), col = "red", font = 2, cex = 0.7, adj = c(0.5, -0.5))
segments(4500, ymax * 0.08, 4500, ymax * 0.01, lty = "dotted", col = "red", lwd = (1.2))
text(x = 4000, ymax * 0.01, labels = "MBA", lty = c("solid"), col = "red", font = 2, cex = 0.7, adj = c(0.5, -0.5))
segments(3550, ymax * 0.08, 3550, ymax * 0.01, lty = "dotted", col = "red", lwd = (1.2))
text(x = 3350, ymax * 0.01, labels = "LBA", lty = c("solid"), col = "red", font = 2, cex = 0.7, adj = c(0.5, -0.5))
segments(3100, ymax * 0.08, 3100, ymax * 0.01, lty = "solid", col = "red", lwd = (1.2))
text(x = 4400, ymax * 0.13, labels = "Bronze Age", lty = c("solid"), col = "red", font = 2, cex = 0.7, adj = c(0.5, 0.8))
segments(5500, ymax * 0.08, 3100, ymax * 0.08, lty = "solid", col = "red", lwd = (1.2))
segments(3100, ymax * 0.08, 3100, ymax * 0.01, lty = "solid", col = "darkolivegreen", lwd = (1.2))
text(x = 2700, 0, labels = "Iron Age", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.7, adj = c(0.5, -3.3))
segments(3100, ymax * 0.08, 2489, ymax * 0.08, lty = "solid", col = "darkolivegreen", lwd = (1.2))
segments(2930, ymax * 0.03, 2489, ymax * 0.03, lty = "solid", col = "darkolivegreen", lwd = (1.2))
text(x = 3000, 0, labels = "I", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.65, adj = c(0.5, -1.4))
segments(2930, ymax * 0.08, 2930, ymax * 0, lty = "solid", col = "darkolivegreen", lwd = (1.2))
text(x = 2710, 0, labels = "II", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.65, adj = c(0.5, -1.4))
text(x = 2850, 0, labels = "A", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.55, adj = c(0.5, -0.2))
segments(2780, ymax * 0.03, 2780, ymax * 0, lty = "solid", col = "darkolivegreen", lwd = (1.2))
text(x = 2720, 0, labels = "B", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.55, adj = c(0.5, -0.2))
segments(2670, ymax * 0.03, 2670, ymax * 0, lty = "solid", col = "darkolivegreen", lwd = (1.2))
text(x = 2580, 0, labels = "C", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.55, adj = c(0.5, -0.2))
segments(2489, ymax * 0.08, 2489, ymax * 0, lty = "solid", col = "darkolivegreen", lwd = (1.2))
segments(2489, ymax * 0.08, 2283, ymax * 0.08, lty = "solid", col = "darkolivegreen", lwd = (1.2))
text(x = 2420, 0, labels = "III", lty = c("solid"), col = "darkolivegreen", font = 2, cex = 0.65, adj = c(0.5, -1.4))
segments(2750, ymax * 0.9, 2350, ymax * 0.9, lty = "solid", col = "grey23", lwd = (1))
text(x = 2550, ymax * 0.91, labels = "Hallstat Plataeu", lty = c("solid"), col = "grey23", font = 2, cex = 0.7, adj = c(0.5, -0.5))

# Add axes
xticks <- seq(6500, 2200, -250)
axis(side = 1, at = xticks, labels = xticks, las = 2, cex.axis = 0.8) # add BP axis
mtext("cal BP", 1, 2.3, at = 6362.5, adj = 0.05, font = 2, cex = 0.6, las = 2)
axis(side = 1, at = xticks - 50, labels = xticks - 2000, las = 2, cex.axis = 0.8, pos = -0.00017) # add BC/AD axis
mtext("BC", 1, 5.5, at = 6362.5, adj = -0.5, font = 2, cex = 0.6, las = 2)

dev.off()


## Comparison plots of logistic fit -----------------------------------------------------------------------

png(file = here::here("output/png/fig-04.png"), width = 1000, height = 1000, res = 150)
layout(matrix(c(1, 2), 2, 1, byrow = TRUE), widths = 8, heights = c(1.5, 2.1))
par(mar = c(0, 1, 0, 1)) # c(bottom, left, top, right)
ymax <- max(logistic_model_samaria$result$PrDens) * 1.1
plot(logistic_model_samaria, ylim = c(0, ymax), xlim = c(6500, 2200), drawaxes = FALSE)
lines(logistic_model_samaria$fit$calBP, logistic_model_samaria$fit$PrDens, col = "black", lty = "dashed", lwd = 0.5)
abline(v = seq(6500, 2500, -250), lty = "dotted", col = "white")
legend(
  x = 6400, ymax * 0.94,
  legend = c("SPD (dates normalised)", "Logistic Model", "95% MC envelope", "positive deviation", "negative deviation"),
  col = c(1, "black", "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)),
  lty = c(1, 2, 1, 1, 1), lwd = c(0.5, 0.5, 5, 5, 5), cex = 0.6, bg = "white", title = ""
)
text(x = 6350, ymax * 0.92, labels = "a. Logistic Fit - Samaria", font = 2, cex = 0.6, adj = c(0, 0.7))
text(x = 6200, ymax * 0.56, cex = 0.7, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(logistic_model_samaria$pval, 4))))
segments(2750, ymax * 0.9, 2350, ymax * 0.9, lty = "solid", col = "grey23", lwd = (1))
text(x = 2550, ymax * 0.91, labels = "Hallstat Plataeu", lty = c("solid"), col = "grey23", font = 2, cex = 0.45, adj = c(0.5, -0.5))
box()
par(mar = c(7, 1, 0, 1)) # c(bottom, left, top, right)
ymax <- max(logistic_model_judah$result$PrDens) * 1.1
plot(logistic_model_judah, ylim = c(0, ymax), xlim = c(6500, 2200), drawaxes = FALSE)
lines(logistic_model_judah$fit$calBP, logistic_model_judah$fit$PrDens, col = "black", lty = "dashed", lwd = 0.5)
abline(v = seq(6500, 2500, -250), lty = "dotted", col = "white")
legend(
  x = 6400, ymax * 0.94,
  legend = c("SPD (dates normalised)", "Logistic Model", "95% MC envelope", "positive deviation", "negative deviation"),
  col = c(1, "black", "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)),
  lty = c(1, 2, 1, 1, 1), lwd = c(0.5, 0.5, 5, 5, 5), cex = 0.6, bg = "white", title = ""
)
text(x = 6350, ymax * 0.92, labels = "b. Logistic Fit - Judah", font = 2, cex = 0.6, adj = c(0, 0.7))
text(x = 6200, ymax * 0.52, cex = 0.7, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(logistic_model_judah$pval, 4))))
segments(2750, ymax * 0.9, 2350, ymax * 0.9, lty = "solid", col = "grey23", lwd = (1))
text(x = 2550, ymax * 0.91, labels = "Hallstat Plataeu", lty = c("solid"), col = "grey23", font = 2, cex = 0.45, adj = c(0.5, -0.5))
box()
xticks <- seq(6500, 2200, -250)
axis(side = 1, at = xticks, labels = xticks, las = 2, cex.axis = 0.8) # add BP axis
mtext("cal BP", 1, 2.3, at = 6362.5, adj = 0, font = 2, cex = 0.65, las = 2)
axis(side = 1, at = xticks - 50, labels = xticks - 2000, las = 2, cex.axis = 0.8, pos = -0.00017) # add BC/AD axis
mtext("BC", 1, 5.3, at = 6362.5, adj = -0.5, font = 2, cex = 0.65, las = 2)
dev.off()

## Plot Regional Permutation Test (6500-2200 BP) -----------------------------------------------------------------------

png(file = here::here("output/png/fig-05.png"), width = 1200, height = 1000, res = 150)
layout(matrix(c(1, 2), 2, 1, byrow = TRUE), widths = 8, heights = c(1.8, 2.2))
par(mar = c(0, 1, 1, 1))
par(yaxs = "i")
par(xaxs = "i")
xticks <- seq(5000, 1000, -250)
ymax <- max(perm_test$envelope[["Samaria"]][, 2], perm_test$observed[["Samaria"]]$PrDens) * 1.1
plot(perm_test, focalm = "Samaria", xlim = c(6500, 2200), col.obs = "brown", lwd.obs = 1, drawaxes = FALSE)
abline(v = seq(6250, 3000, -500), lty = "dotted", col = "white")
legend(x = 6400, y = ymax * 0.90, legend = c("SPD", "95% MC envelope", "positive deviation", "negative deviation"), col = c(1, "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)), lty = c(1, 1, 1, 1), lwd = c(0.5, 5, 5, 5), cex = 0.6, bg = "white", title = expression(bold("a. Samaria")))
text(x = 6400, y = ymax * 0.55, labels = paste("n=", nrow(rc_dates[rc_dates$Region == "Samaria", ]), ", sites=", length(unique(rc_dates$SiteName[rc_dates$Region == "Samaria"])), ", bins=", length(unique(bins[rc_dates$Region == "Samaria"])), sep = ""), font = 1, cex = 0.6, adj = c(0, 0.7))
text(x = 6400, y = ymax * 0.48, cex = 0.6, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(min(perm_test$pValueList["Samaria"], 1 - perm_test$pValueList["Samaria"]), 3))))
segments(8000, ymax * 0.15, 8000, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 7500, ymax * 0.01, labels = "Neolithic", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(6450, ymax * 0.15, 6450, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 5950, ymax * 0.01, labels = "Chalcolithic", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(5500, ymax * 0.15, 5500, ymax * 0.01, lty = "solid", col = "brown", lwd = (1.2))
text(x = 4900, ymax * 0.01, labels = "EBA", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(4500, ymax * 0.12, 4500, ymax * 0.01, lty = "dotted", col = "brown", lwd = (1.2))
text(x = 4000, ymax * 0.01, labels = "MBA", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(3550, ymax * 0.12, 3550, ymax * 0.01, lty = "dotted", col = "brown", lwd = (1.2))
text(x = 3350, ymax * 0.01, labels = "LBA", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(3100, ymax * 0.15, 3100, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 4400, ymax * 0.13, labels = "Bronze Age", lty = c("solid"), col = c("brown"), font = 2, cex = 0.9, adj = c(0.5, -2))
segments(3100, ymax * 0.15, 3100, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 2700, 0, labels = "Iron Age", lty = c("solid"), col = c("brown"), font = 2, cex = 0.8, adj = c(0.5, -5.2))
segments(3100, ymax * 0.15, 2308, ymax * 0.15, lty = "solid", col = "brown", lwd = (1.2))
segments(2930, ymax * 0.09, 2489, ymax * 0.09, lty = "solid", col = "brown", lwd = (1.2))
text(x = 3005, ymax * 0.01, labels = "I", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -2.9))
segments(2930, ymax * 0.15, 2930, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 2723, ymax * 0.01, labels = "II", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -2.9))
text(x = 2850, ymax * 0.01, labels = "A", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.8))
segments(2780, ymax * 0.09, 2780, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 2720, ymax * 0.01, labels = "B", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.8))
segments(2670, ymax * 0.09, 2670, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 2583, ymax * 0.010, labels = "C", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -0.8))
segments(2489, ymax * 0.15, 2489, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
segments(2489, ymax * 0.15, 2280, ymax * 0.15, lty = "solid", col = "brown", lwd = (1.2))
segments(2280, ymax * 0.15, 2280, ymax * 0.02, lty = "solid", col = "brown", lwd = (1.2))
text(x = 2385, 0, labels = "III", lty = c("solid"), col = "brown", font = 2, cex = 0.8, adj = c(0.5, -2.9))
box()
par(mar = c(6, 1, 0, 1))
par(yaxs = "i")
par(xaxs = "i")
xticks <- seq(5000, 1000, -250)
ymax <- max(perm_test$envelope[["Judah"]][, 2], perm_test$observed[["Judah"]]$PrDens) * 1.1
plot(perm_test, focalm = "Judah", xlim = c(6500, 2200), col.obs = "purple", lwd.obs = 1, drawaxes = FALSE)
abline(v = seq(4750, 1600, -500), lty = "dotted", col = "white")
legend(x = 6400, ymax * 0.90, legend = c("SPD", "95% MC envelope", "positive deviation", "negative deviation"), col = c(1, "lightgrey", rgb(0.7, 0, 0, 0.2), rgb(0, 0, 0.7, 0.2)), lty = c(1, 1, 1, 1), lwd = c(0.5, 5, 5, 5), cex = 0.6, bg = "white", title = expression(bold("b. Judah")))
text(x = 6400, y = ymax * 0.52, labels = paste("n=", nrow(rc_dates[rc_dates$Region == "Judah", ]), ", sites=", length(unique(rc_dates$SiteName[rc_dates$Region == "Judah"])), ", bins=", length(unique(bins[rc_dates$Region == "Judah"])), sep = ""), font = 1, cex = 0.6, adj = c(0, 0.7))
text(x = 6400, y = ymax * 0.45, cex = 0.6, adj = c(0, 0.7), labels = substitute(paste(italic(p), "=", x, sep = ""), list(x = round(min(perm_test$pValueList["Judah"], 1 - perm_test$pValueList["Judah"]), 3))))
segments(8000, ymax * 0.15, 8000, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 7500, ymax * 0.01, labels = "Neolithic", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(6450, ymax * 0.15, 6450, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 5950, ymax * 0.01, labels = "Chalcolithic", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(5500, ymax * 0.15, 5500, ymax * 0.01, lty = "solid", col = "purple", lwd = (1.2))
text(x = 4900, ymax * 0.01, labels = "EBA", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(4500, ymax * 0.12, 4500, ymax * 0.01, lty = "dotted", col = "purple", lwd = (1.2))
text(x = 4000, ymax * 0.01, labels = "MBA", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(3550, ymax * 0.12, 3550, ymax * 0.01, lty = "dotted", col = "purple", lwd = (1.2))
text(x = 3350, ymax * 0.01, labels = "LBA", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.5))
segments(3100, ymax * 0.15, 3100, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 4400, ymax * 0.13, labels = "Bronze Age", lty = c("solid"), col = c("purple"), font = 2, cex = 0.9, adj = c(0.5, -0.9))
segments(3100, ymax * 0.15, 3100, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 2700, 0, labels = "Iron Age", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -5))
segments(3100, ymax * 0.15, 2308, ymax * 0.15, lty = "solid", col = "purple", lwd = (1.2))
segments(2930, ymax * 0.09, 2489, ymax * 0.09, lty = "solid", col = "purple", lwd = (1.2))
text(x = 3005, 0, labels = "I", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -2.9))
segments(2930, ymax * 0.15, 2930, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 2723, 0, labels = "II", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -2.9))
text(x = 2850, 0, labels = "A", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.8))
segments(2780, ymax * 0.09, 2780, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 2720, 0, labels = "B", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.8))
segments(2670, ymax * 0.09, 2670, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 2583, 0, labels = "C", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -0.8))
segments(2489, ymax * 0.15, 2489, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
segments(2489, ymax * 0.15, 2280, ymax * 0.15, lty = "solid", col = "purple", lwd = (1.2))
segments(2280, ymax * 0.15, 2280, ymax * 0.02, lty = "solid", col = "purple", lwd = (1.2))
text(x = 2385, 0, labels = "III", lty = c("solid"), col = "purple", font = 2, cex = 0.8, adj = c(0.5, -2.9))
box()
par(mar = c(0, 1, 0, 1))
par(yaxs = "i")
par(xaxs = "i")
xticks <- seq(6500, 2150, -250)
axis(side = 1, at = xticks, labels = xticks, las = 2, cex.axis = 0.6, pos = 0) # add BP axis
mtext("cal BP", 1, 2.5, at = 6362.5, adj = -4.5, font = 2, cex = 0.6, las = 2)
axis(side = 1, at = xticks - 50, labels = xticks - 2000, las = 2, cex.axis = 0.6, pos = -0.00017) # add BC/AD axis
mtext("BC", 1, 5.5, at = 6362.5, adj = -11, font = 2, cex = 0.6, las = 2)
par(yaxs = "r")
par(xaxs = "r")
dev.off()
