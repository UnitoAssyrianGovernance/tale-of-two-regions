## Calculate running correlation between paleoclimate data and archaeodemographic proxies
## Last edited 2025-09-22
## Authors: Andrea Titolo (refactored code), Alessio Palmisano (original code for the running correlation function)
## Packages version at the time of analysis
## here: 1.0.1
## gtools: 3.9.5

# Set-up -----------------------------------------------------------------------

# Load required libraries
required_packages <- c("gtools", "here")

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

# Source functions
source(here::here("code/functions/generate_time_slices.R"))
source(here::here("code/functions/compute_running_correlation.R"))

# Load paleoclimate data
paleoclimate <- read.csv(here::here("data/processed/csv/paleoclimate_all.csv"), header = TRUE)

# load all the SPDs data generated from the radiocarbon script
load(here::here("output/rda/radiocarbon_data.RData"))

# Load archaeological data (all region)
load(here::here(paste0("output/rda/samaria_judah_norm_", 100, ".RData")))
load(here::here(paste0("output/rda/samaria_norm_", 100, ".RData")))
load(here::here(paste0("output/rda/judah_norm_", 100, ".RData")))

# If necessary, we can load the data generated from this script to avoid recomputing everything
if (file.exists(here::here("data/rda/corellation_data.RData"))) {
  load(here::here("output/rda/corellation_data.RData"))
} else {
  print("RData file not found, please run the rest of the script and save the data")
}

# Variables for running correlation
jeita_slice_size <- 50
jeita_blockwidth <- 50
jeita_cor_width <- 10

# Paleoclimate vs SPD -----------------------------------------------------------------------

## Variables -----------------------------------------------------------------------

analysis_window <- c(start = 2201, end = 6500)

# Limit the paleoclimate data from the start to the end of the analysis window
paleoclimate <- paleoclimate[which(paleoclimate$bp >= analysis_window[["start"]] & paleoclimate$bp <= analysis_window[["end"]]), ]

# Create a list with the SPDs data needed for the analysis
spd_l <- list(
  allspd = all_spd_norm,
  samaria = all_spd_norm_samaria,
  judah = all_spd_norm_judah
)

# Create a list of archaeo proxies
aoristic_l <- list(
  all = archaeo_proxies_sum_l$site_blocksum_weights,
  samaria = sam_archaeo_proxies_sum_l$site_blocksum_weights,
  judah = jud_archaeo_proxies_sum_l$site_blocksum_weights
)

sim_df_l <- list(
  all = conf_intervals_l$medcount,
  samaria = sam_conf_intervals_l$medcount,
  judah = jud_conf_intervals_l$medcount
)

count_l <- list(
  all = archaeo_proxies_sum_l$site_blocksum_counts,
  samaria = sam_archaeo_proxies_sum_l$site_blocksum_counts,
  judah = jud_archaeo_proxies_sum_l$site_blocksum_counts
)

area_l <- list(
  all = archaeo_proxies_sum_l$site_blocksum_area,
  samaria = sam_archaeo_proxies_sum_l$site_blocksum_area,
  judah = jud_archaeo_proxies_sum_l$site_blocksum_area
)

## Jeita Cave -----------------------------------------------------------------------

jeita <- paleoclimate[, c("bp", "jeita_cave")]

jeita_list <- list()

# Loop through the spd and apply the gen_time_slices() function
for (i in seq_along(spd_l)) {
  # Apply the gen_time_slices() function to each SPD
  jeita_list[[i]] <- gen_time_slices(
    paleoclim = jeita,
    cname_paleoclim = "jeita_cave",
    cname_bp = "bp",
    spd = spd_l[[i]],
    aorist = aoristic_l[[i]],
    simulated = sim_df_l[[i]],
    area = area_l[[i]],
    counts = count_l[[i]],
    start = analysis_window["start"],
    end = analysis_window["end"],
    range = jeita_slice_size
  )
  names(jeita_list)[i] <- paste0("jeita_", names(spd_l[i])) # TODO use aoristic_l names
}

### Correlation -----------------------------------------------------------------------

# Apply the correlation test on the list elements and extract the correlation coefficient and p-values
result_jeita <- sapply(names(jeita_list), function(x) {
  cor.test(jeita_list[[x]]$spd, jeita_list[[x]]$paleoclim, method = "pearson")[c("estimate", "p.value")]
}, simplify = FALSE)

names(result_jeita) <- names(jeita_list)

result_jeita <- do.call(rbind, lapply(result_jeita, function(x) {
  data.frame(
    estimate = x[1],
    pvalue = x[2]
  )
}))

### Running-correlation  -----------------------------------------------------------------------

jeita_corr_list <- running_cor(
  paleoclim = jeita_list$jeita_allspd$paleoclim,
  spd = jeita_list$jeita_allspd$spd,
  aorist = jeita_list$jeita_allspd$aorist_sums,
  simulated = jeita_list$jeita_allspd$simulated_sums,
  area = jeita_list$jeita_allspd$area_sums,
  counts = jeita_list$jeita_allspd$counts_sums,
  blockwidth = jeita_blockwidth,
  start = analysis_window["start"] -1,
  end = analysis_window["end"],
  method = "pearson",
  cor_width = jeita_cor_width
)

jeita_corr_list_sam <- running_cor(
  paleoclim = jeita_list$jeita_samaria$paleoclim,
  spd = jeita_list$jeita_samaria$spd,
  aorist = jeita_list$jeita_samaria$aorist_sums,
  simulated = jeita_list$jeita_samaria$simulated_sums,
  area = jeita_list$jeita_samaria$area_sums,
  counts = jeita_list$jeita_samaria$counts_sums,
  blockwidth = jeita_blockwidth,
  start = analysis_window["start"] -1,
  end = analysis_window["end"],
  method = "pearson",
  cor_width = jeita_cor_width
)

jeita_corr_list_jud <- running_cor(
  paleoclim = jeita_list$jeita_judah$paleoclim,
  spd = jeita_list$jeita_judah$spd,
  aorist = jeita_list$jeita_judah$aorist_sums,
  simulated = jeita_list$jeita_judah$simulated_sums,
  area = jeita_list$jeita_judah$area_sums,
  counts = jeita_list$jeita_judah$counts_sums,
  blockwidth = jeita_blockwidth,
  start = analysis_window["start"] -1,
  end = analysis_window["end"],
  method = "pearson",
  cor_width = jeita_cor_width
)


# Save data -----------------------------------------------------------------------

save(result_jeita, jeita_corr_list,
  file = here::here("output/rda/corellation_data.RData")
)

# Plotting -----------------------------------------------------------------------

x_ticks <- seq(analysis_window[[2]], analysis_window[[1]], -250)

png(file = here::here(paste0("output/png/fig-08.png")), width = 2000, height = 1400, res = 150)
layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE), widths = 8, heights = c(2.5, 2.5, 2.5))
# layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = TRUE), heights = c(2.5, 2.5, 2.5))
par(mar = c(1, 1, 0, 1))
par(yaxs = "i")
par(xaxs = "i")

par(mar = c(2, 5, 0.5, 1)) # c(bottom, left, top, right)
mean_window <- (rowMeans(jeita_corr_list$spd_all_corr[, 1:2]))
plot(jeita_corr_list$spd_all_corr$run_cor ~ mean_window, lty = "solid", col = "orange", cex.axis = 1, type = "l",
xlab = "", ylab = "Corr. Coefficent", xlim = c(analysis_window[[2]], analysis_window[[1]]), ylim = c(-1, 1), xaxt = "n")
axis(side = 1, at = x_ticks, labels = x_ticks, cex.axis = 0.9)
lines(mean_window, jeita_corr_list$aorist_all_corr$run_cor, col = "black", lty = "dashed")
lines(mean_window, jeita_corr_list$sim_all_corr$run_cor, col = "grey", lty = "solid")
lines(mean_window, jeita_corr_list$area_all_corr$run_cor, col = "red", lty = "solid")
lines(mean_window, jeita_corr_list$counts_all_corr$run_cor, col = "black", lty = "solid")
polygon(
  x = c(analysis_window[[2]], analysis_window[[1]], analysis_window[[1]], analysis_window[[2]]),
  y = c(-1, -1, -0.4, -0.4),
  col = rgb(207, 218, 240, alpha = 60, maxColorValue = 255),
  border = NA
)
polygon(
  x = c(analysis_window[[2]], analysis_window[[1]], analysis_window[[1]], analysis_window[[2]]),
  y = c(1, 1, 0.4, 0.4),
  col = rgb(250, 221, 202, alpha = 60, maxColorValue = 255),
  border = NA
)
abline(h = 0, col = "grey91", lty = "dashed")
legend(x = 6500, y = 0.93,
       legend = c("vs SPD", "vs aoristic weight", "vs randomised start date", "vs total area", "vs raw count"),
       lty = c("solid", "dashed", "solid", "solid", "solid"),
       lwd = c(1, 1, 1, 1, 1),
       col = c("orange", "black", "grey", "red", "black"),
       bty = "n", cex = 1, bg = "transparent", title = "")
# legend(x = 6400, y = 0.93, legend = c("SPD", "aoristic weight", "randomised start date"), lty = c("solid", "solid", "solid"), lwd = c(1, 1, 1), col = c("black", "blue", "red"), bty = "n", cex = 1, bg = "transparent", title = "")
text(x = 6450, y = 0.93, labels = "a. Jeita Cave (All)", font = 2, cex = 1.4, adj = c(0, 0.7))
text(x = 3100, y = 0.85,
  labels = paste0("Bin size: ", jeita_slice_size,
  " years\nMoving window: ", jeita_corr_list$moving_window, " years",
  "\nSteps per window: ", jeita_cor_width, " steps",
  "\nCorrelation step size: ", jeita_blockwidth, " years"),
     font = 1, cex = 1.1, adj = c(0, 0.7)
)

# par(mar = c(2, 5, 0.5, 1)) # c(bottom, left, top, right)
mean_window <- (rowMeans(jeita_corr_list_sam$spd_all_corr[, 1:2]))
plot(jeita_corr_list_sam$spd_all_corr$run_cor ~ mean_window, lty = "solid", col = "orange", cex.axis = 1, type = "l",
xlab = "", ylab = "Corr. Coefficent", xlim = c(analysis_window[[2]], analysis_window[[1]]), ylim = c(-1, 1), xaxt = "n")
axis(side = 1, at = x_ticks, labels = x_ticks, cex.axis = 0.9)
lines(mean_window, jeita_corr_list_sam$aorist_all_corr$run_cor, col = "black", lty = "dashed")
lines(mean_window, jeita_corr_list_sam$sim_all_corr$run_cor, col = "grey", lty = "solid")
lines(mean_window, jeita_corr_list_sam$area_all_corr$run_cor, col = "red", lty = "solid")
lines(mean_window, jeita_corr_list_sam$counts_all_corr$run_cor, col = "black", lty = "solid")
polygon(
  x = c(analysis_window[[2]], analysis_window[[1]], analysis_window[[1]], analysis_window[[2]]),
  y = c(-1, -1, -0.4, -0.4),
  col = rgb(207, 218, 240, alpha = 60, maxColorValue = 255),
  border = NA
)
polygon(
  x = c(analysis_window[[2]], analysis_window[[1]], analysis_window[[1]], analysis_window[[2]]),
  y = c(1, 1, 0.4, 0.4),
  col = rgb(250, 221, 202, alpha = 60, maxColorValue = 255),
  border = NA
)
abline(h = 0, col = "grey91", lty = "dashed")
legend(x = 6500, y = 0.93,
       legend = c("vs SPD", "vs aoristic weight", "vs randomised start date", "vs total area", "vs raw count"),
       lty = c("solid", "dashed", "solid", "solid", "solid"),
       lwd = c(1, 1, 1, 1, 1),
       col = c("orange", "black", "grey", "red", "black"),
       bty = "n", cex = 1, bg = "transparent", title = "")
# text(x = 6400, y = 0.81, labels = "Proxies", font = 2, cex = 1.2, adj = c(0, 0))
text(x = 6450, y = 0.93, labels = "b. Jeita Cave (Samaria)", font = 2, cex = 1.4, adj = c(0, 0.7))
text(x = 3100, y = 0.85,
  labels = paste0("Bin size: ", jeita_slice_size,
  " years\nMoving window: ", jeita_corr_list_sam$moving_window, " years",
  "\nSteps per window: ", jeita_cor_width, " steps",
  "\nCorrelation step size: ", jeita_blockwidth, " years"),
     font = 1, cex = 1.1, adj = c(0, 0.7)
)

# par(mar = c(2, 5, 0.5, 1)) # c(bottom, left, top, right)
mean_window <- (rowMeans(jeita_corr_list_jud$spd_all_corr[, 1:2]))
plot(jeita_corr_list_jud$spd_all_corr$run_cor ~ mean_window, lty = "solid", col = "orange", cex.axis = 1, type = "l",
xlab = "", ylab = "Corr. Coefficent", xlim = c(analysis_window[[2]], analysis_window[[1]]), ylim = c(-1, 1), xaxt = "n")
axis(side = 1, at = x_ticks, labels = x_ticks, cex.axis = 0.9)

lines(mean_window, jeita_corr_list_jud$aorist_all_corr$run_cor, col = "black", lty = "dashed")
lines(mean_window, jeita_corr_list_jud$sim_all_corr$run_cor, col = "grey", lty = "solid")
lines(mean_window, jeita_corr_list_jud$area_all_corr$run_cor, col = "red", lty = "solid")
lines(mean_window, jeita_corr_list_jud$counts_all_corr$run_cor, col = "black", lty = "solid")
polygon(
  x = c(analysis_window[[2]], analysis_window[[1]], analysis_window[[1]], analysis_window[[2]]),
  y = c(-1, -1, -0.4, -0.4),
  col = rgb(207, 218, 240, alpha = 60, maxColorValue = 255),
  border = NA
)
polygon(
  x = c(analysis_window[[2]], analysis_window[[1]], analysis_window[[1]], analysis_window[[2]]),
  y = c(1, 1, 0.4, 0.4),
  col = rgb(250, 221, 202, alpha = 60, maxColorValue = 255),
  border = NA
)
abline(h = 0, col = "grey91", lty = "dashed")
legend(x = 6500, y = 0.93,
       legend = c("vs SPD", "vs aoristic weight", "vs randomised start date", "vs total area", "vs raw count"),
       lty = c("solid", "dashed", "solid", "solid", "solid"),
       lwd = c(1, 1, 1, 1, 1),
       col = c("orange", "black", "grey", "red", "black"),
       bty = "n", cex = 1, bg = "transparent", title = "")
# text(x = 6500, y = 0.81, labels = "Proxies", font = 2, cex = 1.2, adj = c(0, 0))
text(x = 6450, y = 0.93, labels = "c. Jeita Cave (Judah)", font = 2, cex = 1.4, adj = c(0, 0.7))
text(x = 3100, y = 0.85,
  labels = paste0("Bin size: ", jeita_slice_size,
  " years\nMoving window: ", jeita_corr_list_jud$moving_window, " years",
  "\nSteps per window: ", jeita_cor_width, " steps",
  "\nCorrelation step size: ", jeita_blockwidth, " years"),
     font = 1, cex = 1.1, adj = c(0, 0.7)
)
dev.off()
