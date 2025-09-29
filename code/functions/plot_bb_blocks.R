# Plot Boom and Bust Blocks from a CalGrid Object
#
# This function extracts boom and bust periods from a rcarbon SpdModelTest object
# and visualizes them as colored blocks on the current plot.
# This function is a hack to avoid overcomplicating the generation of stacked
# plots of paleoclimate proxies and SPDs as seen in fig-07
#
# Inputs:
#   x: A SpdModelTest object from rcarbon::modelTest() containing observed and
#       expected (envelope) SPD values
#       The function has been tested with a logistic model
#       but should work with an exponential one as well
#
# Details:
#   The function identifies "boom" periods where the observed SPD is above the upper
#   confidence envelope, and "bust" periods where it is below the lower confidence
#   envelope. It then draws these periods as semi-transparent colored blocks on the
#   current plot: red for booms and blue for busts. The blocks span the entire
#   vertical range of the current plot.
#
#   This function must be called after a plot has been created, as it adds elements
#   to the existing plot rather than creating a new one.
#   The function might not be very flexible and has been tested only with the
#   layout set in the 03_paleoclimate script.
#
#   The code is adapted from the rcarbon package's plotting functions (https://github.com/ahb108/rcarbon/blob/master/R/plots.R).
#
# Returns:
#   No return value; the function adds visual elements to the current plot
plot_bb_blocks <- function(x) {
  obs <- x$result[, 1:2]
  envelope <- x$result[, 3:4]
  obs$Years <- obs$calBP

  booms <- which(obs$PrDens > envelope[, 2])
  busts <- which(obs$PrDens < envelope[, 1])
  baseline <- rep(NA, nrow(obs))

  boomPlot <- baseline
  if (length(booms) > 0) {
    boomPlot[booms] <- obs[booms, 2]
  }
  bustPlot <- baseline
  if (length(busts) > 0) {
    bustPlot[busts] <- obs[busts, 2]
  }
  boomBlocks <- vector("list")
  counter <- 0
  state <- "off"
  for (i in 1:length(boomPlot)) {
    if (!is.na(boomPlot[i]) & state == "off") {
      counter <- counter + 1
      boomBlocks <- c(boomBlocks, vector("list", 1))
      boomBlocks[[counter]] <- vector("list", 2)
      boomBlocks[[counter]][[1]] <- boomPlot[i]
      boomBlocks[[counter]][[2]] <- obs[i, "Years"]
      state <- "on"
    }
    if (state == "on") {
      if (!is.na(boomPlot[i])) {
        boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]], boomPlot[i])
        boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]], obs[i, "Years"])
      }
      if (is.na(boomPlot[i])) {
        state <- "off"
      }
    }
  }
  bustBlocks <- vector("list")
  counter <- 0
  state <- "off"
  for (i in 1:length(bustPlot)) {
    if (!is.na(bustPlot[i]) & state == "off") {
      counter <- counter + 1
      bustBlocks <- c(bustBlocks, vector("list", 1))
      bustBlocks[[counter]] <- vector("list", 2)
      bustBlocks[[counter]][[1]] <- bustPlot[i]
      bustBlocks[[counter]][[2]] <- obs[i, "Years"]
      state <- "on"
    }
    if (state == "on") {
      if (!is.na(bustPlot[i])) {
        bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]], bustPlot[i])
        bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]], obs[i, "Years"])
      }
      if (is.na(bustPlot[i])) {
        state <- "off"
      }
    }
  }

  # Get current plot boundaries
  plot_info <- par("usr")
  y_min <- plot_info[3]
  y_max <- plot_info[4]

  if (length(booms) > 0) {
    for (i in 1:length(boomBlocks)) {
      polygon(c(boomBlocks[[i]][[2]], rev(boomBlocks[[i]][[2]])),
              c(rep(y_max, length(boomBlocks[[i]][[1]])), rep(y_min, length(boomBlocks[[i]][[1]]))),
              col = rgb(0.7, 0, 0, 0.2), border = NA)
    }
  }

  if (length(busts) > 0) {
    for (i in 1:length(bustBlocks)) {
      polygon(c(bustBlocks[[i]][[2]], rev(bustBlocks[[i]][[2]])),
              c(rep(y_max, length(bustBlocks[[i]][[1]])), rep(y_min, length(bustBlocks[[i]][[1]]))),
              col = rgb(0, 0, 0.7, 0.2), border = NA)
    }
  }
}
