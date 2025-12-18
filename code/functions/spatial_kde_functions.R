# Functions for Spatial Kernel Density Estimation
# Note: sf is a required dependency because there is no easy way to generate a polygonal
# owin from a terra object.
# If another method arise, the sf dependency will be dropped.

# Generate Kernel Density Estimation
#
# Creates a kernel density estimation from spatial points data
#
# Inputs:
#   x: An sf object containing point data
#   bandwidth: The bandwidth (sigma) for the kernel density estimation
#   cell_res: The cell resolution for the output density surface
#   spat_window: An sf defining the study area window
#   weight_value: (Optional) Name of the column to use for weighted KDE
#
# Returns:
#   A spatstat density object (pixel image)
get_kde <- function(x, bandwidth, cell_res, spat_window, weight_value = NULL) {

  win <- as.owin(spat_window)

  if (is.null(weight_value)) {

    xy <- as.data.frame(sf::st_coordinates(x))
    ppp_obj <- spatstat.geom::ppp(x = xy$X, y = xy$Y, window = win)

    dens <- spatstat.explore::density.ppp(ppp_obj, sigma = bandwidth, eps = cell_res,
      edge = TRUE, diggle = TRUE)

  } else if (!is.null(weight_value)) {

    x <- x[!is.na(x[[weight_value]]), ]
    marks_data <- x[[weight_value]]

    xy <- as.data.frame(sf::st_coordinates(x))
    ppp_obj <- spatstat.geom::ppp(x = xy$X, y = xy$Y, window = win, marks = marks_data)

    dens <- spatstat.explore::density.ppp(ppp_obj, sigma = bandwidth, eps = cell_res,
      edge = TRUE, weights = marks_data, diggle = TRUE)
  }

  return(dens)
}

# Convert Density to Terra SpatRaster
#
# Converts a spatstat density object to a terra SpatRaster for easier plotting
#
# Inputs:
#   dens: A spatstat density object from get_kde()
#   crs_source: A source of coordinate reference system information
#   spat_window: The spatial window used for the density estimation
#
# Returns:
#   A terra SpatRaster object with appropriate CRS and extent
get_terra_dens <- function(dens, crs_source, spat_window) {

  # Extract bandwidth attribute of density map
  sigma <- attributes(dens)$sigma

  dens <- terra::rast(dens)
  terra::crs(dens) <- crs(crs_source)
  terra::ext(dens) <- ext(vect(spat_window))

  # Set metadata to the spatraster object
  metags(dens) <- paste0("sigma=", sigma)

  return(dens)
}

# Generate Multiple KDE Maps
#
# Wrapper function that creates KDE maps for multiple time periods using get_kde()
# the arguments of this functions will be passed to get_kde()
#
# Inputs:
#   sites_list: List of sf objects with site data for each period
#   spat_window: Spatial window for the analysis - must be an sf (see top)
#   bandwidth: Bandwidth parameter for KDE
#   cell_res: Cell resolution for the output
#   periods: Vector of period names
#   weight: Logical; whether to use weighted KDE
#   weight_value: Name of column to use for weights (if weight = TRUE)
#
# Returns:
#   A list of KDE objects, one for each period
generate_kde_maps <- function(sites_list, spat_window, bandwidth, cell_res, periods,
  weight = FALSE, weight_value = NULL) {
  if (bandwidth <= 0 || cell_res <= 0) {
    stop("bandwidth and cell_res must be larger than 0")
  }

  dens_list <- list()
  pb <- txtProgressBar(min = 0, max = length(periods), style = 3)
  for (i in seq_along(periods)) {
    dens_list[[i]] <- get_kde(
      sites_list[[i]],
      bandwidth = bandwidth,
      cell_res = cell_res,
      spat_window = spat_window,
      if (weight) weight_value = weight_value else weight_value = NULL
    )
    print(paste("Calculating KDE for", periods[i], "with bandwidth", bandwidth, "and", nrow(sites_list[[i]]), "sites.."))
    names(dens_list[i]) <- periods[i]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n=====DONE=====\n")
  return(dens_list)
}

# Convert KDE Objects to Terra Format
#
# Converts a list of KDE objects to terra SpatRaster format for plotting with other terra objects
# Wrapper function around get_terra_dens()
#
# Inputs:
#   dens_list: List of spatstat density objects
#   crs_source: Source for coordinate reference system
#   spat_window: Spatial window used for the analysis
#   periods: Vector of period names
#
# Returns:
#   A list of terra SpatRaster objects, one for each period
convert_kde_to_terra <- function(dens_list, crs_source, spat_window, periods) {
  spat_window <- st_transform(spat_window, crs(dem))

  terra_dens_list <- list()
  pb <- txtProgressBar(min = 0, max = length(periods), style = 3)
  for (i in seq_along(periods)) {
    terra_dens_list[[i]] <- get_terra_dens(
      dens = dens_list[[i]],
      crs_source = crs_source,
      spat_window = spat_window
    )
    names(terra_dens_list[i]) <- periods[i]
    print(paste("Converting KDE for", periods[i], "to terra object.."))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(terra_dens_list)
}

# Plot Density Map with DEM Background
#
# Creates a visualization of a density map with a DEM backdrop + other data
#
# Inputs:
#   kde: A terra SpatRaster with density values
#   dem: A terra SpatRaster with elevation data
#   rivers: A vector object with river features
#   lakes: A vector object with lake features
#   dem_pal: Color palette for the DEM
#   dens_pal: Color palette for the density surface
#   zoom_ext: Extent to zoom to (this is passed from the wrapper function below)
#     but is generally equal to spat_window (the study area)
#   buffer_level: Buffer distance around the zoom extent to include a larger area
#   plot_title: Title for the plot
#
# Returns:
#   No return value; creates a plot
# Note: it might be better to make the additional objects optional in the future
plot_dens_with_dem <- function(kde, dem, rivers, lakes, plot_study_area = TRUE, study_area, dem_pal, dens_pal,
  zoom_ext, buffer_level, plot_title) {

  dem_pal <- map.pal(dem_pal, n = 1000)
  dens_pal <- map.pal(dens_pal, n = 1000)

  zoom_extent <- sf::st_transform(zoom_ext, sf::st_crs(dem))
  zoom_extent <- terra::ext(terra::buffer(vect(zoom_extent), buffer_level))

  # crop the dem to the zoom extent
  dem <- terra::crop(dem, zoom_extent)

  # Extract the actual extent of the plots for easier placing of north and scale
  extent <- terra::ext(dem)

  # Extract the bandwidth metadata from the spatraster object
  sigma <- as.numeric(metags(kde, name = "sigma")$value)

  # Crop the rivers and lakes to the zoom extent
  lakes <- terra::crop(vect(lakes), zoom_extent)
  rivers <- terra::crop(vect(rivers), zoom_extent)
  sites_to_plot <- terra::project(vect(sites[sites$Name %in% c("Jerusalem", "Samaria"), ]), crs(rivers))

  plot(dem, col = dem_pal, box = FALSE, axes = FALSE, maxcell = 500000,
    buffer = TRUE, legend = FALSE, main = plot_title, mar = c(1.1, 0.5, 2.1, 3.1)
  )
  lines(rivers, col = "#18b2dc", alpha = 0.9)
  polys(lakes, col = "#1bc0f2", border = NA)
  plot(kde, col = scales::alpha(dens_pal, 0.75), add = TRUE, buffer = TRUE,
    tick = "out", plg = list(nudge = 0.08, format = "g"), mar = c(1.1, 0.5, 2.1, 7.1)
  )
  if (plot_study_area) {
    study_area_vect <- terra::vect(study_area)
    study_area_vect <- terra::project(study_area_vect, terra::crs(dem))
    polys(study_area_vect, border = "yellow", lwd = 0.8)
  }
  text(sites_to_plot, labels = sites_to_plot$Name, pos = 4, offset = 0.2, col = "black", halo = TRUE, hc = "white", hw = 0.25)
  north(type = 1, xy = c(extent$xmin + 0.18, extent$ymax - 0.2), head = 0.05, d = 0.1)
  sbar(d = 25,xy = c(extent$xmin + 0.05, extent$ymax - 0.37),type = "line",
    divs = 2, below = "km", label = c(0, "", 25), cex = 0.9, ticks = FALSE,
    halo = FALSE, adj = c(0.5, -0.5)
  )
  add_mtext(paste0("σ = ", sigma/1000, " km"), side = 2, cex = 1, col = 1, line = -4, srt = 0)
}

# Plot Multiple Density Maps
#
# Creates a multi-panel figure with density maps for multiple periods
# Wraps around plot_dens_with_dem()
#
# Inputs:
#   terra_dens_list: List of terra SpatRaster density objects
#   dem: Digital Elevation Model as terra SpatRaster
#   rivers: Vector object with rivers
#   lakes: Vector object with lakes
#   dem_pal: Color palette for DEM
#   dens_pal: Color palette for density surfaces
#   spat_window: Spatial window for the plots
#   buffer_level: Buffer distance around the spatial window
#   periods: Vector of period names
#   letters_to_use: Vector of letters for labeling panels
#   output_dir: Directory for saving the output
#   output_name: Base name for the output file
#
# Returns:
#   No return value; saves a PNG file with the multi-panel plot
# Notes: in the future, make the number of columns not hardcoded
# same with the output dimensions and resolutions
plot_density_maps <- function(terra_dens_list, dem, rivers, study_area = NULL, lakes, dem_pal, dens_pal,
  spat_window, buffer_level, periods, letters_to_use, output_dir, output_name) {
  n_periods <- length(periods)
  n_cols <- 3
  n_rows <- ceiling(n_periods / n_cols)

  png(
    file = here::here(output_dir, paste0(output_name, ".png")),
    width = 1000 * n_cols,
    height = 1000 * n_rows,
    res = 300
  )
  on.exit(dev.off())

  layout(matrix(1:(n_rows * n_cols), n_rows, n_cols, byrow = TRUE))

  pb <- txtProgressBar(min = 0, max = length(periods), style = 3)
  for (i in seq_along(periods)) {
    period_name <- periods[i]
    letter <- letters_to_use[i]

    plot_dens_with_dem(
      kde = terra_dens_list[[i]],
      dem = dem,
      rivers = rivers,
      lakes = lakes,
      study_area = study_area,
      dem_pal = dem_pal,
      dens_pal = dens_pal,
      zoom_ext = spat_window,
      buffer_level = buffer_level,
      plot_title = paste0(letter, ". ", period_name)
    )
    print(paste("Plotting", periods[i], ".."))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n=====DONE=====\n")
}
