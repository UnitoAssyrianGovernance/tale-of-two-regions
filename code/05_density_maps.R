## Calculate running correlation between paleoclimate data and archaeodemographic proxies
## Last edited 2025-09-22
## Authors: Andrea Titolo, Alessio Palmisano
## Packages version at the time of analysis
## here: 1.0.1
# sf: 1.0-21
# spatstat: 3.4-0
# terra: 1.8-60
# rnaturalearth: 1.1.0
# scales: 1.4.0
# Aim of this script is to carry out point-pattern analysis on settlement data
required_packages <- c("sf", "here", "spatstat", "terra", "rnaturalearth", "scales")

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

# Source functions scripts
source(here::here("code/functions/calculate_bp_dates.R"))
source(here::here("code/functions/filter_sites_by_period.R"))
source(here::here("code/functions/spatial_kde_functions.R"))


# Data -----------------------------------------------------------------------

## Local Data -----------------------------------------------------------------------

# Load archaeological data
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

sites <- st_read(dsn = here("data/dataset.gpkg"), layer = "archaeo_sites")

study_area <- st_read(
  dsn = here("data/context_dataset.gpkg"),
  layer = "context_project_boundaries_buffer"
)

study_area <- study_area[study_area$name %in% c("Samaria", "Judah"), ]

# Raster data
raster_file <- here::here("data/raw/raster/dem_levant_clip.tif")

if (!file.exists(raster_file)) {
  # Create directory if it doesn't exist
  dir.create(dirname(raster_file), recursive = TRUE, showWarnings = FALSE)
  file_link <- "https://ndownloader.figshare.com/files/58300900"
  browseURL(file_link)
  cat(
    "==========\n",
    "Browser download page opened\n",
    "Please download 'dem_levant_clip.tif' and place it in: ", dirname(raster_file), "\n",
    "Then re-run the script\n",
    "If the browser does not open, visit: ", file_link, "\n",
    "==========\n"
  )
  stop("DEM file required - please download and re-run")

} else if(file.exists(raster_file)) {
  dem <- terra::rast(raster_file)
  cat("==========\nDEM file found and loaded\n==========\n")
}

## Online data  -----------------------------------------------------------------------

# Download data from natural earth
countries <- rnaturalearth::ne_download(
  scale = 10,
  type = "admin_0_countries",
  category = "cultural"
)

# select countries, egypt, israel, jordan, lebanon, palestine, syria
countries <- countries[
  countries$NAME_EN %in%
    c("Egypt", "Israel", "Jordan", "Lebanon", "Palestine", "Syria"),
]

countries <- st_make_valid(countries)

countries <- st_union(countries)

# Download rivers and lakes
rivers <- rnaturalearth::ne_download(
  scale = 10,
  type = "rivers_lake_centerlines",
  category = "physical"
)

lakes <- rnaturalearth::ne_download(
  scale = 10,
  type = "lakes",
  category = "physical"
)

# Data Processing -----------------------------------------------------------------------

# Mask dem with countries for plotting without the sea and other artifacts
dem <- terra::mask(dem, vect(countries))

# Make sure we are only using sites from our study area
sites <- sites |> st_intersection(study_area)

# Keep only a handful of colums
sites <- sites[, c("SiteID", "Name", "SizeHa", "Region", "StartDate", "EndDate", "Period")]

sites <- get_bp_dates(
  df = sites,
  start_date = sites$StartDate,
  end_date = sites$EndDate,
  duration = TRUE
)

# Refactor sites_list with new periods
sites_list <- list()
sites_list <- list(
  "Chalcolithic" = filter_sites_by_period(sites, "Chalcolithic", remove_na_size = FALSE),
  "Early Bronze Age I" = filter_sites_by_period(
    sites,
    "Early Bronze Age I[AB]?(?:-II)?$",
    remove_na_size = FALSE
  ),
  "Early Bronze Age II-III" = filter_sites_by_period(
    sites,
    "Early Bronze Age (?:I-)?II(?:-III)?$",
    remove_na_size = FALSE
  ),
  "Int. Bronze Age" = filter_sites_by_period(
    sites,
    "Early Bronze Age IV/Middle Bronze Age I/Int. Bronze",
    remove_na_size = FALSE
  ),
  "Middle Bronze Age" = filter_sites_by_period(
    sites,
    "Middle Bronze",
    exclude_pattern = "Int. Bronze",
    remove_na_size = FALSE
  ),
  "Late Bronze Age" = filter_sites_by_period(sites, "Late Bronze", remove_na_size = FALSE),
  "Iron Age I" = filter_sites_by_period(
    sites,
    paste0(c("Iron Age I$", "Iron Age I-", "Iron Age$"), collapse = "|"),
    remove_na_size = FALSE
  ),
  "Iron Age II" = filter_sites_by_period(
    sites,
    paste0(c("Iron Age II$", "Iron Age II-", "Iron Age II", "Iron Age I-II"), collapse = "|"),
    exclude_pattern = paste0(
      c("Iron Age III-", "Iron Age III \\(Persian\\)", "Iron Age$"),
      collapse = "|"
    ),
    remove_na_size = FALSE
  ),
  "Iron Age IIa" = filter_sites_by_period(
    sites,
    paste0(c("\\bIron Age IIa\\b", "Iron Age IIa-", "Iron Age I-IIa", "Iron Age$"), collapse = "|"),
    remove_na_size = FALSE
  ),
  "Iron Age IIb" = filter_sites_by_period(
    sites,
    paste0(
      c(
        "\\bIron Age IIb\\b",
        "Iron Age IIb-",
        "Iron Age I-IIb",
        "Iron Age IIa-b",
        "Iron Age IIa-c",
        "Iron Age$"
      ),
      collapse = "|"
    ),
    remove_na_size = FALSE
  ),
  "Iron Age IIc" = filter_sites_by_period(
    sites,
    paste0(
      c(
        "\\bIron Age IIc\\b",
        "\\bIron Age I-IIc\\b",
        "Iron Age IIc-",
        "Iron Age II[a-b]?-c",
        "Iron Age$"
      ),
      collapse = "|"
    ),
    remove_na_size = FALSE
  ),
  "Iron Age III" = filter_sites_by_period(
    sites,
    paste0(
      c("Iron Age III-", "Iron Age III \\(Persian\\)", "Iron Age II[a-c]?-III", "Iron Age$"),
      collapse = "|"
    ),
    remove_na_size = FALSE
  )
)

# Count the number of sites per period from the list above
sites_per_period <- lapply(sites_list, function(x) nrow(x))

# Remove Iron Age II subperiods and Hellenistic
sites_list <- sites_list[
  -which(
    names(sites_list) == "Iron Age IIa" |
      names(sites_list) == "Iron Age IIb" |
      names(sites_list) == "Iron Age IIc"
  )
]

plot_titles <- c(
  "Chalcolithic (4500-3800 BCE)",
  "Early Bronze Age I (3800-3050 BCE)",
  "Early Bronze Age II-III (3050-2500 BCE)",
  "Int. Bronze Age (2500-2000 BCE)",
  "Middle Bronze Age (2000-1550 BCE)",
  "Late Bronze Age (1550-1150 BCE)",
  "Iron Age I (1150-980 BCE)",
  "Iron Age II (980-539 BCE)",
  "Iron Age III (539-333 BCE)"
)

# Spatial KDE -----------------------------------------------------------------------

# If necessary, add more bandwidth to the list to experiment
# Uncomment the output name below to programmatically name the figures for each bandwidth
bw_list <- c(10000)
weighted_kde <- TRUE
output_name <- "fig-09"
# output_name <- paste0("site_density", bw_list[j],"w")

for (j in seq_along(bw_list)) {

  if (weighted_kde) {
    output_name <- output_name
    weight  <- TRUE
    weight_value <- "SizeHa"
    plot_title <- paste0(plot_titles, "-weighted")
  } else if (weighted_kde == FALSE) {
    output_name <- output_name
    weight  <- FALSE
    weight_value <- NULL
    plot_title <- plot_titles
  }

    dens_list <- generate_kde_maps(
    sites_list,
    spat_window = study_area,
    bandwidth = bw_list[j],
    cell_res = 100,
    periods = names(sites_list),
    weight = weight,
    weight_value = weight_value
  )

  terra_dens_list <- convert_kde_to_terra(
    dens_list,
    crs_source = dem,
    spat_window = study_area,
    periods = names(sites_list)
  )

  plot_density_maps(
    terra_dens_list,
    dem = dem,
    rivers = rivers,
    lakes = lakes,
    study_area = study_area,
    dem_pal = "grey",
    dens_pal = "plasma",
    spat_window = study_area,
    buffer_level = 8500,
    periods = plot_title,
    letters_to_use = letters[seq(from = 1, length.out = length(names(sites_list)))],
    output_dir = "output/png",
    output_name = output_name
  )
}
