# Initialize renv for reproducibility
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Optionally set a CRAN mirror to avoid prompts in non-interactive sessions
# options(repos = c(CRAN = "https://cloud.r-project.org"))

renv::init()

# List of required packages with their purposes
pkg_list <- c(
  "survival",    # Core survival analysis functions
  "survminer",   # Enhanced visualization for survival analysis
  "ggplot2",     # Data visualization
  "dplyr",       # Data manipulation
  "broom",       # Model tidying
  "splines"      # Spline transformations for modeling non-linear effects
)

# Install missing packages
missing_pkgs <- pkg_list[!pkg_list %in% installed.packages()[,"Package"]]
if (length(missing_pkgs) > 0) {
  message("Installing the following packages: ", paste(missing_pkgs, collapse=", "))
  install.packages(missing_pkgs)
}

# Verify all packages can be loaded
lapply(pkg_list, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "failed to install or load properly"))
  } else {
    message(paste("Package", pkg, "is successfully installed"))
  }
})

# Snapshot the environment
renv::snapshot()

message("Setup completed successfully!")