# Initialize renv for reproducibility
install.packages("renv")
renv::init()

# Install required packages
install.packages(c("survival", "ggplot2", "dplyr", "broom", "survminer"))

# Snapshot the environment
renv::snapshot()