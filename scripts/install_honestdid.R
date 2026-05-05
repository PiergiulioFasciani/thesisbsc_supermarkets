#!/usr/bin/env Rscript

# ==============================================================================
# Install HonestDiD Package
# ==============================================================================
# HonestDiD is not on CRAN, must be installed from GitHub
# Roth & Rambachan (2023) https://github.com/asheshrambachan/HonestDiD

cat("Checking HonestDiD installation...\n\n")

# Check if HonestDiD is installed
if (requireNamespace("HonestDiD", quietly = TRUE)) {
  cat("✓ HonestDiD is already installed\n")
  packageVersion("HonestDiD")
  cat("\n")
} else {
  cat("HonestDiD not found. Installing from GitHub...\n\n")
  
  # Check if devtools is installed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    cat("Installing devtools first...\n")
    install.packages("devtools", repos = "https://cloud.r-project.org/")
  }
  
  # Install HonestDiD
  cat("Installing HonestDiD from GitHub...\n")
  devtools::install_github("asheshrambachan/HonestDiD")
  
  # Verify installation
  if (requireNamespace("HonestDiD", quietly = TRUE)) {
    cat("\n✓ HonestDiD successfully installed\n")
    cat("Version:", as.character(packageVersion("HonestDiD")), "\n\n")
  } else {
    cat("\n✗ Failed to install HonestDiD\n")
    cat("Please install manually with:\n")
    cat("  devtools::install_github('asheshrambachan/HonestDiD')\n\n")
    quit(status = 1)
  }
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("HonestDiD Ready\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
