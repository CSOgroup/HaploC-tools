install_if_needed <- function(package, version = NULL) {
    if (!require(package, character.only = TRUE)) {
        if (!is.null(version)) {
            package <- paste0(package, "_", version)
        }
        install.packages(package, dependencies = TRUE)
    }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

cran_packages <-
    c(
        "R.utils",
        "doParallel",
        "ape",
        "dendextend",
        "fitdistrplus",
        "Matrix",
        "rARPACK",
        "factoextra",
        "data.table",
        "fields",
        "ggplot2",
        "strawr"
    )

## Install igraph with conda as it requires specific libraries
system("conda install conda-forge::r-igraph")

bioconductor_packages <- c("GenomicRanges")

for (pkg in cran_packages) {
    install_if_needed(pkg)
}

for (pkg in bioconductor_packages) {
    BiocManager::install(pkg)
}

# Install calder
install.packages("./HaploC-tools/CALDER2/", repos = NULL, type = "source")

## Check

pkgs <- c(cran_packages, bioconductor_packages, "CALDER")

check_package <- function(package) {
    if (!require(package, character.only = TRUE)) {
        cat(sprintf("Package '%s' is not installed.\n", package))
        return(FALSE)
    } else {
        cat(sprintf("Package '%s' is installed. ", package))
        success <- require(package, character.only = TRUE)
        if (success) {
            cat("Loading successful.\n")
            return(TRUE)
        } else {
            cat("Loading failed.\n")
            return(FALSE)
        }
    }
}

# Check each package
results <- sapply(pkgs, check_package)

# Print summary
cat("\nSummary of package checks:\n")
print(results)

if (all(results)) {
    print("All packages correctly installed!")
} else {
    print("Some packages are not correctly installed. Please check the output above.")
}