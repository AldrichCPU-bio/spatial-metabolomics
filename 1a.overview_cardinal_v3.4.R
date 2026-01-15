# Load required libraries with error handling

#######################################################################################
###### CAUTION: this is an old script, of which Cardinal version <= 3.6 but > 2.xx ####
### if your `Cardinal` version >= 3.8, please refer to `1b.overview_cardinal_v3.8.R` ##
#######################################################################################

load_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it first."))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    futile.logger::flog.info(paste("Loaded package:", pkg))
  }
}

# Define packages
required_packages <- c(
  "fs", "futile.logger", "configr", "stringr",
  "glue", "tidyverse", "dplyr", "Cardinal"
)

# Load packages
load_packages(required_packages)

# Configuration
setup_working_directory <- function(workdir) {
  futile.logger::flog.info(paste("Setting up working directory:", workdir))
  fs::dir_create(workdir)
  setwd(workdir)
  return(workdir)
}

# Read and preprocess imzML data
#######################################################################################
####### CAUTION: importing the `.imzML` file need writable permission, ################
################# please confirm you have the related rights ##########################
#######################################################################################

read_imzml_data <- function(sample_name, data_dir) {
  futile.logger::flog.info(paste("Reading imzML data for sample:", sample_name))
  
  neg_path <- file.path(data_dir, paste0(sample_name, "_neg"))
  pos_path <- file.path(data_dir, paste0(sample_name, "_pos"))
  
  # Check if files exist
  if (!file.exists(paste0(neg_path, ".imzML"))) {
    stop(paste("Negative mode file not found for sample:", sample_name))
  }
  if (!file.exists(paste0(pos_path, ".imzML"))) {
    stop(paste("Positive mode file not found for sample:", sample_name))
  }
  
  # Read negative mode data
  neg_obj <- tryCatch({
    Cardinal::readImzML(name = paste0(sample_name, "_neg"), 
                        folder = data_dir, 
                        units = 'ppm')
  }, error = function(e) {
    stop(paste("Error reading negative mode data for", sample_name, ":", e$message))
  })
  
  # Read positive mode data
  pos_obj <- tryCatch({
    Cardinal::readImzML(name = paste0(sample_name, "_pos"), 
                        folder = data_dir, 
                        units = "ppm")
  }, error = function(e) {
    stop(paste("Error reading positive mode data for", sample_name, ":", e$message))
  })
 
  pData(msi_maldi) = Cardinal::PositionDataFrame(
    coord = Cardinal::coord(msi_maldi)[, 1:2], run = run_info, 
    row.names = rownames(Cardinal::coord(msi_maldi))
  )
  
  # Subset pixels with non-zero sum
  neg_obj <- subsetPixels(neg_obj, pixelApply(neg_obj, sum) > 0)
  pos_obj <- subsetPixels(pos_obj, pixelApply(pos_obj, sum) > 0)
  
  return(list(neg = neg_obj, pos = pos_obj))
}

# Perform PCA analysis
perform_pca_analysis <- function(obj, sample_name, modality, workdir) {
  futile.logger::flog.info(paste("Performing PCA analysis for", 
                                 sample_name, modality, "mode"))
  
  # Create directory for PCA results
  pca_dir <- file.path(workdir, sample_name, "pca", modality)
  fs::dir_create(pca_dir)
  
  # Calculate mean spectrum
  summ_mean <- summarizeFeatures(obj, "mean")
  
  # Plot mean intensity
  mean_intensity_plot <- file.path(pca_dir, paste0(modality, "_mean_intensity.pdf"))
  pdf(mean_intensity_plot, width = 5, height = 3)
  plot(summ_mean)
  dev.off()
  
  # Calculate TIC normalization
  norm_tic <- summarizePixels(obj, c(tic = "sum"))
  tic_plot <- file.path(pca_dir, paste0(modality, "_norm_tic.pdf"))
  pdf(tic_plot, width = 5, height = 3)
  image(norm_tic)
  dev.off()
  
  # Peak picking and alignment
  peak_ref <- summ_mean %>%
    peakPick(SNR = 3) %>%
    peakAlign(ref = "mean", tolerance = 0.5, units = "mz") %>%
    peakFilter(freq.min = 0.1) %>%
    process()
  
  # Peak binning
  peak_peaks <- obj %>%
    normalize(method = "tic") %>%
    peakBin(ref = mz(peak_ref), tolerance = 0.5, units = "mz") %>%
    process()
  
  # Perform PCA
  pca_result <- PCA(peak_peaks, ncomp = 4)
  
  # Plot PCA results
  pca_plot <- file.path(pca_dir, paste0(modality, "_pca2.pdf"))
  pdf(pca_plot, width = 5, height = 4)
  image(pca_result, contrast.enhance = "histogram", normalize.image = "linear")
  plot(pca_result, lwd = 2)
  dev.off()
  
  # Save results
  write_rds(peak_peaks, file.path(pca_dir, "peak_peaks.rds"))
  write_rds(pca_result, file.path(pca_dir, "pca.rds"))
  
  return(list(pca_result = pca_result, peak_peaks = peak_peaks))
}

# Perform spatial shrunken centroids analysis
perform_spatial_shrunken_analysis <- function(obj, sample_name, modality, workdir) {
  futile.logger::flog.info(paste("Performing spatial shrunken centroids analysis for",
                                 sample_name, modality, "mode"))
  
  # Create directory for spatial shrunken results
  ssc_dir <- file.path(workdir, sample_name, "spatial_shrunken", modality)
  fs::dir_create(ssc_dir)
  
  # Set seed for reproducibility
  set.seed(1)
  
  # Perform spatial shrunken centroids
  ssc_result <- spatialShrunkenCentroids(obj, method = "adaptive",
                                         r = 1, s = c(1, 2, 4, 8), k = 6)
  
  # Plot spatial shrunken centroids
  ssc_plot <- file.path(ssc_dir, "spatial_shrunken_centroids.pdf")
  pdf(ssc_plot, width = 8, height = 8)
  image(ssc_result, model = list(s = c(1, 2, 4, 8)))
  plot(ssc_result, model = list(s = c(1, 2, 4, 8)), lwd = 2)
  dev.off()
  
  # Plot shrunken mean spectra
  mean_spectra_plot <- file.path(ssc_dir, "spatial_shrunken_centroids_mean_spectra.pdf")
  pdf(mean_spectra_plot, width = 8, height = 10)
  plot(ssc_result, model = list(s = c(1, 2, 4, 8)), lwd = 2)
  dev.off()
  
  # Plot split segments
  split_plot <- file.path(ssc_dir, "spatial_shrunken_centroids_mean_spectra_split.pdf")
  pdf(split_plot, width = 8, height = 10)
  cols <- discrete.colors(6)
  for (col in 1:6) {
    plot(ssc_result, model = list(s = c(1, 2, 4, 8)), 
         lwd = 2, column = col, col = cols[col])
  }
  dev.off()
  
  # Extract top features
  top_features <- lapply(1:6, function(id) {
    feature_list <- list()
    for (idx in 1:4) {
      sp_value <- c(1, 2, 4, 8)[idx]
      feature_list[[idx]] <- topFeatures(ssc_result, 
                                         model = list(s = sp_value), 
                                         class == id, 
                                         n = Inf) %>% 
        as.data.frame()
    }
    dplyr::bind_rows(feature_list)
  }) %>% dplyr::bind_rows()
  
  # Save results
  write_tsv(top_features, file.path(ssc_dir, "top_features.tsv"))
  write_rds(ssc_result, file.path(workdir, sample_name, 
                                  paste0(modality, "_cardinal_res.rds")))
  
  return(list(ssc_result = ssc_result, top_features = top_features))
}

# Process a single sample
process_sample <- function(sample_name, project, dataset, species, workdir) {
  futile.logger::flog.info(paste("Processing sample:", sample_name))
  
  # Create sample directory
  sample_dir <- file.path(workdir, sample_name)
  fs::dir_create(sample_dir)
  
  # Define data directory
  data_dir <- glue("~/projects/{project}/data/{dataset}/{species}/metabolism/{sample_name}")
  
  # Read imzML data
  obj_list <- read_imzml_data(sample_name, data_dir)
  
  # Process each modality
  modalities <- c("neg", "pos")
  results <- list()
  
  for (modality in modalities) {
    futile.logger::flog.info(paste("Processing", modality, "mode for", sample_name))
    
    # Perform PCA analysis
    pca_results <- perform_pca_analysis(obj_list[[modality]], 
                                        sample_name, 
                                        modality, 
                                        workdir)
    
    # Perform spatial shrunken centroids analysis
    ssc_results <- perform_spatial_shrunken_analysis(obj_list[[modality]], 
                                                     sample_name, 
                                                     modality, 
                                                     workdir)
    
    results[[modality]] <- list(pca = pca_results, ssc = ssc_results)
  }
  
  futile.logger::flog.info(paste("Completed processing for sample:", sample_name))
  return(results)
}

# Main execution function
main <- function() {
  # Configuration
  workdir <- "your path of working directory"
  project <- "your_project_name"
  dataset <- "your_dataset_name"
  species <- "your_species_name"
  
  # List of samples to process
  samples <- c("sample1", "sample2", "sample3", 
               "sample4", "sample5", "sample6")
  
  # Set up working directory
  workdir <- setup_working_directory(workdir)
  
  # Process each sample
  all_results <- lapply(samples, function(sample) {
    tryCatch({
      process_sample(sample, project, dataset, species, workdir)
    }, error = function(e) {
      futile.logger::flog.error(paste("Error processing sample", sample, ":", e$message))
      return(NULL)
    })
  })
  
  # Name the results list
  names(all_results) <- samples
  
  # Remove NULL results (failed samples)
  all_results <- all_results[!sapply(all_results, is.null)]
  
  # Summary
  futile.logger::flog.info(paste("Processing completed. Successfully processed", 
                                 length(all_results), "out of", length(samples), "samples."))
  
  return(all_results)
}

# Execute the main function if script is run directly
if (!interactive() && sys.nframe() == 0) {
  results <- main()
}
