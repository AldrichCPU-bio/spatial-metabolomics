#######################################################################
## CAUTION1: this script depends on Cardianl version = 3.8 or higher ##
########################### please check it ###########################
#### CAUTION2: do you have the writable permission `.imzML` file? #####
#######################################################################

# ===================== 自定义主题定义 =====================
create_custom_theme <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      legend.key.size = unit(3, "mm"),
      axis.text = element_text(color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      panel.border = element_rect(linewidth = .5, color = "black", fill = NA)
    )
}

# ===================== 数据读取与预处理 =====================
read_and_preprocess_imzml <- function(imzml_file_path) {
  message(paste("📂 Reading imzML file:", basename(imzml_file_path)))
  
  # 读取imzML数据
  maldi_data <- tryCatch({
    Cardinal::readImzML(imzml_file_path)
  }, error = function(e) {
    stop(paste("Error reading imzML file:", e$message))
  })
  
  # 转换为MSI实验对象
  msi_experiment <- convertMSImagingArrays2Experiment(maldi_data)
  
  # 汇总像素数据
  mse_summary <- summarizePixels(msi_experiment)
  
  # 计算并显示TIC统计
  tic_stats <- pData(mse_summary) %>% 
    as.data.frame() %>% 
    .$TIC %>% 
    summary()
  
  message("📊 TIC Summary:")
  print(tic_stats)
  
  return(mse_summary)
}

# ===================== 基础可视化函数 =====================
visualize_basic_msi_data <- function(mse_object, output_dir = ".") {
  message("🎨 Creating basic visualizations...")
  
  # 1. TIC图像
  tic_plot_path <- file.path(output_dir, "tic_image.pdf")
  pdf(tic_plot_path, width = 8, height = 5)
  image(mse_object, "tic", main = "TIC Image")
  dev.off()
  message(paste("✓ Saved TIC image:", tic_plot_path))
  
  # 2. 特定m/z图像（示例：m/z 198.7 ± 0.2）
  mz_plot_path <- file.path(output_dir, "mz198.7_image.pdf")
  pdf(mz_plot_path, width = 8, height = 5)
  image(
    mse_object,
    mz = 198.7,
    tolerance = 0.2,
    units = "mz",
    main = "m/z 198.7 Image"
  )
  dev.off()
  message(paste("✓ Saved m/z image:", mz_plot_path))
  
  # 3. 质谱图（示例：特定像素的光谱）
  spectra_plot_path <- file.path(output_dir, "spectra_examples.pdf")
  pdf(spectra_plot_path, width = 12, height = 3)
  # 随机选择几个像素进行可视化，或者使用特定的像素索引
  if (nrow(coord(mse_object)) >= 2) {
    plot(
      mse_object,
      i = c(12, 45),  # 示例像素索引
      xlim = c(195, 198),
      main = "Example Spectra"
    )
  } else {
    plot(mse_object, i = 1, main = "Single Pixel Spectrum")
  }
  dev.off()
  message(paste("✓ Saved spectra examples:", spectra_plot_path))
  
  # 4. 坐标分布图
  coord_data <- tibble(
    x = pData(mse_object)$x,
    y = pData(mse_object)$y
  )
  
  density_x <- ggplot(coord_data, aes(x = x)) +
    geom_density() +
    create_custom_theme() +
    labs(x = "X Coordinate", y = "Density")
  
  density_y <- ggplot(coord_data, aes(x = y)) +
    geom_density() +
    create_custom_theme() +
    labs(x = "Y Coordinate", y = "Density")
  
  coord_plot_path <- file.path(output_dir, "coordinate_density.pdf")
  ggsave(coord_plot_path, density_x | density_y, width = 12, height = 4)
  message(paste("✓ Saved coordinate density plots:", coord_plot_path))
}

# ===================== 数据验证与检查 =====================
validate_msi_data <- function(mse_object) {
  message("🔍 Validating MSI data...")
  
  validation_results <- list(
    total_pixels = nrow(coord(mse_object)),
    mz_range = range(mz(mse_object)),
    tic_range = range(pData(mse_object)$tic),
    x_range = range(pData(mse_object)$x),
    y_range = range(pData(mse_object)$y)
  )
  
  message(paste("  Total pixels:", validation_results$total_pixels))
  message(paste("  m/z range:", paste(round(validation_results$mz_range, 2), collapse = " - ")))
  message(paste("  TIC range:", paste(round(validation_results$tic_range, 2), collapse = " - ")))
  message(paste("  X coordinate range:", paste(round(validation_results$x_range, 2), collapse = " - ")))
  message(paste("  Y coordinate range:", paste(round(validation_results$y_range, 2), collapse = " - ")))
  
  return(validation_results)
}

# ===================== 数据保存函数 =====================
save_msi_data <- function(mse_object, output_dir = ".") {
  message("💾 Saving MSI data...")
  
  # 创建RDS目录
  rds_dir <- file.path(output_dir, "rds")
  if (!dir.exists(rds_dir)) {
    dir.create(rds_dir, recursive = TRUE)
  }
  
  # 保存数据
  data_path <- file.path(rds_dir, "msi_data.rds")
  write_rds(mse_object, data_path)
  message(paste("✓ Saved data to:", data_path))
  
  # 可选：保存元数据为CSV
  metadata_path <- file.path(rds_dir, "msi_metadata.csv")
  metadata <- data.frame(
    pixel_id = rownames(coord(mse_object)),
    x = pData(mse_object)$x,
    y = pData(mse_object)$y,
    tic = pData(mse_object)$tic
  )
  write_csv(metadata, metadata_path)
  message(paste("✓ Saved metadata to:", metadata_path))
  
  return(list(data_path = data_path, metadata_path = metadata_path))
}

# ===================== 主分析函数 =====================
analyze_imzml_data <- function(imzml_file_path, output_dir = ".", 
                              create_visualizations = TRUE) {
  message("🚀 Starting imzML data analysis pipeline")
  message(paste("Input file:", basename(imzml_file_path)))
  message(paste("Output directory:", output_dir))
  
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. 读取和预处理数据
  mse_data <- read_and_preprocess_imzml(imzml_file_path)
  
  # 2. 数据验证
  validation_results <- validate_msi_data(mse_data)
  
  # 3. 创建可视化
  if (create_visualizations) {
    visualize_basic_msi_data(mse_data, output_dir)
  }
  
  # 4. 保存数据
  save_paths <- save_msi_data(mse_data, output_dir)
  
  # 返回结果
  results <- list(
    mse_object = mse_data,
    validation = validation_results,
    save_paths = save_paths,
    metadata = list(
      input_file = imzml_file_path,
      output_directory = output_dir,
      timestamp = Sys.time()
    )
  )
  
  message("✅ Analysis completed successfully!")
  return(results)
}

# ===================== 快速检查函数 =====================
quick_check_imzml <- function(imzml_file_path) {
  message(paste("🔎 Quick check of:", basename(imzml_file_path)))
  
  # 仅读取元数据，不加载全部光谱数据
  maldi_data <- tryCatch({
    Cardinal::readImzML(imzml_file_path, as = "MSImagingExperiment")
  }, error = function(e) {
    Cardinal::readImzML(imzml_file_path)
  })
  
  # 基本信息
  cat("IMZML FILE INFORMATION:\n")
  cat("=======================\n")
  cat(paste("File:", basename(imzml_file_path), "\n"))
  cat(paste("Total pixels:", nrow(coord(maldi_data)), "\n"))
  cat(paste("m/z features:", length(mz(maldi_data)), "\n"))
  cat(paste("X range:", paste(range(coord(maldi_data)$x), collapse = " - "), "\n"))
  cat(paste("Y range:", paste(range(coord(maldi_data)$y), collapse = " - "), "\n"))
  cat(paste("Data size:", format(object.size(maldi_data), units = "auto"), "\n"))
  
  return(maldi_data)
}

# ===================== 使用示例 =====================
example_analysis <- function() {
  # 设置输入文件路径和输出目录
  imzml_file <- "path/to/your/imzml_file.imzML"
  output_directory <- "msi_analysis_results"
  
  # 运行分析
  results <- analyze_imzml_data(
    imzml_file_path = imzml_file,
    output_dir = output_directory,
    create_visualizations = TRUE
  )
  
  return(results)
}

# 辅助函数：检查多个imzML文件
check_multiple_imzml_files <- function(imzml_files) {
  file_info <- list()
  
  for (i in seq_along(imzml_files)) {
    imzml_file <- imzml_files[i]
    cat(sprintf("\n[%d/%d] %s\n", i, length(imzml_files), basename(imzml_file)))
    
    tryCatch({
      info <- quick_check_imzml(imzml_file)
      file_info[[basename(imzml_file)]] <- list(
        pixels = nrow(coord(info)),
        mz_features = length(mz(info)),
        x_range = range(coord(info)$x),
        y_range = range(coord(info)$y)
      )
    }, error = function(e) {
      cat(paste("  Error:", e$message, "\n"))
      file_info[[basename(imzml_file)]] <- paste("Error:", e$message)
    })
  }
  
  return(file_info)
}

# 如果直接运行此脚本，执行示例
if (!interactive() && sys.nframe() == 0) {
  # 加载必要包
  required_packages <- c("Cardinal", "ggplot2", "patchwork", "tibble", "readr")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  message("⚠️  This is an example. Please modify file paths as needed.")
  
  # 示例：快速检查文件
  # imzml_file <- "path/to/your/imzml_file.imzML"
  # quick_check_imzml(imzml_file)
  
  # 示例：完整分析
  # results <- example_analysis()
}
