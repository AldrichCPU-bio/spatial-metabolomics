# ===================== 加载必要的包 =====================
load_required_packages <- function() {
  required_packages <- c("Cardinal", "BiocParallel", "glue")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it first."))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message(paste("✓ Loaded package:", pkg))
  }
}

# ===================== 可视化辅助函数 =====================
create_visualization <- function(data, filename, title, 
                                 mz_low = NULL, mz_high = NULL, 
                                 i = 1, linewidth = 0.5, 
                                 width = 8, height = 3) {
  # 创建质量窗口图
  pdf(filename, width = width, height = height)
  
  if (!is.null(mz_low) && !is.null(mz_high)) {
    # 如果指定了mz范围，则进行子集化
    subset_data <- Cardinal::subsetFeatures(data, mz < mz_high, mz > mz_low)
    plot(subset_data, i = i, linewidth = linewidth, main = title)
  } else {
    # 如果没有指定mz范围，直接绘图
    plot(data, i = i, linewidth = linewidth, main = title)
  }
  
  dev.off()
  message(paste("✓ Created visualization:", filename))
}

# ===================== 数据预处理函数 =====================
preprocess_ms_data <- function(mse_data, sample_id, workers = 8) {
  message(paste("\n🔬 Processing sample:", sample_id))
  message(paste("   Using", workers, "workers for parallel processing"))
  
  # 1. 平滑处理
  message("   Step 1: Smoothing data...")
  mse_smoothed <- smooth(mse_data, method = "adaptive")
  
  # 可视化平滑结果
  create_visualization(
    data = mse_smoothed,
    filename = glue("smooth_{sample_id}.pdf"),
    title = paste("Smoothed Spectrum -", sample_id),
    mz_low = 198, mz_high = 200,
    i = 16, linewidth = 0.5
  )
  
  # 2. 基线校正
  message("   Step 2: Baseline correction...")
  mse_baseline <- reduceBaseline(mse_smoothed, method = "locmin")
  
  # 可视化基线校正结果
  create_visualization(
    data = mse_baseline,
    filename = glue("baseline_{sample_id}.pdf"),
    title = paste("Baseline Corrected -", sample_id),
    mz_low = 195, mz_high = 200,
    i = 16, linewidth = 0.5
  )
  
  # 3. 处理数据
  message("   Step 3: Processing data...")
  mse_processed <- process(
    mse_baseline, 
    BPPARAM = MulticoreParam(workers = workers, progressbar = TRUE)
  )
  
  # 4. 估计参考峰
  message("   Step 4: Estimating reference peaks...")
  mse_ref <- estimateReferencePeaks(mse_processed)
  
  # 5. 重新校准
  message("   Step 5: Recalibration...")
  mse_recalibrated <- recalibrate(
    mse_processed, 
    ref = mse_ref, 
    method = "locmax", 
    tolerance = 1500, 
    units = "ppm"
  ) %>%
    process(BPPARAM = MulticoreParam(workers = workers, progressbar = TRUE))
  
  # 可视化重新校准结果
  create_visualization(
    data = mse_recalibrated,
    filename = glue("recalibrate_{sample_id}.pdf"),
    title = paste("Recalibrated -", sample_id),
    mz_low = 197, mz_high = 200,
    i = 186:187, linewidth = 0.5
  )
  
  # 6. 峰检测
  message("   Step 6: Peak picking...")
  mse_peakpicked <- peakPick(
    mse_recalibrated, 
    method = "diff", 
    SNR = 3, 
    units = "mz", 
    type = "height"
  ) %>%
    process(BPPARAM = MulticoreParam(workers = workers, progressbar = TRUE))
  
  # 可视化峰检测结果（对比）
  pdf(glue("peakpick_{sample_id}.pdf"), width = 8, height = 6)
  par(mfrow = c(2, 1))
  
  # 原始数据的峰检测结果
  mse_recalibrated %>% 
    Cardinal::subsetFeatures(mz < 200, mz > 197) %>% 
    plot(., i = 16, linewidth = 0.5, 
         main = paste("Before Peak Picking -", sample_id))
  
  # 峰检测后的结果
  mse_peakpicked %>% 
    Cardinal::subsetFeatures(mz < 200, mz > 197) %>% 
    plot(., i = 16, linewidth = 0.5, 
         main = paste("After Peak Picking -", sample_id))
  
  dev.off()
  message(paste("✓ Created peak picking comparison:", glue("peakpick_{sample_id}.pdf")))
  
  # 7. 峰对齐
  message("   Step 7: Peak alignment...")
  mse_aligned <- peakAlign(
    mse_peakpicked, 
    tolerance = 0.05, 
    units = "mz",
    BPPARAM = MulticoreParam(workers = workers, progressbar = TRUE)
  )
  
  # 可视化峰对齐结果
  create_visualization(
    data = mse_aligned,
    filename = glue("peakalign_{sample_id}.pdf"),
    title = paste("Peak Aligned -", sample_id),
    mz_low = 197, mz_high = 200,
    i = c(16, 48), linewidth = 0.5
  )
  
  message(paste("✅ Completed processing for sample:", sample_id))
  return(mse_aligned)
}

# ===================== 强度图分析函数 =====================
analyze_intensity_profiles <- function(mse_data, sample_id) {
  message(paste("\n📊 Analyzing intensity profiles for:", sample_id))
  
  # 使用两种不同的峰检测方法
  p1 <- peakPick(mse_data, method = "diff", SNR = 3) |>
    plot(i = c(1, 100), linewidth = 2, 
         main = paste("Derivative-based SNR -", sample_id))
  
  p2 <- peakPick(mse_data, method = "filter", SNR = 3) |>
    plot(i = c(1, 100), linewidth = 2, 
         main = paste("Dynamic filtering-based SNR -", sample_id))
  
  # 合并两个图
  combined_plot <- matter::as_facets(
    list(p1, p2), 
    nrow = 2,
    labels = c("Derivative-based SNR", "Dynamic filtering-based SNR")
  )
  
  return(combined_plot)
}

# ===================== 批量处理函数 =====================
batch_process_ms_data <- function(mse_list, workers = 8) {
  message("🚀 Starting batch processing of MS data")
  message(paste("Number of samples:", length(mse_list)))
  
  # 处理每个样本
  processed_results <- lapply(names(mse_list), function(sample_name) {
    mse_data <- mse_list[[sample_name]]
    
    # 获取样本ID
    sample_id <- unique(mse_data$sample_id)
    if (length(sample_id) == 0) {
      sample_id <- sample_name  # 如果没有sample_id，使用列表名
    }
    
    # 预处理数据
    processed_data <- tryCatch({
      preprocess_ms_data(mse_data, sample_id, workers)
    }, error = function(e) {
      warning(paste("Error processing sample", sample_name, ":", e$message))
      return(NULL)
    })
    
    return(processed_data)
  })
  
  # 命名结果列表
  names(processed_results) <- names(mse_list)
  
  # 移除失败的处理
  successful_processing <- !sapply(processed_results, is.null)
  message(paste("✅ Successfully processed:", sum(successful_processing), 
                "out of", length(mse_list), "samples"))
  
  return(processed_results)
}

# ===================== 主执行函数 =====================
run_ms_analysis_pipeline <- function(mse_list, intensity_analysis = TRUE, 
                                     output_dir = ".", workers = 8) {
  
  # 加载必要的包
  load_required_packages()
  
  # 设置工作目录
  original_dir <- getwd()
  setwd(output_dir)
  on.exit(setwd(original_dir))  # 确保函数结束后恢复原目录
  
  message(paste("📁 Output directory:", output_dir))
  
  # 批量处理数据
  processed_data <- batch_process_ms_data(mse_list, workers)
  
  # 如果需要强度图分析
  if (intensity_analysis) {
    message("\n📈 Generating intensity profile analysis")
    
    # 创建强度图
    pdf("mse_intensity_profiles.pdf", width = 10, height = 6)
    
    # 示例：绘制特定样本的强度图
    if ("mse_e6.5_rp1" %in% names(mse_list)) {
      # 原始数据
      mse_list[["mse_e6.5_rp1"]] %>% 
        Cardinal::subsetFeatures(mz < 199, mz > 197) %>% 
        plot(., i = 16, main = "Original Data - mse_e6.5_rp1")
      
      # 处理后的数据（如果可用）
      if (!is.null(processed_data[["mse_e6.5_rp1"]])) {
        processed_data[["mse_e6.5_rp1"]] %>% 
          Cardinal::subsetFeatures(mz < 199, mz > 197) %>% 
          plot(., i = 16, main = "Processed Data - mse_e6.5_rp1")
      }
    }
    
    dev.off()
    message("✓ Created intensity profile analysis: mse_intensity_profiles.pdf")
    
    # 为每个样本生成峰检测方法对比图
    intensity_plots <- lapply(names(mse_list), function(sample_name) {
      mse_data <- mse_list[[sample_name]]
      sample_id <- unique(mse_data$sample_id)
      if (length(sample_id) == 0) sample_id <- sample_name
      
      analyze_intensity_profiles(mse_data, sample_id)
    })
  }
  
  # 返回处理结果
  return(list(
    processed_data = processed_data,
    summary = list(
      total_samples = length(mse_list),
      successful_samples = sum(!sapply(processed_data, is.null)),
      failed_samples = sum(sapply(processed_data, is.null))
    )
  ))
}

# ===================== 使用示例 =====================
example_usage <- function() {
  # 假设 mse_lst 是你的数据列表
  # mse_lst <- list(...)
  
  # 运行分析管道
  results <- run_ms_analysis_pipeline(
    mse_list = mse_lst,
    intensity_analysis = TRUE,
    output_dir = "./ms_analysis_results",
    workers = 8  # 根据你的系统调整
  )
  
  # 打印摘要
  message("\n📋 Analysis Summary:")
  message(paste("  Total samples processed:", results$summary$total_samples))
  message(paste("  Successful:", results$summary$successful_samples))
  message(paste("  Failed:", results$summary$failed_samples))
  
  return(results)
}

# 如果直接运行此脚本，执行示例
if (!interactive() && sys.nframe() == 0) {
  message("⚠️  This is an example. Please modify the function calls as needed.")
  # example_usage()
}
