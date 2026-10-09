#' Linear normalization: 1~10
#' @param x A value to be normalized
#' @param x_max Queried maximum in historical values
#' @param x_min Queried minimum in historical values
#' @param decreasing if True, a decreased linear normalization will be applied.
#' @export
qc_linear_norm <- function(x, x_min, x_max, decreasing = F) {
  if (length(x) == 1) {
    if (x <= x_min) {
      x <- x_min
    } else if (x >= x_max) {
      x <- x_max
    }
  }
  if (decreasing == F) {
    x_norm <- sapply(x, function(a) 1 + (a - x_min) * 9 / (x_max - x_min))
  } else {
    x_norm <- sapply(x, function(a) 10 - (a - x_min) * 9 / (x_max - x_min))
  }
  
  return(x_norm <- round(x_norm, digits = 3))
}

#' Label the performance: bad, fair, good, great
#' @param x The value of the metric to be labeled
#' @param x_ref The reference values of the metric
#' @param cutoff The cut-offs of the labeling
#' @import stats
#' @importFrom dplyr between
#' @export
qc_performance <- function(x, x_ref, cutoff = c(0, 0.2, 0.5, 0.8, 1)) {
  ref_perc <- quantile(x_ref, cutoff, na.rm = T)
  if (between(x, ref_perc[1], ref_perc[2])) {
    x_class <- "Bad"
  } else if (between(x, ref_perc[2], ref_perc[3])) {
    x_class <- "Fair"
  } else if (between(x, ref_perc[3], ref_perc[4])) {
    x_class <- "Good"
  } else if (between(x, ref_perc[4], ref_perc[5])) {
    x_class <- "Great"
  }
  
  return(x_class)
}

#' Rank in historical performances
#' @param x The value of the metric to be labeled
#' @param x_hist The reference values of the metric
#' @export
qc_rank <- function(x, x_hist) {
  x_all <- c(x, x_hist[!is.na(x_hist)])
  x_pos <- floor(rank(-x_all)[1])
  x_rank <- c(paste(x_pos, "/", length(x_all), sep = ""))
  
  return(x_rank)
}

#' Calculating all QC metrics
#' @param pro_dt A expression table file (at protein level)
#' @param pep_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory for results
#' @param plot if True, a plot will be output.
#' @export
qc_allmetrics <- function(pro_dt, meta_dt, pep_dt = NULL,
                          output_dir = NULL, plot = TRUE) {
  # Basic information ---------------------------
  pro_info <- qc_info(pro_dt, meta_dt)
  
  # SNR -----------------------------------------
  snr_results <- qc_snr(pro_dt, meta_dt, output_dir, plot)
  snr_value <- snr_results$SNR
  
  # 3. Recall (自动使用内置 reference_dataset_quali)
  if (!is.null(pep_dt)) {
    recall_value <- qc_recall(pep_dt, meta_dt = meta_dt)
  } else {
    recall_value <- NA
  }
  
  # 4. PCC (自动使用内置 reference_dataset_quant)
  if (!is.null(pep_dt)) {
    cor_results <- qc_cor(pep_dt, meta_dt, output_dir, plot)
    cor_value <- cor_results$COR
  } else {
    cor_results <- NULL
    cor_value <- NA
  }
  
  # QC results ----------------------------------
  metrics <- c(
    "Number of features",
    "Missing percentage (%)",
    "Absolute Correlation",
    "Coefficient of variantion (CV, %)",
    "Recall",
    "Signal-to-Noise Ratio (SNR)",
    "PCC"
  )
  # qc_values <- c(pro_info, snr_value, cor_value)
  qc_values <- c(pro_info, recall_value, snr_value, cor_value)
  
  # Output --------------------------------------
  output_table <- data.table(
    "Quality Metrics" = metrics,
    "Value" = qc_values
  )
  
  all_results <- list(
    snr_results = snr_results,
    cor_results = cor_results,
    output_table = output_table
  )
  return(all_results)
}


#' Calculating: Total Score
#' @param allmetrics_dt the output table from qc_allmetrics
#' @param ref_qc historical data set
#' @param ref_qc_norm historical data set
#' @param ref_qc_stat historical data set
#' @param normalized if True, the qc values will be linearly normalized to 1~10.
#' @export
qc_total <- function(allmetrics_dt,
                     ref_qc, ref_qc_norm, ref_qc_stat, normalized = T) {
  # Only align metric names; this does not recompute historical baselines.
  legacy_pcc_name <- "Relative Correlation with Reference Datasets (RC)"
  ref_qc <- as.data.frame(ref_qc, check.names = FALSE)
  ref_qc_norm <- as.data.frame(ref_qc_norm, check.names = FALSE)
  ref_qc_stat <- as.data.frame(ref_qc_stat, check.names = FALSE)
  if (!"PCC" %in% names(ref_qc)) {
    names(ref_qc)[names(ref_qc) == legacy_pcc_name] <- "PCC"
  }
  if (!"PCC" %in% names(ref_qc_norm)) {
    names(ref_qc_norm)[names(ref_qc_norm) == legacy_pcc_name] <- "PCC"
  }
  if (!"PCC" %in% ref_qc_stat$`Quality Metrics`) {
    idx <- which(ref_qc_stat$`Quality Metrics` == legacy_pcc_name)
    ref_qc_stat$`Quality Metrics` <- as.character(ref_qc_stat$`Quality Metrics`)
    ref_qc_stat$`Quality Metrics`[idx] <- "PCC"
  }
  
  # Normalize & Rank: All metrics ----------------------
  output_class <- c()
  output_norm <- c()
  output_rank <- c()
  metrics <- allmetrics_dt$`Quality Metrics`
  for (m in metrics) {
    x <- allmetrics_dt$Value[allmetrics_dt$`Quality Metrics` %in% m]
    
    # --- 修改开始：增加检查，看该指标是否存在于历史数据中 ---
    # 检查 m 是否在 ref_qc 的列名中
    if (m %in% colnames(ref_qc)) {
      # 如果历史数据里有这个指标 (比如 SNR, PCC)，正常计算
      x_ref_norm <- as.numeric(ref_qc_norm[, colnames(ref_qc_norm) %in% m])
      x_ref <- ref_qc[, colnames(ref_qc) %in% m]
      
      if (!is.na(x)) {
        x_max <- max(x_ref, na.rm = T)
        x_min <- min(x_ref, na.rm = T)
        if (m %in% c(
          "Coefficient of variantion (CV, %)",
          "Missing percentage (%)"
        )) {
          x_norm <- qc_linear_norm(x, x_min, x_max, decreasing = T)
        } else {
          x_norm <- qc_linear_norm(x, x_min, x_max)
        }
        x_rank <- qc_rank(x_norm, x_ref_norm)
        x_class <- qc_performance(x_norm, x_ref_norm)
      } else {
        x_norm <- NA
        x_rank <- NA
        x_class <- NA
      }
    } else {
      # --- 如果历史数据里没有这个指标 (比如 Recall) ---
      # 跳过排名计算，赋予默认值，防止报错
      x_norm <- NA    # 无法归一化
      x_rank <- "-"   # 无法排名
      x_class <- "-"  # 无法评级 (Bad/Good)
    }
    # --- 修改结束 ---
    
    output_norm <- c(output_norm, x_norm)
    output_rank <- c(output_rank, x_rank)
    output_class <- c(output_class, x_class)
  }
  
  # Normalize & Rank: Total score ----------------------
  # 计算 Total Score 时，需要排除掉那些没有归一化值的指标 (Recall)
  valid_norms <- as.numeric(output_norm[!is.na(output_norm)])
  
  total_ref <- as.numeric(ref_qc_norm$Total)
  
  # 使用几何平均数计算 Total Score (只基于 valid_norms)
  if (length(valid_norms) > 0) {
    total_value <- round(geometric.mean(valid_norms), 3)
    total_norm <- qc_linear_norm(total_value, min(total_ref), max(total_ref))
  } else {
    total_value <- 0
    total_norm <- 0
  }
  
  total_ref_norm <- as.numeric(ref_qc_norm$Total_norm)
  total_rank <- qc_rank(total_norm, total_ref_norm)
  total_c <- qc_performance(total_norm, total_ref_norm)
  
  # Output ---------------------------------------------
  allmetrics <- c(metrics, "Total", "Total_norm")
  allnorm_dt <- data.table(
    "Quality Metrics" = allmetrics,
    "Value" = c(output_norm, total_value, total_norm)
  )
  allmetrics_dt <- rbind(
    allmetrics_dt,
    data.table(
      "Quality Metrics" = "Total Score",
      "Value" = total_norm
    )
  )
  # Merge 时使用 all.x = TRUE 防止因为指标缺失导致行丢失
  output_table <- merge(
    allmetrics_dt,
    ref_qc_stat,
    by = "Quality Metrics",
    all.x = TRUE 
  )
  
  # 因为 merge 可能会打乱顺序，或者有了新指标，简单的行索引 (output_table[c(4,3...)]) 可能会出错
  # 建议根据 metrics 的顺序重新排列，或者保持 merge 后的顺序
  # 这里为了保持你的代码风格，我们只修复报错
  
  # 重新构建最终表格
  # 注意：Rank 和 Performance 的长度现在和 output_table 的行数对应
  # 我们通过匹配 Quality Metrics 来赋值，而不是假定行号
  
  output_table$Rank <- NA
  output_table$Performance <- NA
  
  # 填入计算好的 Rank 和 Performance (除 Total 外)
  for(i in 1:length(metrics)) {
    row_idx <- which(output_table$`Quality Metrics` == metrics[i])
    if(length(row_idx) > 0) {
      output_table$Rank[row_idx] <- output_rank[i]
      output_table$Performance[row_idx] <- output_class[i]
    }
  }
  
  # 填入 Total Score 的 Rank 和 Performance
  total_row <- which(output_table$`Quality Metrics` == "Total Score")
  if(length(total_row) > 0) {
    output_table$Rank[total_row] <- total_rank
    output_table$Performance[total_row] <- total_c
  }
  
  output_list <- list(
    Normalized = allnorm_dt,
    Raw = output_table
  )
  return(output_list)
}
