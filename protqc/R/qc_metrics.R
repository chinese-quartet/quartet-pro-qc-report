#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @export
qc_info <- function(expr_dt, meta_dt) {
  # Load data --------------------------------
  m <- meta_dt
  d <- expr_dt
  s <- meta_dt$sample
  # Infer protein scale from the original protein matrix, independently of peptides.
  protein_matrix <- as.matrix(as.data.frame(expr_dt)[, -1L, drop = FALSE])
  storage.mode(protein_matrix) <- "double"
  protein_scale <- pcc_detect_data_scale(protein_matrix, level_label = "Protein CV")
  
  # Replace zero by NA -----------------------
  d[d == 0] <- NA
  
  # Statistics: number of features -----------
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)
  
  # Statistics: missing percentage -----------
  d_all_num <- nrow(d) * (ncol(d) - 1)
  d_missing_num <- length(which(is.na(d)))
  prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)
  
  # Check if replicates available ------------
  samples <- table(s)
  rep_samples <- samples[samples > 1]
  rep_num <- length(rep_samples)
  if (rep_num == 0) {
    stop("No replicates are available.")
  } else {
    # Calculating: absolute correlation ------
    d_mtx <- d[, 2:ncol(d)]
    d_cortest <- corr.test(d_mtx, method = "pearson", adjust = "fdr")
    d_pmtx <- d_cortest$p
    d_cormtx <- d_cortest$r
    d_cormtx[d_pmtx > 0.05] <- 0
    d_cordf <- melt(d_cormtx)
    d_cordf <- d_cordf[d_cordf$Var2 != d_cordf$Var1, ]
    d_cordf <- merge(d_cordf, m, by.x = "Var1", by.y = "library")
    d_cordf <- merge(d_cordf, m, by.x = "Var2", by.y = "library")
    d_cordf <- d_cordf[d_cordf$sample.x == d_cordf$sample.y, ]
    cor_value <- median(d_cordf$value)
    stat_acor <- round(cor_value, 3)
  }
  
  # Calculating: CV --------------------------
  d_long <- melt(d)
  d_long <- na.omit(d_long)
  # CV is calculated on linear protein abundance.
  if (protein_scale == "log2") {
    d_long$value <- 2^(d_long$value)
  }
  d_long <- merge(d_long, m, by.x = "variable", by.y = "library")
  colnames(d_long) <- c("library", "feature", "value", "sample")
  d_cv <- aggregate(
    value ~ feature + sample,
    data = d_long,
    FUN = function(x) sd(x) / mean(x)
  )
  stat_cv <- round(median(d_cv$value, na.rm = T) * 100, 3)
  
  # Output -----------------------------------
  stat_all <- c(stat_num, stat_missing, stat_acor, stat_cv)
  
  return(stat_all)
}

#' Get color mapping for samples
#' @param samples Vector of sample names
#' @return Named vector of colors
#' @export
get_sample_colors <- function(samples) {
  # Define fixed color palette
  color_palette <- c(
    "D5" = "#4CC3D9", # Blue
    "D6" = "#7BC8A4", # Green
    "F7" = "#FFC65D", # Yellow
    "M8" = "#F16745" # Red
  )
  
  # Return colors only for existing samples
  available_colors <- color_palette[samples]
  return(available_colors)
}

#' Calculating SNR value; Plotting a PCA panel
#' @param expr_dt A expression profile (at protein level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @import stats
#' @import utils
#' @importFrom rlang :=
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 ggsave
#' @importFrom ggthemes theme_few
#' @export
qc_snr <- function(expr_dt, meta_dt, output_dir = NULL, plot = TRUE) {
  # Load data --------------------------------------
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])
  # Infer protein scale independently; PCA/SNR always uses log2 abundance.
  protein_matrix <- as.matrix(expr_df)
  storage.mode(protein_matrix) <- "double"
  protein_scale <- pcc_detect_data_scale(protein_matrix, level_label = "Protein SNR")
  expr_df <- pcc_prepare_abundance(
    native = protein_matrix, scale = protein_scale,
    level_label = "Protein SNR", zero_is_missing = TRUE
  )$log2
  
  # Replace NA by zero -----------------------------
  expr_df[is.na(expr_df)] <- 0
  
  # Label the grouping info ------------------------
  ids <- colnames(expr_df)
  # Match groups by the exact library ID, independent of metadata row order.
  if (!all(c("library", "sample") %in% names(meta_dt))) {
    stop("SNR metadata must contain library and sample columns.")
  }
  libraries <- as.character(meta_dt$library)
  if (anyNA(libraries) || any(!nzchar(trimws(libraries))) ||
      anyDuplicated(libraries)) {
    stop("SNR metadata library IDs must be non-missing and unique.")
  }
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("SNR expression sample column names must be non-missing and unique.")
  }
  meta_index <- match(ids, libraries)
  if (anyNA(meta_index)) {
    stop("SNR expression samples missing from metadata$library: ",
         paste(ids[is.na(meta_index)], collapse = ", "))
  }
  group <- as.character(meta_dt$sample[meta_index])
  if (anyNA(group) || any(!nzchar(trimws(group)))) {
    stop("SNR matched sample groups must be non-missing.")
  }
  ids_group_mat <- data.table(id = ids, group = group)
  
  # PCA --------------------------------------------
  expr_df_t <- t(expr_df)
  
  # Remove constant/zero variance features before PCA
  col_vars <- apply(expr_df_t, 2, var, na.rm = TRUE)
  zero_var_cols <- which(col_vars == 0 | is.na(col_vars))
  
  # 记录原始特征数和过滤后特征数
  n_features_original <- ncol(expr_df_t)
  n_features_removed <- length(zero_var_cols)
  
  if (n_features_removed > 0) {
    message(sprintf(
      "Removing %d features with zero variance before PCA analysis.",
      n_features_removed
    ))
    expr_df_t <- expr_df_t[, -zero_var_cols, drop = FALSE]
  }
  
  n_features_used <- ncol(expr_df_t)
  
  # Check if enough features remain
  if (n_features_used < 2) {
    stop("Not enough features with non-zero variance for PCA analysis. At least 2 features required.")
  }
  
  pca_prcomp <- prcomp(expr_df_t, retx = T, scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$sample_id <- rownames(pcs)
  pcs$sample <- group[match(pcs$sample_id, ids)]
  
  # Calculating: SNR -------------------------------
  dt_perc_pcs <- data.table(
    PCX = 1:nrow(pcs),
    Percent = summary(pca_prcomp)$importance[2, ],
    AccumPercent = summary(pca_prcomp)$importance[3, ]
  )
  
  dt_dist <- data.table(
    id_a = rep(ids, each = length(ids)),
    id_b = rep(ids, time = length(ids))
  )
  
  dt_dist$group_a <- ids_group_mat[match(dt_dist$id_a, ids_group_mat$id)]$group
  dt_dist$group_b <- ids_group_mat[match(dt_dist$id_b, ids_group_mat$id)]$group
  
  dt_dist[, type := ifelse(id_a == id_b, "Same",
                           ifelse(group_a == group_b, "Intra", "Inter")
  )]
  
  dt_dist[, dist := (dt_perc_pcs[1]$Percent * (pcs[id_a, 1] - pcs[id_b, 1])^2 +
                       dt_perc_pcs[2]$Percent * (pcs[id_a, 2] - pcs[id_b, 2])^2)]
  
  dt_dist_stats <- dt_dist[, list(avg_dist = mean(dist)), by = list(type)]
  setkey(dt_dist_stats, type)
  signoise <- dt_dist_stats["Inter"]$avg_dist / dt_dist_stats["Intra"]$avg_dist
  signoise_db <- round(10 * log10(signoise), 3)
  
  # Plot -------------------------------------------
  p <- NULL
  if (plot) {
    # Keep plot labels/colors aligned with the matched expression samples.
    unique_samples <- unique(group)
    colors_custom <- get_sample_colors(unique_samples)
    
    text_custom_theme <- element_text(
      size = 16,
      face = "plain",
      color = "black",
      hjust = 0.5
    )
    scale_axis_x <- c(min(pcs$PC1), max(pcs$PC1))
    scale_axis_y <- c(min(pcs$PC2), max(pcs$PC2))
    
    pc1_prop <- summary(pca_prcomp)$importance[2, 1]
    pc2_prop <- summary(pca_prcomp)$importance[2, 2]
    text_axis_x <- sprintf("PC1(%.2f%%)", pc1_prop * 100)
    text_axis_y <- sprintf("PC2(%.2f%%)", pc2_prop * 100)
    limit_x <- c(1.1 * scale_axis_x[1], 1.1 * scale_axis_x[2])
    limit_y <- c(1.1 * scale_axis_y[1], 1.1 * scale_axis_y[2])
    
    
    # 修改图表标题，显示实际使用的特征数
    p_title <- paste("SNR = ", signoise_db, sep = "")
    p_subtitle <- paste("(Proteins used in PCA = ", n_features_used,
                        "/", n_features_original, ")",
                        sep = ""
    )
    # p_title <- paste("SNR = ", signoise_db, sep = "")
    # p_subtitle <- paste("(Number of proteins = ", nrow(expr_dt), ")", sep = "")
    p <- ggplot(pcs, aes(x = .data$PC1, y = .data$PC2)) +
      geom_point(aes(color = sample), size = 8) +
      theme_few() +
      theme(
        plot.title = text_custom_theme,
        plot.subtitle = text_custom_theme,
        axis.title = text_custom_theme,
        axis.text = text_custom_theme,
        legend.title = text_custom_theme,
        legend.text = element_text(size = 16, color = "gray40")
      ) +
      labs(
        x = text_axis_x,
        y = text_axis_y,
        title = p_title,
        subtitle = p_subtitle
      ) +
      scale_color_manual(values = colors_custom) +
      scale_x_continuous(limits = limit_x) +
      scale_y_continuous(limits = limit_y) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }
  
  pc_num <- ncol(pcs)
  output <- data.table(pcs[, c((pc_num - 1):pc_num, 1:(pc_num - 2))])
  
  # Save & Output ----------------------------------
  if (!is.null(output_dir)) {
    if (plot) {
      output_dir_final1 <- file.path(output_dir, "pca_plot.png")
      ggsave(output_dir_final1, p, width = 6, height = 5.5)
    }
    output_dir_final2 <- file.path(output_dir, "pca_table.tsv")
    write.table(output, output_dir_final2, sep = "\t", row.names = F)
  }
  
  # return(list(table = output, SNR = signoise_db)
  return(list(table = output, SNR = signoise_db, snr_plot = p))
}

# PCC helpers migrated from the second supplied text.
# Reference value is FC and is converted to log2FC.
pcc_normalize_colname <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(x)))

pcc_find_column <- function(dt, candidates, label = paste(candidates, collapse = "/"), required = TRUE) {
  original_names <- names(dt)
  hit <- which(pcc_normalize_colname(original_names) %in% pcc_normalize_colname(candidates))
  if (!length(hit)) {
    if (!required) return(NULL)
    stop("Missing ", label, " column. Expected one of: ", paste(candidates, collapse = ", "), ". Current columns: ", paste(original_names, collapse = ", "))
  }
  if (length(hit) > 1L) stop("Multiple columns match ", label, ": ", paste(original_names[hit], collapse = ", "))
  original_names[hit]
}

# Automatically infer raw/log2 protein or peptide abundance without interactive input.
# Thresholds are heuristics: ambiguous values default to log2.
pcc_detect_data_scale <- function(x, level_label = "Expression") {
  values <- x[is.finite(x)]
  if (!length(values)) {
    message(level_label, ": no finite values; default to log2.")
    return("log2")
  }
  if (any(values < 0)) {
    message(level_label, ": negative values detected; use log2.")
    return("log2")
  }
  q99 <- as.numeric(stats::quantile(values, probs = 0.99, names = FALSE))
  if (q99 <= 50) {
    data_scale <- "log2"
    rule <- "99th percentile <= 50"
  } else if (q99 >= 1000) {
    data_scale <- "raw"
    rule <- "99th percentile >= 1000"
  } else {
    data_scale <- "log2"
    rule <- "ambiguous scale (50 < 99th percentile < 1000); default to log2"
  }
  message(level_label, ": use ", data_scale, "; ", rule,
          " (99th percentile = ", signif(q99, 6), ").")
  data_scale
}

pcc_prepare_abundance <- function(
    native,
    scale,
    level_label,
    zero_is_missing = TRUE) {
  if (!scale %in% c("raw", "log2")) stop(level_label, " scale must be 'raw' or 'log2'.")
  native_clean <- native
  native_clean[!is.finite(native_clean)] <- NA_real_
  n_zero_reclassified <- 0L
  if (scale == "raw") {
    n_nonpositive <- sum(native_clean <= 0, na.rm = TRUE)
    if (n_nonpositive > 0L) {
      message(
        level_label, ": converted ", n_nonpositive,
        " nonpositive linear-intensity values to NA."
      )
    }
    native_clean[native_clean <= 0] <- NA_real_
    raw_mtx <- native_clean
    log2_mtx <- log2(raw_mtx)
  } else {
    if (isTRUE(zero_is_missing)) {
      n_zero_reclassified <- sum(native_clean == 0, na.rm = TRUE)
      if (n_zero_reclassified > 0L) {
        message(
          level_label, ": converted ", n_zero_reclassified,
          " exact log2 zeros to NA because zero_is_missing = TRUE."
        )
      }
      native_clean[native_clean == 0] <- NA_real_
    }
    log2_mtx <- native_clean
    raw_mtx <- 2^log2_mtx
    raw_mtx[!is.finite(raw_mtx)] <- NA_real_
  }
  list(native = native, raw = raw_mtx, log2 = log2_mtx,
       detected = is.finite(log2_mtx), missing = !is.finite(log2_mtx),
       scale = scale, zero_is_missing = zero_is_missing,
       N.Zero.Reclassified.Missing = n_zero_reclassified)
}

standardize_quant_pcc_reference <- function(reference_dt, diagnostic_dir) {
  ref <- data.table::as.data.table(reference_dt)
  sequence_name <- pcc_find_column(
    ref, c("Sequence", "peptide_sequence", "peptide"),
    "PCC reference peptide sequence"
  )
  pair_name <- pcc_find_column(
    ref, c("Sample.Pair", "Sample.pair", "sample_pair", "group", "dataset"),
    "PCC reference sample pair"
  )
  value_name <- pcc_find_column(ref, "value", "PCC quantitative reference FC value")
  protein_name <- pcc_find_column(
    ref, c("protein_id", "Protein.ID", "protein_accession"),
    "PCC reference protein ID", required = FALSE
  )
  gene_name <- pcc_find_column(
    ref, c("gene_symbol", "Gene.Symbol", "gene"),
    "PCC reference gene symbol", required = FALSE
  )
  
  # Retain annotations until record selection; do not discard them during deduplication.
  out <- data.table::data.table(
    Sequence = trimws(as.character(ref[[sequence_name]])),
    Sample.Pair = toupper(gsub("\\s+", "", trimws(as.character(ref[[pair_name]])))),
    FC.Reference = suppressWarnings(as.numeric(ref[[value_name]])),
    Protein.ID = if (is.null(protein_name)) rep(NA_character_, nrow(ref)) else trimws(as.character(ref[[protein_name]])),
    Gene.Symbol = if (is.null(gene_name)) rep(NA_character_, nrow(ref)) else trimws(as.character(ref[[gene_name]]))
  )
  out <- out[
    !is.na(Sequence) & nzchar(Sequence) &
      !is.na(Sample.Pair) & nzchar(Sample.Pair) & is.finite(FC.Reference)
  ]
  if (any(out$FC.Reference <= 0)) {
    stop("reference_dataset_quant$value must contain positive FC values.")
  }
  
  value_panel <- unique(out[, .(Sequence, Sample.Pair, FC.Reference)])
  if (nrow(value_panel) < nrow(out)) {
    message(
      "Collapsed ", nrow(out) - nrow(value_panel),
      " duplicate reference rows with identical peptide, pair and value."
    )
  }
  conflict_keys <- value_panel[, .N, by = .(Sequence, Sample.Pair)][N > 1L]
  
  if (nrow(conflict_keys)) {
    decisions <- vector("list", nrow(conflict_keys))
    for (i in seq_len(nrow(conflict_keys))) {
      peptide_i <- conflict_keys$Sequence[i]
      pair_i <- conflict_keys$Sample.Pair[i]
      candidates <- out[Sequence == peptide_i & Sample.Pair == pair_i]
      protein_ids <- candidates$Protein.ID
      genes <- candidates$Gene.Symbol
      base_ids <- sub("-[0-9]+$", "", protein_ids)
      selected_id <- NA_character_
      selected_fc <- NA_real_
      reason <- "selected unique unsuffixed protein record"
      
      if (anyNA(protein_ids) || any(!nzchar(protein_ids)) ||
          anyNA(genes) || any(!nzchar(genes))) {
        reason <- "missing protein ID or gene symbol"
      } else if (data.table::uniqueN(genes) != 1L) {
        reason <- "different gene symbols"
      } else if (any(!grepl("^[A-Za-z0-9]+(-[0-9]+)?$", protein_ids)) ||
                 data.table::uniqueN(base_ids) != 1L) {
        reason <- "different or ambiguous protein accession bases"
      } else {
        # This priority depends only on reference annotations, never on test results.
        main_rows <- candidates[Protein.ID == base_ids[1L]]
        main_values <- unique(main_rows$FC.Reference)
        if (!nrow(main_rows)) {
          reason <- "no unsuffixed protein record"
        } else if (length(main_values) != 1L) {
          reason <- "multiple FC values for the unsuffixed protein record"
        } else {
          selected_id <- base_ids[1L]
          selected_fc <- main_values[1L]
        }
      }
      decisions[[i]] <- data.table::data.table(
        Sequence = peptide_i,
        Sample.Pair = pair_i,
        Resolved = !is.na(selected_fc),
        Selected.Protein.ID = selected_id,
        Selected.FC.Reference = selected_fc,
        Reason = reason
      )
    }
    decisions <- data.table::rbindlist(decisions)
    resolved <- decisions[Resolved == TRUE]
    unresolved <- decisions[Resolved == FALSE]
    
    # Each peptide-pair contributes one supplied reference value; no averaging.
    value_panel <- data.table::rbindlist(list(
      value_panel[!conflict_keys, on = .(Sequence, Sample.Pair)],
      resolved[, .(Sequence, Sample.Pair, FC.Reference = Selected.FC.Reference)]
    ))
    
    resolved_file <- conflict_file <- NULL
    if (!is.null(diagnostic_dir)) {
      if (!dir.exists(diagnostic_dir)) {
        dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
      }
      audit <- merge(
        out, decisions, by = c("Sequence", "Sample.Pair"),
        all = FALSE, sort = TRUE
      )
      audit[, Selection.Status := "excluded: unresolved conflict"]
      audit[Resolved == TRUE, Selection.Status := "not selected: alternative annotation"]
      audit[Resolved == TRUE & Protein.ID == Selected.Protein.ID &
              FC.Reference == Selected.FC.Reference,
            Selection.Status := "selected"]
      resolved_file <- file.path(diagnostic_dir, "pcc_reference_resolved_values.tsv")
      conflict_file <- file.path(diagnostic_dir, "pcc_reference_conflicting_values.tsv")
      data.table::fwrite(audit[Resolved == TRUE], resolved_file, sep = "\t", na = "NA")
      # Write an empty table when none remain, so a prior diagnostic is not reused.
      data.table::fwrite(audit[Resolved == FALSE], conflict_file, sep = "\t", na = "NA")
    }
    if (nrow(resolved)) {
      message(
        "Resolved ", nrow(resolved),
        " peptide-pair keys by selecting the unique unsuffixed protein record. ",
        "Supplied FC values were retained; no averaging was performed.",
        if (!is.null(resolved_file)) paste0(" Details: ", resolved_file) else ""
      )
    }
    if (nrow(unresolved)) {
      warning(
        "Excluded ", nrow(unresolved),
        " peptide-pair keys without an unambiguous representative reference record. ",
        "No averaging was performed.",
        if (!is.null(conflict_file)) paste0(" Details: ", conflict_file) else ""
      )
    }
  }
  if (!nrow(value_panel)) {
    stop("No unambiguous quantitative reference peptide-pairs remain.")
  }
  
  value_panel[, log2FC.Reference := log2(FC.Reference)]
  data.table::setorder(value_panel, Sample.Pair, Sequence)
  value_panel
}

dep_analysis_log2 <- function(
    expr_log2,
    group,
    min_replicates = 2L,
    trend = TRUE,
    robust = TRUE) {
  expr_log2 <- as.matrix(expr_log2)
  storage.mode(expr_log2) <- "double"
  expr_log2[!is.finite(expr_log2)] <- NA_real_
  
  group <- droplevels(factor(group, ordered = FALSE))
  if (ncol(expr_log2) != length(group)) {
    stop("Expression columns do not match the group vector.")
  }
  if (nlevels(group) != 2L) {
    stop("Exactly two sample groups are required by dep_analysis_log2().")
  }
  if (any(table(group) < min_replicates)) return(NULL)
  
  ## The legacy combined-missing rule is applied by qc_cor_differential().
  ## Retain the minimal estimability safeguard of at least one observation in
  ## each group; this does not impose two observations per group.
  keep <- vapply(seq_len(nrow(expr_log2)), function(i) {
    z <- expr_log2[i, ]
    all(vapply(levels(group), function(g) {
      sum(is.finite(z[group == g])) >= 1L
    }, logical(1)))
  }, logical(1))
  expr_log2 <- expr_log2[keep, , drop = FALSE]
  if (!nrow(expr_log2)) return(NULL)
  
  ## The first group level is the denominator; the second is the numerator.
  design <- stats::model.matrix(
    ~ group,
    contrasts.arg = list(
      group = stats::contr.treatment(levels(group), base = 1)
    )
  )
  fit <- limma::lmFit(expr_log2, design)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  result <- limma::topTable(
    fit,
    coef = ncol(design),
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  result$Sequence <- rownames(result)
  result$Sequence.Number <- nrow(result)
  result$Sample1 <- levels(group)[1]
  result$Sample2 <- levels(group)[2]
  result$Sample.Pair <- paste(levels(group)[2], levels(group)[1], sep = "/")
  result
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  stats::cor(x[ok], y[ok], method = method)
}

qc_cor_differential <- function(
    pep_dt,
    meta_dt,
    abundance,
    reference_pcc,
    output_dir,
    plot = TRUE,
    show_sample_pairs = TRUE,
    fdr_cutoff = 0.05,
    min_replicates = 2L,
    trend = TRUE,
    robust = TRUE) {
  sample_colors <- c(D5 = "#4CC3D9", D6 = "#7BC8A4", F7 = "#FFC65D", M8 = "#F16745")
  ref_dt <- standardize_quant_pcc_reference(reference_pcc, output_dir)
  
  ## Direct limma input: normalized log2 MS abundance. Missing values remain NA.
  expr_matrix <- abundance$log2
  rownames(expr_matrix) <- trimws(as.character(pep_dt[[1]]))
  colnames(expr_matrix) <- meta_dt$library
  expr_matrix[!is.finite(expr_matrix)] <- NA_real_
  
  samples <- unique(as.character(meta_dt$sample))
  if (length(samples) < 2L) {
    stop("At least two sample groups are required for PCC calculation.")
  }
  if (!"D6" %in% samples) {
    stop("D6 is required as the denominator for Quartet PCC/RC.")
  }
  reference_sample <- "D6"
  numerators <- setdiff(samples, reference_sample)
  
  result_list <- vector("list", length(numerators))
  tested_list <- vector("list", length(numerators))
  flow_list <- vector("list", length(numerators))
  result_index <- 0L
  
  for (i in seq_along(numerators)) {
    numerator <- numerators[i]
    sample_pair <- paste(numerator, reference_sample, sep = "/")
    ref_tmp <- ref_dt[Sample.Pair == sample_pair]
    col_num <- which(meta_dt$sample == numerator)
    col_ref <- which(meta_dt$sample == reference_sample)
    
    flow_row <- data.table::data.table(
      Sample.Pair = sample_pair,
      N.Numerator.Replicates = length(col_num),
      N.Reference.Replicates = length(col_ref),
      N.Reference.Peptides = nrow(ref_tmp),
      N.Reference.In.Expression = sum(rownames(expr_matrix) %in% ref_tmp$Sequence),
      N.Passing.Missing.Rule = 0L,
      N.Reference.Passing.Missing.Rule = 0L,
      N.Limma.Tested = 0L,
      N.FDR.Significant = 0L,
      N.FDR.Significant.Matched = 0L,
      Status = "not evaluated"
    )
    
    if (!nrow(ref_tmp)) {
      flow_row[, Status := "no reference rows for this sample pair"]
      flow_list[[i]] <- flow_row
      next
    }
    if (length(col_num) < min_replicates || length(col_ref) < min_replicates) {
      flow_row[, Status := "insufficient group replicates"]
      flow_list[[i]] <- flow_row
      next
    }
    
    e_tmp <- expr_matrix[, c(col_num, col_ref), drop = FALSE]
    ## Legacy rule from the original qc_cor(): count missing observations
    ## across both groups together. For n1 vs n2, retain a peptide when
    ## N.missing < min(n1, n2). This does not require two valid values in each
    ## group (for 3 vs 3, a 1-vs-3 observed pattern can pass).
    smaller_group_size <- min(length(col_num), length(col_ref))
    keep_missing <- rowSums(!is.finite(e_tmp)) < smaller_group_size
    e_tmp <- e_tmp[keep_missing, , drop = FALSE]
    expr_grouped <- e_tmp[
      rownames(e_tmp) %in% ref_tmp$Sequence,
      ,
      drop = FALSE
    ]
    flow_row[, `:=`(
      N.Passing.Missing.Rule = nrow(e_tmp),
      N.Reference.Passing.Missing.Rule = nrow(expr_grouped)
    )]
    if (!nrow(e_tmp)) {
      flow_row[, Status := "no peptides pass the missing-value rule"]
      flow_list[[i]] <- flow_row
      next
    }
    if (!nrow(expr_grouped)) {
      flow_row[, Status := "no reference-matched peptides pass the missing-value rule"]
      flow_list[[i]] <- flow_row
      next
    }
    
    sample_pairs <- factor(
      rep(
        c(numerator, reference_sample),
        times = c(length(col_num), length(col_ref))
      ),
      levels = c(reference_sample, numerator),
      ordered = FALSE
    )
    
    result_all <- dep_analysis_log2(
      ## Only reference-matched peptides are tested. BH adjustment is
      ## performed within the quantitative reference panel for this pair.
      expr_log2 = expr_grouped,
      group = sample_pairs,
      min_replicates = min_replicates,
      trend = trend,
      robust = robust
    )
    if (is.null(result_all) || !nrow(result_all)) {
      flow_row[, Status := "limma returned no tested peptides"]
      flow_list[[i]] <- flow_row
      next
    }
    result_all <- data.table::as.data.table(result_all)
    tested_list[[i]] <- result_all
    flow_row[, N.Limma.Tested := nrow(result_all)]
    
    result_sig <- result_all[
      is.finite(adj.P.Val) & adj.P.Val < fdr_cutoff
    ]
    flow_row[, N.FDR.Significant := nrow(result_sig)]
    if (!nrow(result_sig)) {
      flow_row[, Status := "no peptides pass the FDR cutoff"]
      flow_list[[i]] <- flow_row
      next
    }
    
    result_index <- result_index + 1L
    result_list[[result_index]] <- result_sig
    flow_row[, Status := "FDR-significant peptides available"]
    flow_list[[i]] <- flow_row
  }
  
  tested_all <- data.table::rbindlist(tested_list, fill = TRUE)
  flow_summary <- data.table::rbindlist(flow_list, fill = TRUE)
  result_list <- result_list[seq_len(result_index)]
  
  if (!length(result_list)) {
    return(list(
      DEPs = NULL,
      All.Tests = tested_all,
      logfc = NULL,
      COR = NA_real_,
      cor_plot = NULL,
      summary = data.table::data.table(),
      flow = flow_summary,
      FDR.Cutoff = fdr_cutoff,
      Min.Replicates = min_replicates,
      Method = "reference-panel limma on normalized log2 abundance"
    ))
  }
  
  result_final <- data.table::rbindlist(result_list, fill = TRUE)
  result_trim <- result_final[, .(
    Sequence,
    Sample.Pair,
    logFC.Test = as.numeric(logFC),
    FC.Test = 2^as.numeric(logFC),
    AveExpr = as.numeric(AveExpr),
    P.Value = as.numeric(P.Value),
    adj.P.Val = as.numeric(adj.P.Val),
    t = as.numeric(t),
    B = as.numeric(B)
  )]
  result_withref <- merge(
    result_trim,
    ref_dt,
    by = c("Sequence", "Sample.Pair"),
    all = FALSE,
    sort = FALSE
  )
  result_withref <- result_withref[
    is.finite(logFC.Test) & is.finite(log2FC.Reference)
  ]
  result_withref[, Difference := logFC.Test - log2FC.Reference]
  
  if (!nrow(result_withref)) {
    return(list(
      DEPs = result_final,
      All.Tests = tested_all,
      logfc = result_withref,
      COR = NA_real_,
      cor_plot = NULL,
      summary = data.table::data.table(),
      flow = flow_summary,
      FDR.Cutoff = fdr_cutoff,
      Min.Replicates = min_replicates,
      Method = paste0(
        "reference-panel limma on normalized log2 abundance; BH within matched reference peptides; eBayes trend=",
        trend,
        ", robust=",
        robust
      )
    ))
  }
  
  if (nrow(flow_summary)) {
    matched_counts <- result_withref[, .(
      N.FDR.Significant.Matched = .N
    ), by = Sample.Pair]
    flow_summary[matched_counts,
                 N.FDR.Significant.Matched := i.N.FDR.Significant.Matched,
                 on = "Sample.Pair"]
    flow_summary[N.FDR.Significant.Matched > 0L,
                 Status := "PCC-ready significant matched peptides"]
  }
  
  cor_value <- safe_cor(
    result_withref$logFC.Test,
    result_withref$log2FC.Reference,
    method = "pearson"
  )
  cor_value <- if (is.finite(cor_value)) round(cor_value, 3) else NA_real_
  
  pair_summary <- result_withref[, {
    fit_pair <- if (.N >= 2L) {
      stats::lm(logFC.Test ~ log2FC.Reference)
    } else {
      NULL
    }
    .(
      N.Significant.Matched.Peptides = .N,
      PCC = safe_cor(logFC.Test, log2FC.Reference, "pearson"),
      Spearman = safe_cor(logFC.Test, log2FC.Reference, "spearman"),
      Direction.Concordance = mean(
        sign(logFC.Test) == sign(log2FC.Reference),
        na.rm = TRUE
      ),
      Reference.SD = sd(log2FC.Reference),
      Test.SD = sd(logFC.Test),
      MAE = mean(abs(Difference)),
      RMSE = sqrt(mean(Difference^2)),
      Intercept = if (is.null(fit_pair)) NA_real_ else unname(coef(fit_pair)[1]),
      Slope = if (is.null(fit_pair)) NA_real_ else unname(coef(fit_pair)[2])
    )
  }, by = Sample.Pair]
  
  pooled_summary <- result_withref[, .(
    Sample.Pair = "OVERALL_POOLED",
    N.Significant.Matched.Peptides = .N,
    PCC = safe_cor(logFC.Test, log2FC.Reference, "pearson"),
    Spearman = safe_cor(logFC.Test, log2FC.Reference, "spearman"),
    Direction.Concordance = mean(
      sign(logFC.Test) == sign(log2FC.Reference),
      na.rm = TRUE
    ),
    Reference.SD = sd(log2FC.Reference),
    Test.SD = sd(logFC.Test),
    MAE = mean(abs(Difference)),
    RMSE = sqrt(mean(Difference^2)),
    Intercept = unname(coef(lm(logFC.Test ~ log2FC.Reference))[1]),
    Slope = unname(coef(lm(logFC.Test ~ log2FC.Reference))[2])
  )]
  pcc_summary <- data.table::rbindlist(list(pair_summary, pooled_summary), fill = TRUE)
  
  p <- NULL
  if (plot && nrow(result_withref)) {
    limit <- max(abs(c(
      result_withref$log2FC.Reference,
      result_withref$logFC.Test
    )), na.rm = TRUE)
    if (!is.finite(limit) || limit == 0) limit <- 1
    limit_axis <- c(-limit * 1.05, limit * 1.05)
    
    unique_comps <- unique(result_withref$Sample.Pair)
    pair_colors <- vapply(unique_comps, function(comp_name) {
      numerator_sample <- strsplit(comp_name, "/", fixed = TRUE)[[1]][1]
      if (numerator_sample %in% names(sample_colors)) {
        unname(sample_colors[numerator_sample])
      } else {
        "gray"
      }
    }, character(1))
    names(pair_colors) <- unique_comps
    
    p <- ggplot2::ggplot(
      result_withref,
      ggplot2::aes(x = .data$log2FC.Reference, y = .data$logFC.Test)
    ) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray55") +
      ggthemes::theme_few() +
      ggplot2::labs(
        y = "log2FC (Test; limma logFC)",
        x = "log2FC (Reference; log2(value))",
        title = paste0(
          "PCC = ",
          ifelse(is.na(cor_value), "NA", sprintf("%.3f", cor_value))
        ),
        subtitle = paste0(
          "Significant matched peptide-pairs (FDR < ",
          fdr_cutoff,
          ") = ",
          nrow(result_withref)
        )
      ) +
      ggplot2::coord_fixed(xlim = limit_axis, ylim = limit_axis) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, color = "black", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 12, color = "black", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 14, color = "black"),
        axis.text = ggplot2::element_text(size = 12, color = "black"),
        legend.title = ggplot2::element_text(size = 13, color = "black"),
        legend.text = ggplot2::element_text(size = 12, color = "gray30")
      )
    if (show_sample_pairs) {
      p <- p +
        ggplot2::geom_point(ggplot2::aes(color = .data$Sample.Pair), size = 2.5, alpha = 0.55) +
        ggplot2::scale_color_manual(values = pair_colors, name = "Sample Pair")
    } else {
      p <- p + ggplot2::geom_point(color = "steelblue4", size = 2.5, alpha = 0.4)
    }
  }
  
  list(
    DEPs = result_final,
    All.Tests = tested_all,
    logfc = result_withref,
    COR = cor_value,
    cor_plot = p,
    summary = pcc_summary,
    flow = flow_summary,
    FDR.Cutoff = fdr_cutoff,
    Min.Replicates = min_replicates,
    Method = paste0(
      "reference-panel limma on normalized log2 abundance; BH within matched reference peptides; eBayes trend=",
      trend,
      ", robust=",
      robust
    )
  )
}

#' Calculating RC value; Plotting a scatterplot
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @param show_sample_pairs if True, samples in plot will be labeled.
#' @import stats
#' @import utils
#' @importFrom rlang .data
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggthemes theme_few
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggsave
#' @export
qc_cor <- function(expr_dt, meta_dt,
                   output_dir = NULL, plot = FALSE, show_sample_pairs = TRUE) {
  # Adapt the second text's PCC implementation to the original package API.
  # Only PCC uses these defaults; the other QC functions are unchanged.
  fdr_cutoff <- 0.05
  min_replicates <- 2L
  trend <- TRUE
  robust <- TRUE
  if (robust && !requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required for limma::eBayes(robust = TRUE).")
  }
  
  # Use the second text's quantitative reference: value is FC, not log2FC.
  ref_env <- new.env(parent = emptyenv())
  utils::data("reference_dataset_quant", package = "protqc", envir = ref_env)
  if (!exists("reference_dataset_quant", envir = ref_env, inherits = FALSE)) {
    stop("Built-in reference_dataset_quant was not found in protqc.")
  }
  reference_pcc <- get("reference_dataset_quant", envir = ref_env, inherits = FALSE)
  
  # Automatically infer raw/log2 scale; supply log2 abundance to limma.
  expr_df <- as.data.frame(expr_dt, check.names = FALSE)
  expr_matrix <- as.matrix(expr_df[, -1L, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- expr_df[[1L]]
  expr_matrix <- expr_matrix[, meta_dt$library, drop = FALSE]
  expr_matrix[!is.finite(expr_matrix)] <- NA_real_
  data_scale <- pcc_detect_data_scale(expr_matrix, level_label = "Peptide PCC")
  # As in the second text, zeros encode missing observations and NA stays NA.
  abundance <- pcc_prepare_abundance(
    native = expr_matrix, scale = data_scale,
    level_label = "Peptide PCC", zero_is_missing = TRUE
  )
  
  # Needed only for the second text's conflicting-reference diagnostic.
  diagnostic_dir <- if (is.null(output_dir)) tempdir() else output_dir
  if (!dir.exists(diagnostic_dir)) {
    dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
  }
  result <- qc_cor_differential(
    pep_dt = expr_df, meta_dt = meta_dt, abundance = abundance,
    reference_pcc = reference_pcc, output_dir = diagnostic_dir,
    plot = plot, show_sample_pairs = show_sample_pairs,
    fdr_cutoff = fdr_cutoff, min_replicates = min_replicates,
    trend = trend, robust = robust
  )
  # Preserve the former column alias for callers using qc_cor()$logfc.
  if (!is.null(result$logfc) && "log2FC.Reference" %in% names(result$logfc)) {
    result$logfc[, logFC.Reference := log2FC.Reference]
  }
  result
}

# Recall helpers migrated from the newly supplied text.
standardize_quali_reference <- function(reference_dt) {
  ref <- data.table::as.data.table(reference_dt)
  sample_name <- pcc_find_column(ref, c("sample", "sample_type", "group"), "qualitative reference sample")
  sequence_name <- pcc_find_column(ref, c("peptide_sequence", "Sequence", "peptide"), "qualitative reference peptide sequence")
  out <- data.table::data.table(Sample = as.character(ref[[sample_name]]), Sequence = as.character(ref[[sequence_name]]))
  out[, Sample := toupper(gsub("^QUARTET[ _-]*", "", trimws(Sample), ignore.case = TRUE))]
  out[, Sequence := trimws(Sequence)]
  unique(out[!is.na(Sample) & nzchar(Sample) & !is.na(Sequence) & nzchar(Sequence)])
}

qc_qualitative_recall <- function(pep_dt, meta_dt, abundance, reference_quali) {
  ref <- standardize_quali_reference(reference_quali)
  feature_id <- trimws(as.character(pep_dt[[1]]))
  detected <- abundance$detected
  rownames(detected) <- feature_id
  common_groups <- intersect(unique(meta_dt$sample), unique(ref$Sample))
  if (!length(common_groups)) stop("No common samples between metadata and qualitative reference.")
  details <- data.table::rbindlist(lapply(common_groups, function(g) {
    cols <- which(meta_dt$sample == g)
    n_replicates <- length(cols)
    min_required <- floor(n_replicates / 2) + 1L
    n_detected <- rowSums(detected[, cols, drop = FALSE])
    names(n_detected) <- feature_id
    rg <- ref[Sample == g]
    counts <- unname(n_detected[match(rg$Sequence, names(n_detected))])
    counts[is.na(counts)] <- 0L
    data.table::data.table(Sample = g, Sequence = rg$Sequence, N.Replicates = n_replicates, Min.Replicates.Required = min_required,
                           N.Detected.Replicates = as.integer(counts), Detected.AnyReplicate = counts >= 1L, Detected.MajorityReplicates = counts >= min_required)
  }))
  by_group <- details[, .(N.Reference.Peptides = .N, N.Detected.AnyReplicate = sum(Detected.AnyReplicate),
                          N.Detected.MajorityReplicates = sum(Detected.MajorityReplicates), Recall.AnyReplicate = mean(Detected.AnyReplicate),
                          Recall.MajorityReplicates = mean(Detected.MajorityReplicates), N.Replicates = unique(N.Replicates),
                          Min.Replicates.Required = unique(Min.Replicates.Required)), by = Sample]
  overall_micro <- details[, .(Sample = "OVERALL_MICRO", N.Reference.Peptides = .N,
                               N.Detected.AnyReplicate = sum(Detected.AnyReplicate), N.Detected.MajorityReplicates = sum(Detected.MajorityReplicates),
                               Recall.AnyReplicate = mean(Detected.AnyReplicate), Recall.MajorityReplicates = mean(Detected.MajorityReplicates),
                               N.Replicates = NA_integer_, Min.Replicates.Required = NA_integer_)]
  overall_macro <- data.table::data.table(Sample = "OVERALL_MACRO", N.Reference.Peptides = sum(by_group$N.Reference.Peptides),
                                          N.Detected.AnyReplicate = sum(by_group$N.Detected.AnyReplicate), N.Detected.MajorityReplicates = sum(by_group$N.Detected.MajorityReplicates),
                                          Recall.AnyReplicate = mean(by_group$Recall.AnyReplicate), Recall.MajorityReplicates = mean(by_group$Recall.MajorityReplicates),
                                          N.Replicates = NA_integer_, Min.Replicates.Required = NA_integer_)
  list(details = details, by_group = by_group, summary = data.table::rbindlist(list(by_group, overall_micro, overall_macro), fill = TRUE),
       Recall.AnyReplicate = overall_micro$Recall.AnyReplicate, Recall.MajorityReplicates = overall_micro$Recall.MajorityReplicates,
       Macro.Recall.AnyReplicate = overall_macro$Recall.AnyReplicate, Macro.Recall.MajorityReplicates = overall_macro$Recall.MajorityReplicates)
}

#' Calculating Recall
#' @description Uses the supplied text's majority-replicate micro Recall.
#' A reference peptide must be detected in floor(n / 2) + 1 replicates.
#' Overall Recall is total detected reference peptide-sample pairs divided
#' by total reference peptide-sample pairs for the sample groups present.
#' Reference pairs are deduplicated; missing expression peptides count as undetected.
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file (to map library to sample type)
#' @export
qc_recall <- function(expr_dt, meta_dt) {
  # Same package entry point and numeric return type; only Recall changes.
  utils::data("reference_dataset_quali", package = "protqc", envir = environment())
  reference_quali <- reference_dataset_quali
  
  expr_df <- as.data.frame(expr_dt, check.names = FALSE)
  meta <- as.data.frame(meta_dt, stringsAsFactors = FALSE)
  meta$library <- as.character(meta$library)
  meta$sample <- toupper(gsub("^QUARTET[ _-]*", "",
                              trimws(as.character(meta$sample)), ignore.case = TRUE))
  missing_libraries <- setdiff(meta$library, names(expr_df)[-1L])
  if (length(missing_libraries)) {
    stop("Libraries absent from the peptide matrix: ",
         paste(missing_libraries, collapse = ", "))
  }
  expr_matrix <- as.matrix(expr_df[, meta$library, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  expr_matrix[!is.finite(expr_matrix)] <- NA_real_
  
  # Use the same automatic scale rule as PCC; ambiguous values default to log2.
  # Exact zeros remain missing placeholders; negative log2 values are valid.
  data_scale <- pcc_detect_data_scale(expr_matrix, level_label = "Peptide Recall")
  abundance <- pcc_prepare_abundance(
    native = expr_matrix, scale = data_scale,
    level_label = "Peptide Recall", zero_is_missing = TRUE
  )
  ref <- standardize_quali_reference(reference_quali)
  if (!length(intersect(unique(meta$sample), unique(ref$Sample)))) {
    warning("No common sample types found between metadata and reference dataset.")
    return(NA_real_)
  }
  recall_results <- qc_qualitative_recall(
    pep_dt = expr_df, meta_dt = meta,
    abundance = abundance, reference_quali = ref
  )
  # Weighted micro Recall, using majority detection as the main definition.
  # sum_g(N.detected.majority[g]) / sum_g(N.reference[g])
  # Any-replicate and macro summaries remain available in the helper result.
  as.numeric(recall_results$Recall.MajorityReplicates)
}
