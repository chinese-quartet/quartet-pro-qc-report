# ---------------------------------------------------------------------------- #
#' @title Generate Quartet Proteomics report
#'
#' @description Use calculated Met result to generate report
#'
#' @param qc_result list
#' @param report_template character
#' @param report_dir character
#' @param report_name character
#'
#' @importFrom dplyr %>%
#' @importFrom flextable flextable
#' @importFrom flextable theme_vanilla
#' @importFrom flextable color
#' @importFrom flextable set_caption
#' @importFrom flextable align
#' @importFrom flextable width
#' @importFrom flextable bold
#' @importFrom flextable theme_box
#' @importFrom flextable bg
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom officer body_add_par
#' @importFrom flextable body_add_flextable
#' @importFrom officer body_add_gg
#' @importFrom officer body_add_break
#' @importFrom officer read_docx
#'
#' @examples
#' # 加载示例数据
#' example_prot_data_path <- system.file("extdata", "proteomics_pipeline_data_template.csv", package = "protqc")
#' example_prot_metadata_path <- system.file("extdata", "proteomics_pipeline_meta_template.csv", package = "protqc")
#' qc_result <- qc_conclusion(example_prot_data_path, example_prot_metadata_path)
#'
#' # 指定包内文档的路径
#' report_template <- system.file("extdata", "quartet_template.docx", package = "protqc")
#'
#' # 运行函数
#' generate_protein_report(qc_result = qc_result, report_template = report_template)
#'
#' @export
generate_protein_report <- function(qc_result,
                                    report_template,
                                    report_dir = NULL,
                                    report_name = NULL) {
  if (is.null(qc_result) || is.null(report_template)) {
    stop("All arguments (qc_result, report_template) are required.")
  }
  
  if (is.null(report_dir)) {
    path <- getwd()
    report_dir <- file.path(path, "output")
    dir.create(report_dir, showWarnings = FALSE)
  }
  
  ### 读取quarter报告模板并生成报告
  if (is.null(report_name)) {
    report_name <- "Quartet_protein_report.docx"
  }
  output_file <- file.path(report_dir, report_name)
  
  
  # 有绘图对象时插入图片；空图改为说明文字，继续生成报告。
  add_report_plot <- function(doc, plot_obj, empty_text) {
    if (is.null(plot_obj)) {
      return(body_add_par(doc, value = empty_text, style = "Normal"))
    }
    body_add_gg(doc, value = plot_obj, style = "centered")
  }
  
  # --- 1. 定义中文文本内容 (基于 ShenKang-Quartet-Protein-Report_v0.2.docx) ---
  
  # 摘要
  text_abstract <- "本报告基于多项组学关键质量控制指标，总结了 Quartet Protein参考物质所生成数据的质量情况。质量控制流程从用户输入肽段表达矩阵开始，分别计算外部质控品的检测特征数（Number of features）、标称特性灵敏度 (Recall of nominal characteristics)、信噪比 (Signal-to-Noise Ratio, SNR)、与参考数据集的Pearson相关系数（Pearson correlation coefficient, PCC）及整体质量判断。"
  
  # 指标定义
  def_features_title <- "检测特征数（Number of features）"
  def_features_desc <- "该指标反映在蛋白组检测中成功鉴定并映射到基因符号的蛋白数量。通常期望获得尽可能多的蛋白特征，以支持后续的生物学分析。"
  
  def_recall_title <- "标称特性灵敏度（Recall of nominal characteristics）"
  def_recall_desc <- "定义为标称特性肽段在测试数据集中被成功检测到的比例。"
  
  def_snr_title <- "信噪比（Signal-to-Noise Ratio, SNR）"
  def_snr_desc <- "用于刻画某一检测平台、实验室或批次区分不同生物样本组之间内在生物学差异（“信号”）与同一样本组技术重复变异（“噪声”）的能力。SNR 越高，表明区分多组样本差异的能力越强。"
  
  def_rc_title <- "与参考数据集的Pearson相关系数（Pearson correlation coefficient, PCC）"
  def_rc_desc <- "用于评估测试数据在相对定量层面与参考数据集之间的一致性。在 shotgun 蛋白组学中，肽段层面的定量在理论上更为可靠。因此，参考数据集基于历史数据中各样本对 (D5/D6、F7/D6) 在肽段水平的相对表达值（log2FC）构建。在评估过程中，仅选取在测试数据与参考数据集中均覆盖、且满足统计学阈值（BH 校正后的 FDR < 0.05）的肽段 log2FC 作为输入，计算测试数据与参考数据之间的 Pearson 相关系数，作为 PCC 指标。"
  
  # 参考文献
  text_refs <- c(
    "1. Zheng Y, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology, 2024",
    "2. Tian S, et al. Quartet protein reference materials and datasets for multi-platform assessment of label-free proteomics. Genome Biology, 2023.",
    "3. 上海临床队列组学检测工作指引（征求意见稿）, 2025/11/26"
  )
  
  # 免责声明
  text_disclaimer <- "本数据质量报告仅针对所评估的特定数据集提供分析结果，仅供信息参考之用。尽管已尽最大努力确保分析结果的准确性和可靠性，但本报告按“现状（AS IS）”提供，不附带任何形式的明示或暗示担保。报告作者及发布方不对基于本报告内容所采取的任何行动承担责任。本报告中的结论不应被视为对任何产品或流程质量的最终判定，也不应用于关键应用场景、商业决策或法规合规用途，除非经过专业核查和独立验证。对于分析结果的正确性、准确性、可靠性或适用性，不作任何明示或暗示的保证。"
  
  # --- 2. 构建指标数据框 ---
  
  # 假设 qc_result$conclusion 是结果表格
  # 假设 qc_result$results 包含 snr_plot 和 cor_plot
  raw_table <- qc_result$conclusion
  
  # 辅助函数：安全提取 Value
  get_val <- function(df, pattern) {
    # 找到包含 Metric 名称的行
    # 兼容 "Quality Metrics" 或第一列
    col_idx <- grep("Quality.*Metrics", names(df), ignore.case = TRUE)
    if (length(col_idx) == 0) col_idx <- 1
    
    # 找到 Value 对应的列
    val_idx <- grep("^Value$", names(df), ignore.case = TRUE)
    if (length(val_idx) == 0) val_idx <- 2
    
    # 提取
    row_idx <- grep(pattern, df[[col_idx]], ignore.case = TRUE)
    if (length(row_idx) > 0) {
      return(as.numeric(df[[val_idx]][row_idx[1]]))
    }
    return(NA)
  }
  
  val_feat <- get_val(raw_table, "Number of features")
  val_recall <- get_val(raw_table, "Recall")
  val_snr <- get_val(raw_table, "Signal-to-Noise Ratio")
  val_rc <- get_val(raw_table, "^PCC$")
  
  batch_name <- "Query_data" # 可改为从参数传入
  
  # 格式化与阈值判断
  # 1. Features
  str_feat <- ifelse(is.na(val_feat), "-", sprintf("%.0f", val_feat))
  
  # 2. Recall (>= 0.90)
  if (is.na(val_recall)) {
    str_recall <- "-"
  } else {
    str_recall <- sprintf("%.2f", val_recall)
    # if (val_recall < 0.90) str_recall <- paste0(str_recall, " ↓")
  }
  
  # 3. SNR (>= 10)
  if (is.na(val_snr)) {
    str_snr <- "-"
  } else {
    str_snr <- sprintf("%.2f", val_snr)
    if (val_snr < 10) str_snr <- paste0(str_snr, " ↓")
  }
  
  # 4. PCC (>= 0.80)
  if (is.na(val_rc)) {
    str_rc <- "-"
  } else {
    str_rc <- sprintf("%.2f", val_rc)
    # if (val_rc < 0.80) str_rc <- paste0(str_rc, " ↓")
  }
  
  # 整体判断：Recall、SNR 和 PCC 均达到推荐标准才通过
  if (!is.na(val_recall) && !is.na(val_snr) && !is.na(val_rc)) {
    pass <- (val_recall >= 0.90) &&
      (val_snr >= 10) &&
      (val_rc >= 0.80)
    str_qual <- ifelse(pass, "Yes", "No")
  } else {
    str_qual <- "No" # 任一指标缺失时不通过
  }
  
  # 构建 2行 x 6列 数据框
  df_table <- data.frame(
    "样本组" = c("推荐质量标准", batch_name),
    "检测特征数" = c("-", str_feat),
    "标称特性灵敏度" = c("≥0.90", str_recall),
    "信噪比" = c("≥10", str_snr),
    "Pearson相关系数" = c("≥0.80", str_rc),
    "是否通过" = c("-", str_qual),
    check.names = FALSE
  )
  
  # 创建 Flextable
  ft1 <- flextable(df_table) %>%
    theme_box() %>%
    flextable::font(part = "all", fontname = "Times New Roman") %>%
    align(align = "center", part = "all") %>%
    width(width = 1.1) %>% # 适当调整列宽
    bold(part = "header") %>%
    bold(i = 1, part = "body") %>%
    bg(part = "header", bg = "#EFEFEF") %>%
    color(i = 2, j = "是否通过", color = ifelse(str_qual == "No", "#B80D0D", "black"))
  
  # --- 3. 生成报告流程 (Pipeline) ---
  
  doc <- read_docx(report_template) %>%
    # 1. 标题
    body_add_par(value = "Quartet蛋白肽段组质量报告", style = "heading 1") %>%
    # 2. 摘要
    body_add_par(value = "摘要", style = "heading 2") %>%
    body_add_par(value = text_abstract, style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>% # 空行
    
    # 3. 表格
    body_add_flextable(ft1) %>%
    # 4. 质量控制指标定义
    body_add_par(value = "质量控制指标", style = "heading 2") %>%
    body_add_par(value = def_features_title, style = "heading 3") %>%
    body_add_par(value = def_features_desc, style = "Normal") %>%
    body_add_par(value = def_recall_title, style = "heading 3") %>%
    body_add_par(value = def_recall_desc, style = "Normal") %>%
    body_add_par(value = def_snr_title, style = "heading 3") %>%
    body_add_par(value = def_snr_desc, style = "Normal") %>%
    body_add_par(value = def_rc_title, style = "heading 3") %>%
    body_add_par(value = def_rc_desc, style = "Normal") %>%
    # # 5. 图片 (SNR 和 Correlation)
    # # 确保 qc_result 中有这些 plot 对象
    # body_add_par(value = "Signal-to-Noise Ratio", style = "heading 2") %>%
    # body_add_gg(value = qc_result$results$snr_results$snr_plot, style = "centered") %>%
    # body_add_break() %>%
    #
    # body_add_par(value = "Correlation with Reference Datasets", style = "heading 2") %>%
    # body_add_gg(value = qc_result$results$cor_results$cor_plot, style = "centered") %>%
    # body_add_break() %>%
    
    # 5. 参考文献
    body_add_par(value = "参考文献", style = "heading 2") %>%
    body_add_par(value = text_refs[1], style = "Normal") %>%
    body_add_par(value = text_refs[2], style = "Normal") %>%
    body_add_par(value = text_refs[3], style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>%
    # 6. 免责声明
    body_add_par(value = "免责声明", style = "heading 3") %>%
    body_add_par(value = text_disclaimer, style = "Normal") %>%
    body_add_break() %>%
    # 7. 图片 (SNR 和 PCC)
    body_add_par(value = "Signal-to-Noise Ratio", style = "heading 2") %>%
    add_report_plot(
      plot_obj = qc_result$results$snr_results$snr_plot,
      empty_text = "未生成 SNR 图：未开启绘图或没有可用的绘图对象。"
    ) %>%
    body_add_break() %>%
    body_add_par(value = "Pearson Correlation Coefficient", style = "heading 2") %>%
    add_report_plot(
      plot_obj = qc_result$results$cor_results$cor_plot,
      empty_text = if (is.na(val_rc)) {
        "当前 PCC 无法计算，因此未生成相关性图；可能没有满足筛选条件的匹配肽段，或有效数据不足。"
      } else {
        "未生成 PCC 图：未开启绘图或没有可用的绘图对象。"
      }
    )
  
  # 输出文件
  print(doc, target = output_file)
}
