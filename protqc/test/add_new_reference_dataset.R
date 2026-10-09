# 步骤 1: 读取 CSV 文件
# 确保你的 R 工作目录中有 7_quantprop_final.csv 文件
raw_data <- read.csv("./test/7_quantprop_final.csv", stringsAsFactors = FALSE)

# (可选) 如果需要，可以在这里对 raw_data 进行预处理
# 例如：重命名列，转换数据类型等
reference_dataset_quant <- raw_data 

# 步骤 2: 将其保存到你的 R 包中
# 方法 A: 使用 devtools/usethis (推荐，会自动处理路径和压缩)
# 确保你当前的 R 工作目录是 protqc 包的根目录
usethis::use_data(reference_dataset_quant, overwrite = TRUE)


# 添加定性数据集

raw_data <- read.csv("./test/7_qualiprop_final.csv", stringsAsFactors = FALSE)

reference_dataset_quali <- raw_data 

usethis::use_data(reference_dataset_quali, overwrite = TRUE)
