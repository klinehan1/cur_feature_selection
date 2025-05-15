library(readxl)

df <- read_excel("GSE10072_series_matrix.xlsx", sheet = 1)

# find patients with tumor
idx <- grepl("Tumor", df[47,2:108], fixed = TRUE)
tumor_idx <- which(idx == TRUE)

write.csv(tumor_idx, file = "tumor_idx.csv", row.names = FALSE)
