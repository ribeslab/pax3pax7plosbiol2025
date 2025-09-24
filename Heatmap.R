library(pheatmap)
library(viridis)
data <- read.table("p01.txt", header = TRUE, row.names = 1)

#data_log <- log2(data + 0.001)

# Calculate Z-scores across all genes (normalizes across samples)
#data_zscore <- t(scale(t(data)))  # Transpose, scale, and transpose back

# Alternatively, use a power transformation to compress the range
#data_power <- data^(1/3)  # Cube root transformation

# Créer le graphique
pdf(file = "HEATMAP_p01.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 6)

pheatmap(data, 
         scale = "none",                # Row scaling for better contrast in gene expression
         #color = viridis(n = 256, alpha = 1, 
         #               begin = 0, end = 1, option = "viridis"), # Color gradient
         color = colorRampPalette(c("gray85", "#FDE725"))(100),
         main = "Gene Expression Heatmap",
         cellwidth = 6,
         cellheight = 15,
         fontsize_row = 6,              # Adjust row font size for readability
         fontsize_col = 8,              # Adjust column font size
         cluster_rows = FALSE,           # Clustering genes
         cluster_cols = FALSE)           # Clustering samples

dev.off()


library(pheatmap)

# Log-transform data (adding 1 to avoid log(0))
data_log <- log2(data + 1)

# Define a custom color scale
# Choose a color gradient where low values are closer to white and high values are a gradual red
color_scale <- colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100)

# Generate the heatmap
pheatmap(data_log,
         color = color_scale,
         breaks = seq(min(data_log), max(data_log), length.out = 100), # Custom breaks for color scaling
         main = "Gene Expression Heatmap (Log-Scaled)",
         fontsize_row = 6,
         fontsize_col = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         legend_breaks = c(0, 2, 4, 6, 8, 10), # Adjust based on the log scale
         legend_labels = c("0", "2", "4", "6", "8", "10"))
