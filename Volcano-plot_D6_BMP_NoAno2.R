library(ggplot2)
library(dplyr)
library(ggrepel)


# Read data
data <- read.table("D6_BMP.txt", header = TRUE, row.names = 1)

neuron_associated_genes <- read.table("Dev_genes_associated_list.txt", header = TRUE)

# Calculate the mean FPKM for DKO and WT
data$mean_fpkm_DKO <- rowMeans(data[, c("DKO_F12_B_D6", "DKO_I28_B_D6", "DKO_I35_B_D6")])
data$mean_fpkm_WT <- rowMeans(data[, c("HM1_1_B_D6", "HM1_2_B_D6", "HM1_3_B_D6")])

# Filter genes based on your criteria
data$de <- ifelse((data$mean_fpkm_DKO >= 5 | data$mean_fpkm_WT >= 5) &
                    abs(data$diffexp_log2fc_HM1_BMP_D6_vs_DKO_BMP_D6) > 0.56 &
                    data$diffexp_deseq2_qvalue_HM1_BMP_D6_vs_DKO_BMP_D6 < 0.01, "DE", "NonDE")

data$Description <- "C_Non_DE"
data$Description[data$de == "DE" & data$diffexp_log2fc_HM1_BMP_D6_vs_DKO_BMP_D6 > 0] <- "B_Upregulated in DKO"
data$Description[data$de == "DE" & data$diffexp_log2fc_HM1_BMP_D6_vs_DKO_BMP_D6 < 0] <- "A_Downregulated in DKO"



# Filtrer les données pour n'inclure que les gènes étiquetés comme "Downregulated in DKO"
df_sorted <- data %>%
  arrange(desc(Description))

# Tracer le graphique avec uniquement les gènes "Downregulated in DKO"
# Calculer le nombre de gènes DE avec un FC inférieur à 0 et supérieur à 0
group_counts <- df_sorted %>%
  filter(Description %in% c("A_Downregulated in DKO", "B_Upregulated in DKO")) %>%
  group_by(Description) %>%
  summarise(count = n())

df_sorted$Description <- factor(df_sorted$Description)

#Pour les qvalue = 0, les remplacer par une qvalue plotable
df_sorted$diffexp_deseq2_qvalue_HM1_BMP_D6_vs_DKO_BMP_D6[df_sorted$diffexp_deseq2_qvalue_HM1_BMP_D6_vs_DKO_BMP_D6 == 0] <- 4.51501634660052E-299


# Ajouter une colonne "alpha" au data frame
#df_sorted$transparence <- 0.25
#df_sorted$transparence <- ifelse(as.character(row.names(df_sorted)) %in% as.character(neuron_associated_genes$Genes) & df_sorted$Description %in% c("A_Downregulated in DKO", "B_Upregulated in DKO"), 1, 0.25)

#highlight_genes_vector <- c("Olig3", "Msx1", "Atoh1", "Zic2", "Dbx1", "Dbx2", "Prdm12", "Cdon", "FoxD3", "Sox21", "Prdm13", "Irx3", "Irx5", "Irx1")
#highlight_genes_rows <- rownames(df_sorted)[row.names(df_sorted) %in% highlight_genes_vector]
 
# Mettre à jour highlight_genes_vector avec les noms de lignes correspondants
#highlight_genes_vector <- highlight_genes_rows
 
# Créer le data frame highlight_genes
#highlight_genes <- data.frame(Genes = highlight_genes_vector)
 
# Sélectionner les lignes correspondant aux gènes d'intérêt
#highlighted_points <- df_sorted[row.names(df_sorted) %in% highlight_genes_vector & df_sorted$de == "DE", ]


# Créer le graphique
pdf(file = "Volcano_plot_Day6BMP_noAno2.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 4)

# Tracer le graphique avec les valeurs alpha mises à jour
ggplot(df_sorted, aes(x = -diffexp_log2fc_HM1_BMP_D6_vs_DKO_BMP_D6, y = -log10(diffexp_deseq2_qvalue_HM1_BMP_D6_vs_DKO_BMP_D6), color = Description)) +
  coord_cartesian(xlim = c(-11, 11), ylim = c(0, 320)) +# définir les limites du graphique pour inclure toutes les étiquettes de texte
  geom_point(size = 2) +
  labs(x = "Log2 fold change", y = "-Log10 q-value") +
  geom_vline(xintercept = c(-0.56, 0.56), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linetype = "solid"),
        axis.ticks.length = unit(0.25, "cm")) +
  scale_alpha_identity(name = "Alpha", guide = "none")  +
  scale_color_manual(values = c("#266799", "#F68532", "gray"),
                     labels = c("Downregulated in DKO", "Upregulated in DKO", "Non DE")) +
  theme(legend.position="none")


dev.off()








