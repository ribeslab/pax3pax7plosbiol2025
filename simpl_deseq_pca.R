

library(ggplot2)
library(DESeq2)

raw_counts <- read.table("rc_WTd3_WTBMP_D6.txt", header = TRUE)
row.names(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]
raw_counts <- as.matrix(raw_counts)
data <- raw_counts[rowSums(raw_counts) >= 0, ]

meta <- read.table("info_WTd3_WTBMP_D6.txt", header = TRUE)
rownames(meta) <- meta$labels

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData = round(data), colData = meta, design = ~ group)

#if full library 
dds <- estimateSizeFactors(dds) 
#sf <- sizeFactors(dds)

# Supprimez les noms des facteurs de taille et imprimez les valeurs
#print(unname(sf))

#if subset library, pre-check your size factor on the full, and input it there
#sizeFactors(dds) <- c(0.4961683, 1.2883907, 0.7661422 ,2.2231175)

rld <- rlog(dds, blind=TRUE)

dds <- DESeq(dds)
res <- results(dds)

  write.csv(res, file="DE_WTD3_WTD6BMP.csv")

pdf(file = "PCA_plot_all_day.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5)


plotPCA(rld, intgroup="group", ntop=500)

dev.off()

# Plot dispersion

plotDispEsts(dds)

# Build results table

res <- results(dds)
#Output DESeq table
write.csv(res, file="DE_results_D5.csv")


summary(res$log2FoldChange, -log10(res$padj))

# plotMA
# Pay attention, if something like plotMA is available on other packages, it'll go for these other packages.
# here we wanted the deseq one, we could indicate DESeq2::plotMA...
DESeq2::plotMA(res, ylim=c(-2,2))

#volcano plot
#plot(re)


# --------------------------------
# Functional analysis with clusterProfiler

# Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(rownames(res))

# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_genes <- as.character(rownames(signif_res))

# Run GO enrichment analysis 
ego <- enrichGO(gene = signif_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                "org.Mm.eg.db", 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

# Dotplot 
dotplot(ego)

# Enrichmap
ego2 <- pairwise_termsim(ego)

#pdf(file = "emap_plot.pdf",   # The directory you want to save the file in
#    width = 10, # The width of the plot in inches
#    height = 10)

emapplot(ego2)

#dev.off()
# Category Netplot
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         vertex.label.font=6)




