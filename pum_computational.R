# Loading libraries
library(tximport)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(tidyr)

# Setting working directory
setwd("/Users/josemartinez/pum_computational")

# Upload annotation from GRCm39 gencode v34
tx2gene <- read.csv("./annotation/tx2gene.csv")
g2gene <-  read.csv("./annotation/g2gene.csv")

# Analysis of Pum Knockout RNA-Seq (GSE95102) with limma-voom

config <- read.table("config_knockdown.txt", header=TRUE)
config$condition <- factor(config$condition, levels = c("WT","DCKO","P1KO","P2KO"))
rownames(config) <- config$sample
config[,c("sample","condition")]
files <- file.path("salmon", config$sample, "quant.sf")
names(files) <- config$sample
file.exists(files)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

y <- DGEList(counts = txi$counts)  
design <- model.matrix(~ condition, data = config)

keep <- rowSums(y$counts >= 20) >= 3  # Genes must have at least 20 counts in at least 3 samples (minimum number of biological replicates)
y$counts <- y$counts[keep, ]  # Filter counts based on 'keep'
y <- calcNormFactors(y)  # Recalculate normalization factors after filtering

v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)

data <- v$E  # This extracts the logCPM matrix for heatmap
condition_mapping <- with(config, setNames(condition, sample))
colnames(data) <- condition_mapping[colnames(data)]

pheatmap(data, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", 
         color = colorRampPalette(c("blue", "white", "red"))(255),
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 10
)

config <- config[config$sample %in% colnames(y), ]
colors <- setNames(rainbow(length(levels(config$condition))), levels(config$condition))
plotMDS(v, col=colors[config$condition], pch=20)
legend("topleft", legend=names(colors), col=colors, pch=20, cex=0.8) # PCA plot

# results for DCKO
contrast.matrix <- makeContrasts(DCKO_vs_WT = conditionDCKO, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust = "BH", number = Inf)

res_DCKO <- results
res_DCKO$gene_id <- row.names(res_DCKO)
res_DCKO <- merge(res_DCKO, g2gene, by="gene_id")
res_DCKO <- res_DCKO[order(res_DCKO[,6],decreasing=FALSE),]

# results for P1KO
contrast.matrix <- makeContrasts(P1KO_vs_WT = conditionP1KO, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust = "BH", number = Inf)

res_P1KO <- results
res_P1KO$gene_id <- row.names(res_P1KO)
res_P1KO <- merge(res_P1KO, g2gene, by="gene_id")
res_P1KO <- res_P1KO[order(res_P1KO[,6],decreasing=FALSE),]

# results for P2KO
contrast.matrix <- makeContrasts(P2KO_vs_WT = conditionP2KO, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust = "BH", number = Inf)

res_P2KO <- results
res_P2KO$gene_id <- row.names(res_P2KO)
res_P2KO <- merge(res_P2KO, g2gene, by="gene_id")
res_P2KO <- res_P2KO[order(res_P2KO[,6],decreasing=FALSE),]

# Writing results table
write.csv(res_DCKO,"./results_voom/voom_DCKO.csv")
write.csv(res_P1KO,"./results_voom/voom_P1KO.csv")
write.csv(res_P2KO,"./results_voom/voom_P2KO.csv")

# Volcano Plots
df <- res_DCKO
df$color = "grey"
df$color[df$adj.P.Val < 0.05 & df$logFC > 1] = "blue"
df$color[df$adj.P.Val < 0.05 & df$logFC < -1] = "red"
count_red <- sum(df$color == "red")
count_blue <- sum(df$color == "blue")

p <- ggplot(df, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c("grey"="grey", "blue"="blue", "red"="red")) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 7)) + 
  labs(title = "Pum1/2 dCKO", x="log2 Fold Change", y="-log10 padj") +
  theme_classic() +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size=14, color="black"), 
        axis.title = element_text(size=14, color="black"),
        legend.position = "none") +
  annotate("text", x = -1.5, y = 7, size=5, label = paste("Down:", count_red), color="red", hjust=1) +
  annotate("text", x = 1.5, y = 7, size=5, label = paste("Up:", count_blue), color="blue", hjust=0)

print(p)

ggsave("./figures_voom/volcano_CDKO.pdf", plot = p, width = 4, height = 4, device = cairo_pdf)

df <- res_P1KO
df$color = "grey"
df$color[df$adj.P.Val < 0.05 & df$logFC > 1] = "blue"
df$color[df$adj.P.Val < 0.05 & df$logFC < -1] = "red"
count_red <- sum(df$color == "red")
count_blue <- sum(df$color == "blue")

p <- ggplot(df, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c("grey"="grey", "blue"="blue", "red"="red")) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 7)) + 
  labs(title = "Pum1 KO", x="log2 Fold Change", y="-log10 padj") +
  theme_classic() +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size=14, color="black"), 
        axis.title = element_text(size=14, color="black"),
        legend.position = "none") +
  annotate("text", x = -1.5, y = 7, size=5, label = paste("Down:", count_red), color="red", hjust=1) +
  annotate("text", x = 1.5, y = 7, size=5, label = paste("Up:", count_blue), color="blue", hjust=0)

print(p)

ggsave("./figures_voom/volcano_P1KO.pdf", plot = p, width = 4, height = 4, device = cairo_pdf)

df <- res_P2KO
df$color = "grey"
df$color[df$adj.P.Val < 0.05 & df$logFC > 1] = "blue"
df$color[df$adj.P.Val < 0.05 & df$logFC < -1] = "red"
count_red <- sum(df$color == "red")
count_blue <- sum(df$color == "blue")

p <- ggplot(df, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c("grey"="grey", "blue"="blue", "red"="red")) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 7)) + 
  labs(title = "Pum2 KO", x="log2 Fold Change", y="-log10 padj") +
  theme_classic() +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size=14, color="black"), 
        axis.title = element_text(size=14, color="black"),
        legend.position = "none") +
  annotate("text", x = -1.5, y = 7, size=5, label = paste("Down:", count_red), color="red", hjust=1) +
  annotate("text", x = 1.5, y = 7, size=5, label = paste("Up:", count_blue), color="blue", hjust=0)

print(p)

ggsave("./figures_voom/volcano_P2KO.pdf", plot = p, width = 4, height = 4, device = cairo_pdf)

# GSEA on log2FC ranked list
res_DCKO <- res_DCKO[order(res_DCKO[, 2], decreasing = TRUE),]
gene_list <- res_DCKO[, 2]
names(gene_list) <- res_DCKO[, 1]
gene_list <- sort(gene_list, decreasing = TRUE)
head(gene_list)

gsea <- gseGO(gene_list,
              ont = "BP",
              keyType = "ENSEMBL",
              OrgDb = "org.Mm.eg.db",
              eps = 1e-300)

gsea_df <- data.frame(gsea)

write.csv(gsea_df,"./results_voom/DCKO_gsea.csv")

categories <- c("axon guidance", "axonogenesis","synapse assembly", "regulation of neurogenesis", 
                "oxidative phosphorylation", "mitochondrial translation",
                "cytoplasmic translation", "translation at synapse")

p <- dotplot(gsea, showCategory = categories, split=".sign") + 
  facet_grid(. ~ .sign, scales = "free_x", space = "free_x") +  # Allow free space adjustment
  scale_x_continuous(breaks = seq(0.35, 0.6, by = 0.05), 
                     limits = c(0.35, 0.6)) +  # Adjust according to your specific data ranges
  theme(
    text = element_text(family = "Arial"),  # Set global text to Arial
    axis.title = element_text(family = "Arial", size = 16),  # Set axis titles to Arial and size 16
    axis.text = element_text(family = "Arial", size = 12),  # Set axis text to Arial and size 12
    axis.text.y = element_text(family = "Arial", size = 16),  # Set y-axis text (GO terms) to Arial and size 16
    strip.text = element_text(family = "Arial", size = 20),  # Set facet title text to Arial and size 20
    panel.spacing = unit(1, "lines")) # Increase spacing between panels
  
print(p)

ggsave("./figures_voom/gsea_pum_cdko.pdf", plot = p, width = 8, height = 8, device = cairo_pdf)

# Incorporating iClip data from 

pum1_genes <- read.csv("./iclip/Pum1_iClip.txt", header = FALSE)
colnames(pum1_genes) <- "gene_name"
pum2_genes <- read.csv("./iclip/Pum2_iClip.txt", header = FALSE)
colnames(pum2_genes) <- "gene_name"
pum_genes <- merge(pum1_genes,pum2_genes) # Intersection of Pum1 and Pum2 targets

write.csv(pum_genes,"./iclip/Pum_iClip_genes.csv")

Pum1_GO <- enrichGO(pum1_genes$gene_name,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    universe = res_DCKO$gene_name)

p <- dotplot(Pum1_GO, showCategory=5) + 
  theme(
    text = element_text(family = "Arial"),  # Set global text to Arial
    axis.title = element_text(family = "Arial", size = 16),  # Set axis titles to Arial and size 14
    axis.text = element_text(family = "Arial", size = 12),  # Set axis text to Arial and size 12
    axis.text.y = element_text(family = "Arial", size = 16)  # Set y-axis text (GO terms) to Arial and size 16
  )

print(p)

ggsave("./figures_voom/go_pum1_iclip.pdf", plot = p, width = 8, height = 8, device = cairo_pdf)

Pum2_GO <- enrichGO(pum2_genes$gene_name,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    universe = res_DCKO$gene_name)

p <- dotplot(Pum2_GO, showCategory=5) + 
  theme(
    text = element_text(family = "Arial"),  # Set global text to Arial
    axis.title = element_text(family = "Arial", size = 16),  # Set axis titles to Arial and size 14
    axis.text = element_text(family = "Arial", size = 12),  # Set axis text to Arial and size 12
    axis.text.y = element_text(family = "Arial", size = 16)  # Set y-axis text (GO terms) to Arial and size 16
  )

print(p)
ggsave("./figures_voom/go_pum2_iclip.pdf", plot = p, width = 8, height = 8, device = cairo_pdf)

Pum_GO <- enrichGO(pum_genes$gene_name,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   universe = res_DCKO$gene_name)

p <- dotplot(Pum_GO, showCategory=5) + 
  theme(
    text = element_text(family = "Arial"),  # Set global text to Arial
    axis.title = element_text(family = "Arial", size = 16),  # Set axis titles to Arial and size 14
    axis.text = element_text(family = "Arial", size = 12),  # Set axis text to Arial and size 12
    axis.text.y = element_text(family = "Arial", size = 16)  # Set y-axis text (GO terms) to Arial and size 16
  )

print(p)
ggsave("./figures_voom/go_pum_iclip.pdf", plot = p, width = 8, height = 8, device = cairo_pdf)

# INCORPORTATING ANALYSIS OF CORTICAL NEURON GENE EXPRESSION ACROSS DEVELOPMENT FROM Weyn-Vanhentenryck et al.

config <- read.table("config_brain.txt", header=TRUE)
config$condition <- factor(config$condition, levels = c("E14.5","E16.5","P0","P4","P7","P15","P30","P110","21M"))
rownames(config) <- config$sample
config[,c("sample","condition")]
files <- file.path("salmon", config$sample, "quant.sf")
names(files) <- config$sample
file.exists(files)

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

y <- DGEList(counts = txi$counts)  # Initialize DGEList with counts
design <- model.matrix(~ condition, data = config)

keep <- rowSums(y$counts >= 20) >= 2  
y$counts <- y$counts[keep, ]

y <- calcNormFactors(y)  # Recalculate normalization factors after filtering
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)

data <- v$E  # This extracts the logCPM matrix for heatmap
condition_mapping <- with(config, setNames(condition, sample))
colnames(data) <- condition_mapping[colnames(data)]

heatmap_plot <- pheatmap(data, 
                         clustering_distance_rows = "euclidean", 
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete", 
                         color = colorRampPalette(c("blue", "white", "red"))(255),
                         scale = "row",
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         fontsize_row = 10,
                         fontsize_col = 12,
                         border_color = NA)

print(heatmap_plot)
ggsave("./figures_voom/heatmap_cortical_neurons.pdf", plot = heatmap_plot$gtable, width = 6, height = 8, device = cairo_pdf)

# Plotting expression of Pum targets
norm_counts <- v$E
col_data <- config

# Combine normalized counts with sample information
count_data <- as.data.frame(norm_counts)
count_data$gene_id <- rownames(count_data)
count_data_long <- gather(count_data, key = "sample", value = "expression", -gene_id)
count_data_long <- inner_join(count_data_long, col_data, by = "sample")

# Calculate mean expression by gene and condition
normalized_mean_counts <- count_data_long %>%
  group_by(gene_id, condition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE))

# Reshape data into wide format and merge with gene-to-gene mapping
normalized_counts_matrix <- spread(normalized_mean_counts, key = condition, value = mean_expression)
normalized_counts_matrix <- merge(normalized_counts_matrix, g2gene, by="gene_id")
rownames(normalized_counts_matrix) <- normalized_counts_matrix$gene_id

# Only Pum1 and Pum2
genes_of_interest <- c("Pum1", "Pum2")
normalized_counts_matrix <- filter(normalized_counts_matrix, gene_name %in% genes_of_interest)

# Reorder the factor levels for the condition column
long_data <- gather(normalized_counts_matrix, key = "condition", value = "expression", -gene_id, -gene_name)
condition_order <- c("E14.5", "E16.5", "P0", "P4", "P7", "P15", "P30", "P110", "21M")
long_data$condition <- factor(long_data$condition, levels = condition_order)

# Define your color palette for clarity and easy adjustments
color_palette <- c("Pum1" = "green", "Pum2" = "salmon")

p <- ggplot(long_data, aes(x = condition, y = expression, color = gene_name, group = gene_name)) +
  geom_line(linewidth=1) +  # Use linewidth for line thickness
  geom_point(size=2, shape=1) +  # Using shape 1, which is an open circle
  scale_color_manual(values = color_palette) +
  theme_classic(base_family="Arial") +  # Set Arial as the base font for the plot
  labs(
    title = "Expression of Pum1 and Pum2 Across Neurodevelopment",
    x = "developmental stage", 
    y = expression(log[2] ~ " mean expression"),  # Use expression() for the subscript
    color = ""  # No label for the color legend
  ) +
  theme(
    axis.text = element_text(size=14, color="black"),  # Text size and color for axes
    axis.title = element_text(size=16, color="black", margin = margin(t = 15, r = 15)),  # Text size, color, and margin for axis titles
    legend.position = "top",  # Move the legend to the top of the plot
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(color="transparent", fill="transparent"),  # Transparent legend background and border
    legend.key = element_blank()  # Remove legend key border
  )

# Print the plot
print(p)

ggsave("./figures_voom/expression_pum1_pum2_neurodevelopment.pdf", plot = p, width = 6, height = 6, device = cairo_pdf)

# Plot expression of Pum1 and Pum2 shared iClip target genes

normalized_mean_counts <- count_data_long %>%
  group_by(gene_id, condition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE))

# Reshape data into wide format and merge with gene-to-gene mapping
normalized_counts_matrix <- spread(normalized_mean_counts, key = condition, value = mean_expression)
normalized_counts_matrix <- merge(normalized_counts_matrix, g2gene, by = "gene_id")
rownames(normalized_counts_matrix) <- normalized_counts_matrix$gene_id

# Filter for genes of interest
genes_of_interest <- pum_genes$gene_name
normalized_counts_matrix <- filter(normalized_counts_matrix, gene_name %in% genes_of_interest)

# Reshape data into long format for all conditions
long_data <- gather(normalized_counts_matrix, key = "condition", value = "expression", -gene_id, -gene_name)

# Define and set the factor levels for the condition column
condition_order <- c("E14.5", "E16.5", "P0", "P4", "P7", "P15", "P30", "P110", "21M")
long_data$condition <- factor(long_data$condition, levels = condition_order)

# Calculate the mean expression for each condition
mean_expression_by_condition <- long_data %>%
  group_by(condition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>%
  mutate(condition = factor(condition, levels = condition_order))

# Find and normalize the mean expression values relative to E14.5
mean_expression_E14_5 <- mean_expression_by_condition %>%
  filter(condition == "E14.5") %>%
  pull(mean_expression)

mean_expression_by_condition <- mean_expression_by_condition %>%
  mutate(norm_expression = mean_expression / mean_expression_E14_5)

# Define a dummy variable for legend creation
mean_expression_by_condition$legend_label <- "Pum1/2 shared mRNA targets"

# Updated ggplot
p <- ggplot(mean_expression_by_condition, aes(x = condition, y = mean_expression, group = 1, color = legend_label)) +
  geom_line(linewidth=1) +  # Use linewidth for line thickness
  geom_point(size=2, shape=1) +  # Using shape 1, which is an open circle
  theme_classic(base_family="Arial") +  # Set Arial as the base font for the plot
  labs(
    title = "Expression of Pum1 and Pum2 Targets Across Neurodevelopment",
    x = "developmental stage", 
    y = expression(log[2] ~ " mean expression"),  # Use expression() for the subscript
    color = ""  # Label for the color legend
  ) +
  scale_color_manual(values = c("Pum1/2 shared mRNA targets" = "black")) +  # Color for the legend label
  theme(
    axis.text = element_text(size=14, color="black"),  # Text size and color for axes
    axis.title = element_text(size=16, color="black", margin = margin(t = 15, r = 15)),  # Text size, color, and margin for axis titles
    legend.position = "top",  # Move the legend to the top of the plot
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(color="transparent", fill="transparent"),  # Transparent legend background and border
    legend.key = element_blank()  # Remove legend key border
  )

# Print the plot
print(p)
ggsave("./figures_voom/expression_pum_targets_neurodevelopment.pdf", plot = p, width = 6, height = 6, device = cairo_pdf)


save(Pum1_GO, Pum2_GO, Pum_GO, categories, gsea, data, res_P1KO, res_P2KO, res_DCKO, file = "pum_computational.RData")
