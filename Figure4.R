###########Figure 4: Single-Cell Analysis###########
#Figure 4A: UMAP Clustering Plot
# Load required libraries
library(ggplot2)
library(Seurat)
library(viridis)
library(patchwork)

# Load single-cell clustering data (assuming a pre-processed Seurat object or similar data frame)
umap_coords <- read.csv("data/processed/umap_coordinates.csv")
cell_annotations <- read.csv("data/processed/cell_annotations.csv")

# Combine data
plot_data <- merge(umap_coords, cell_annotations, by="cell_id")

# Create UMAP plot
p <- ggplot(plot_data, aes(x=UMAP_1, y=UMAP_2, color=cluster)) +
  geom_point(size=0.5, alpha=0.7) +
  scale_color_manual(values=colorRampPalette(viridis(17))(17)) +
  labs(title="Cell Clustering by UMAP",
       color="Cell Type") +
  theme_minimal() +
  theme(legend.position="right")

# Save plot
ggsave("figures/figure4A_umap_clustering.png", p, width=12, height=10, dpi=300)

#Figure 4B: Cell Type Differences in Sensitivity Scores
# Load sensitivity score data by cell type
sensitivity_scores <- read.csv("data/processed/cell_type_sensitivity.csv")

# Create plot
p <- ggplot(sensitivity_scores, aes(x=cell_type, y=sensitivity_score, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=c("Sensitive"="blue", "Resistant"="red")) +
  labs(title="NPC-RSS Scores by Cell Type",
       x="Cell Type", 
       y="NPC-RSS Score",
       fill="Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot
ggsave("figures/figure4B_cell_type_sensitivity.png", p, width=14, height=8, dpi=300)

#Figure 4C: UMAP with Sensitivity Scores
# Load sensitivity score data mapped to cells
cell_scores <- read.csv("data/processed/cell_sensitivity_scores.csv")

# Merge with UMAP coordinates
plot_data <- merge(umap_coords, cell_scores, by="cell_id")

# Create UMAP plot colored by sensitivity score
p <- ggplot(plot_data, aes(x=UMAP_1, y=UMAP_2, color=sensitivity_score)) +
  geom_point(size=0.5, alpha=0.7) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="NPC-RSS Scores on UMAP",
       color="NPC-RSS Score") +
  theme_minimal()

# Save plot
ggsave("figures/figure4C_umap_sensitivity.png", p, width=12, height=10, dpi=300)

#Figure 4D: NPC-RSS Gene Expression in Cell Subpopulations
# Load gene expression data across cell subpopulations
gene_expr_subpop <- read.csv("data/processed/gene_expression_subpopulations.csv")

# Reshape for heatmap
expr_matrix <- reshape2::dcast(gene_expr_subpop, gene ~ cell_subpopulation, value.var="expression")
row.names(expr_matrix) <- expr_matrix$gene
expr_matrix$gene <- NULL

# Filter for NPC-RSS genes
npc_rss_genes <- c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1', 'C1QTNF3', 'CA11', 
                   'CDK5RAP3', 'CLDN1', 'EYA1', 'IFI44L', 'KREMEN1', 'NT5DC2', 
                   'NTRK3', 'RFX4', 'RHOBTB3', 'SLC1A2', 'TRIM58')
expr_matrix <- expr_matrix[npc_rss_genes, ]

# Create heatmap
png("figures/figure4D_gene_expression_heatmap.png", width=14, height=10, units="in", res=300)
pheatmap(as.matrix(expr_matrix), 
         scale="row", 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         main="NPC-RSS Model Gene Expression in Cell Subpopulations",
         fontsize_row=10,
         fontsize_col=10)
dev.off()

#Figure 4E: Cell Type Composition by Group
# Load cell type composition data
cell_composition <- read.csv("data/processed/cell_type_composition.csv")

# Create stacked bar chart
p <- ggplot(cell_composition, aes(x=group, y=percentage, fill=cell_type)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis_d() +
  labs(title="Cell Type Composition by Group",
       x="Group", 
       y="Percentage",
       fill="Cell Type") +
  theme_minimal()

# Save plot
ggsave("figures/figure4E_cell_composition.png", p, width=10, height=8, dpi=300)

#Figure 4F: Marker Gene Expression
# Load marker gene expression data
marker_genes <- read.csv("data/processed/marker_gene_expression.csv")

# Reshape for heatmap
marker_matrix <- reshape2::dcast(marker_genes, gene ~ cell_type, value.var="expression")
row.names(marker_matrix) <- marker_matrix$gene
marker_matrix$gene <- NULL

# Create heatmap
png("figures/figure4F_marker_gene_expression.png", width=14, height=12, units="in", res=300)
pheatmap(as.matrix(marker_matrix), 
         scale="row", 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         main="Marker Gene Expression",
         fontsize_row=10,
         fontsize_col=10)
dev.off()

