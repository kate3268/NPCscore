###########Figure 3: Immune-Related Analysis###########
#Figure 3A: Immune Cell Infiltration Comparison
# Load required libraries
library(ggplot2)
library(reshape2)

# Load immune cell data
immune_cell_data <- read.csv("data/processed/immune_cell_infiltration.csv")

# Reshape for plotting
immune_long <- melt(immune_cell_data, 
                    id.vars="cell_type", 
                    measure.vars=c("sensitive", "resistant"),
                    variable.name="group", 
                    value.name="infiltration")

# Calculate p-values for each cell type
p_values <- sapply(unique(immune_cell_data$cell_type), function(cell) {
  cell_data <- immune_cell_data[immune_cell_data$cell_type == cell, ]
  t.test(cell_data$sensitive, cell_data$resistant)$p.value
})

# Create significance labels
sig_labels <- sapply(p_values, function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
})

# Add significance data
sig_data <- data.frame(
  cell_type = names(p_values),
  p_value = p_values,
  sig = sig_labels
)

# Merge with long format data
immune_long$cell_type <- factor(immune_long$cell_type, levels=unique(immune_cell_data$cell_type))

# Create plot
p <- ggplot(immune_long, aes(x=cell_type, y=infiltration, fill=group)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("sensitive"="blue", "resistant"="red"),
                    labels=c("Sensitive", "Resistant")) +
  labs(title="Immune Cell Infiltration in NPC Tissue Groups",
       x="", y="Infiltration Score",
       fill="Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Add significance stars for cell types with p < 0.05
sig_data_filtered <- sig_data[sig_data$sig != "ns", ]
if(nrow(sig_data_filtered) > 0) {
  max_heights <- aggregate(infiltration ~ cell_type, data=immune_long, max)
  sig_data_filtered <- merge(sig_data_filtered, max_heights, by="cell_type")
  
  p <- p + geom_text(data=sig_data_filtered, 
                     aes(x=cell_type, y=infiltration*1.1, label=sig),
                     position=position_dodge(width=0.9))
}

# Save plot
ggsave("figures/figure3A_immune_infiltration.png", p, width=14, height=8, dpi=300)

#Figure 3B: Correlation between NPC-RSS key genes and immune cells
# Load required libraries
library(ggplot2)
library(reshape2)

# Load correlation data between genes and immune cells
gene_immune_corr <- read.csv("data/processed/gene_immune_correlation.csv")

# Convert to long format for plotting
corr_long <- melt(gene_immune_corr, 
                  id.vars="gene", 
                  variable.name="immune_cell", 
                  value.name="correlation")

# Load p-values
p_values <- read.csv("data/processed/gene_immune_pvalues.csv")
p_long <- melt(p_values, 
               id.vars="gene", 
               variable.name="immune_cell", 
               value.name="pvalue")

# Combine correlation and p-values
plot_data <- merge(corr_long, p_long, by=c("gene", "immune_cell"))

# Filter for top 5 genes
plot_data <- plot_data[plot_data$gene %in% c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1'), ]

# Add significance indicator
plot_data$significant <- plot_data$pvalue < 0.05
plot_data$point_size <- -log10(plot_data$pvalue)

# Create bubble plot
p <- ggplot(plot_data, aes(x=immune_cell, y=gene, 
                           size=point_size, 
                           color=correlation)) +
  geom_point(alpha=0.7) +
  scale_color_gradient2(low="blue", mid="white", high="orange", midpoint=0) +
  scale_size_continuous(range=c(0, 8), guide=FALSE) +
  labs(title="Correlation between NPC-RSS key genes and immune cells",
       x="Immune Cell Type", 
       y="Gene",
       color="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot
ggsave("figures/figure3B_gene_immune_correlation.png", p, width=14, height=8, dpi=300)

#Figure 3C: Immune cell interactions
# Load required libraries
library(ggplot2)
library(reshape2)
library(corrplot)

# Load immune cell correlation data
immune_corr <- read.csv("data/processed/immune_cell_correlations.csv", row.names=1)

# Create correlation plot
png("figures/figure3C_immune_correlations.png", width=10, height=10, units="in", res=300)
corrplot(as.matrix(immune_corr), 
         method="circle", 
         type="upper", 
         order="hclust",
         tl.col="black", 
         tl.srt=45,
         tl.cex=0.7,
         col=colorRampPalette(c("blue", "white", "red"))(100),
         p.mat = matrix(0.05, nrow=ncol(immune_corr), ncol=ncol(immune_corr)),
         sig.level=0.05,
         insig="blank",
         main="Immune Cell Type Correlations")
dev.off()

#Figure 3D: Correlation of NPC-RSS key genes with immune features
# Load required libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Load TISIDB correlation data
tisidb_data <- read.csv("data/processed/tisidb_correlations.csv")

# Reshape data for heatmap
tisidb_matrix <- reshape2::dcast(tisidb_data, gene ~ feature, value.var="correlation")
row.names(tisidb_matrix) <- tisidb_matrix$gene
tisidb_matrix$gene <- NULL

# Filter for top 5 genes
top_genes <- c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1')
tisidb_matrix <- tisidb_matrix[top_genes, ]

# Load p-values
p_values <- read.csv("data/processed/tisidb_pvalues.csv")
p_matrix <- reshape2::dcast(p_values, gene ~ feature, value.var="pvalue")
row.names(p_matrix) <- p_matrix$gene
p_matrix$gene <- NULL
p_matrix <- p_matrix[top_genes, ]

# Create annotation for significance
stars_matrix <- matrix("", nrow=nrow(p_matrix), ncol=ncol(p_matrix))
rownames(stars_matrix) <- rownames(p_matrix)
colnames(stars_matrix) <- colnames(p_matrix)

for(i in 1:nrow(p_matrix)) {
  for(j in 1:ncol(p_matrix)) {
    if(p_matrix[i,j] < 0.001) stars_matrix[i,j] <- "***"
    else if(p_matrix[i,j] < 0.01) stars_matrix[i,j] <- "**"
    else if(p_matrix[i,j] < 0.05) stars_matrix[i,j] <- "*"
  }
}

# Create annotation for types of immune features
feature_types <- data.frame(
  Feature_Type = rep(c("Chemokine", "MHC", "Immunostimulator", "Immunoinhibitor", "Receptor"), 
                     each=ceiling(ncol(tisidb_matrix)/5))
)
feature_types <- feature_types[1:ncol(tisidb_matrix), , drop=FALSE]
rownames(feature_types) <- colnames(tisidb_matrix)

# Create heatmap
png("figures/figure3D_tisidb_correlations.png", width=14, height=8, units="in", res=300)
pheatmap(as.matrix(tisidb_matrix), 
         scale="none", 
         cluster_rows=FALSE, 
         cluster_cols=TRUE,
         annotation_col=feature_types,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),
         main="Correlation of NPC-RSS key genes with immune features",
         fontsize_row=10,
         fontsize_col=8,
         display_numbers=stars_matrix)
dev.off()


