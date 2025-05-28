###########Figure 2: NPC-RSS Model Construction and Validation###########
#Figure 2A: Model Performance Heatmap
# Load required libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Load model performance data
model_performance <- read.csv("data/processed/model_performance.csv", row.names=1)

# Create heatmap of model performance
p <- ggplot(melt(model_performance), aes(X2, X1, fill=value)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  geom_text(aes(label=sprintf("%.3f", value)), size=3) +
  theme_minimal() +
  labs(title="Model Performance (AUC) for Different Algorithm Combinations",
       x="Classification Methods", 
       y="Feature Selection Methods",
       fill="AUC") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot
ggsave("figures/figure2A_model_heatmap.png", p, width=16, height=12, dpi=300)

#Figure 2B: Gene Weight Coefficients
# Load required libraries
library(ggplot2)

# Gene weights from NPC-RSS model
gene_weights <- data.frame(
  gene = c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1', 'C1QTNF3', 'CA11', 
           'CDK5RAP3', 'CLDN1', 'EYA1', 'IFI44L', 'KREMEN1', 'NT5DC2', 
           'NTRK3', 'RFX4', 'RHOBTB3', 'SLC1A2', 'TRIM58'),
  weight = c(22.6640, 14.3984, 7.6856, 8.1306, 10.0473, 2.7369, 3.8292, 
             -15.9445, -11.2475, -6.2474, 1.6897, 0.9464, -6.5371, 
             -10.4192, -0.5755, -2.1728, -7.1597, -0.3858)
)

# Filter genes with absolute weight > 10
high_weight_genes <- gene_weights[abs(gene_weights$weight) > 10, ]

# Sort genes by absolute weight
high_weight_genes <- high_weight_genes[order(-abs(high_weight_genes$weight)), ]

# Create bar chart
p <- ggplot(high_weight_genes, aes(x=reorder(gene, -abs(weight)), y=weight, fill=weight > 0)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("red", "green"), guide="none") +
  geom_hline(yintercept=0, linetype="solid", size=0.5, color="black") +
  labs(title="NPC-RSS Genes with Weight Coefficients > 10",
       x="Genes", 
       y="Weight Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot
ggsave("figures/figure2B_gene_weights.png", p, width=12, height=8, dpi=300)

#Figure 2C: ROC Curves
# Load required libraries
library(ggplot2)
library(pROC)

# Load validation data
training_results <- read.csv("data/processed/training_roc_data.csv")
validation_results <- read.csv("data/processed/validation_roc_data.csv")

# Create a data frame for plotting
roc_data <- rbind(
  data.frame(fpr=training_results$fpr, tpr=training_results$tpr, Dataset="Training"),
  data.frame(fpr=validation_results$fpr, tpr=validation_results$tpr, Dataset="Validation")
)

# Plot ROC curves
p <- ggplot(roc_data, aes(x=fpr, y=tpr, color=Dataset)) +
  geom_line(size=1.2) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray") +
  scale_color_manual(values=c("Training"="blue", "Validation"="red"),
                     labels=c(paste0("Training Set (AUC = ", round(training_results$auc[1], 3), ")"),
                              paste0("Validation Set (AUC = ", round(validation_results$auc[1], 3), ")"))) +
  labs(title=paste0("ROC Curves (Weighted AUC = 0.932)"),
       x="False Positive Rate", 
       y="True Positive Rate") +
  theme_minimal() +
  theme(legend.position="bottom") +
  coord_fixed(ratio=1, xlim=c(0, 1), ylim=c(0, 1))

# Save plot
ggsave("figures/figure2C_roc_curves.png", p, width=10, height=8, dpi=300)

#Figure 2D: Volcano Plot
# Load required libraries
library(ggplot2)
library(ggrepel)

# Load differential expression data
deg_data <- read.csv("data/processed/differential_expression.csv")

# Add significance column
deg_data$significant <- deg_data$pvalue < 0.05 & abs(deg_data$log2FoldChange) > 0.585

# Create volcano plot
p <- ggplot(deg_data, aes(x=log2FoldChange, y=-log10(pvalue), color=significant)) +
  geom_point(alpha=0.6, size=3) +
  scale_color_manual(values=c("FALSE"="gray", "TRUE"="red"), guide=FALSE) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") +
  geom_vline(xintercept=c(-0.585, 0.585), linetype="dashed", color="gray") +
  labs(title="Volcano Plot of Differentially Expressed Genes",
       x="Log2 Fold Change", 
       y="-Log10 P-value") +
  theme_minimal()

# Highlight key genes
key_genes <- c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1')
key_data <- deg_data[deg_data$gene %in% key_genes, ]
if(nrow(key_data) > 0) {
  p <- p + geom_text_repel(data=key_data, 
                           aes(label=gene), 
                           size=4,
                           box.padding=0.5,
                           point.padding=0.5)
}

# Save plot
ggsave("figures/figure2D_volcano_plot.png", p, width=12, height=10, dpi=300)

#Figure 2F-I: In vitro Validation
# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# Load gene expression data from CNE2 cell lines
cne2_data <- read.csv("data/processed/cne2_expression.csv")

# Figure A: Top 5 genes expression 
top_genes <- c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1')

# Prepare data for plotting
plot_list <- list()

for(gene in top_genes) {
  # Extract gene expression data
  gene_data <- data.frame(
    Expression = cne2_data[[gene]],
    Group = cne2_data$group
  )
  
  # Calculate means and standard errors
  gene_stats <- aggregate(Expression ~ Group, data=gene_data, 
                          FUN=function(x) c(mean=mean(x), se=sd(x)/sqrt(length(x))))
  gene_stats <- do.call(data.frame, gene_stats)
  names(gene_stats) <- c("Group", "Mean", "SE")
  
  # Perform t-test
  t_test <- t.test(Expression ~ Group, data=gene_data)
  p_val <- t_test$p.value
  
  # Determine significance level
  if(p_val < 0.001) sig <- "***"
  else if(p_val < 0.01) sig <- "**"
  else if(p_val < 0.05) sig <- "*"
  else sig <- "ns"
  
  # Create plot
  p <- ggplot(gene_stats, aes(x=Group, y=Mean, fill=Group)) +
    geom_bar(stat="identity", width=0.6) +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2) +
    scale_fill_manual(values=c("CNE2-P"="green", "CNE2-RS"="red")) +
    labs(title=gene, y="Expression Level") +
    theme_minimal() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1)) +
    annotate("text", x=1.5, y=max(gene_stats$Mean + gene_stats$SE) * 1.1, 
             label=sig, size=5)
  
  plot_list[[gene]] <- p
}

# Arrange plots in a grid
combined_plot <- grid.arrange(grobs=plot_list, ncol=5)

# Save plot
ggsave("figures/figure2F_top_genes_expression.png", combined_plot, width=20, height=6, dpi=300)

# Figure 2G: Heatmap of DEGs
# Load DEG expression data
deg_expr <- read.csv("data/processed/deg_expression_matrix.csv", row.names=1)

# Create heatmap
png("figures/figure2G_deg_heatmap.png", width=12, height=10, units="in", res=300)
pheatmap(as.matrix(deg_expr), 
         scale="none", 
         cluster_rows=TRUE, 
         cluster_cols=FALSE,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),
         fontsize_row=10,
         fontsize_col=10)
dev.off()

# Figure 2H: Heatmap of NPC-RSS genes (z-score)
# Get NPC-RSS genes expression data
gene_list <- c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1', 'C1QTNF3', 'CA11', 
               'CDK5RAP3', 'CLDN1', 'EYA1', 'IFI44L', 'KREMEN1', 'NT5DC2', 
               'NTRK3', 'RFX4', 'RHOBTB3', 'SLC1A2', 'TRIM58')
rss_expr <- deg_expr[gene_list, ]

# Z-score normalize
rss_expr_z <- t(scale(t(rss_expr)))

# Create heatmap
png("figures/figure2H_rss_zscore_heatmap.png", width=12, height=8, units="in", res=300)
pheatmap(as.matrix(rss_expr_z), 
         scale="none", 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),
         main="Z-scores of NPC-RSS Genes in CNE2 Cell Lines",
         fontsize_row=10,
         fontsize_col=10)
dev.off()

# Figure 2I: NPC-RSS scores comparison
# Load NPC-RSS scores
score_data <- read.csv("data/processed/npc_rss_scores_cne2.csv")

# Rename groups
score_data$group <- ifelse(score_data$group == "CNE2-P", "Sensitive", "Resistant")

# Perform t-test
t_test_result <- t.test(npc_rss_score ~ group, data=score_data)
p_val <- t_test_result$p.value

# Create significance annotation
if(p_val < 0.001) sig <- "***"
else if(p_val < 0.01) sig <- "**"
else if(p_val < 0.05) sig <- "*"
else sig <- "ns"

# Create boxplot
p <- ggplot(score_data, aes(x=group, y=npc_rss_score, fill=group)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.1, size=3, alpha=0.7) +
  scale_fill_manual(values=c("Sensitive"="green", "Resistant"="red")) +
  labs(title="Comparison of NPC-RSS Scores in CNE2 Cell Lines",
       y="NPC-RSS Score") +
  theme_minimal() +
  theme(legend.position="none",
        axis.title.x=element_blank()) +
  annotate("segment", x=1, xend=2, 
           y=max(score_data$npc_rss_score)*1.1, 
           yend=max(score_data$npc_rss_score)*1.1) +
  annotate("text", x=1.5, y=max(score_data$npc_rss_score)*1.15, 
           label=sig, size=5)

# Save plot
ggsave("figures/figure2I_npc_rss_scores.png", p, width=6, height=8, dpi=300)
