###########Figure 5: Pathway Analysis###########
#Figure 5A-B: GSVA Analysis
# Load required libraries
library(ggplot2)

# Load GSVA results for SMARCA2
gsva_smarca2 <- read.csv("data/processed/gsva_smarca2.csv")

# Create barplot for SMARCA2
ggplot(gsva_smarca2, aes(x=reorder(pathway, NES), y=NES, fill=NES > 0)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("FALSE"="blue", "TRUE"="red"), guide=FALSE) +
  labs(title="GSVA of SMARCA2",
       x="", 
       y="Normalized Enrichment Score") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5)) +
  coord_flip() +
  ggsave("figures/figure5A_gsva_smarca2.png", width=12, height=8, dpi=300)

# Load GSVA results for DMC1
gsva_dmc1 <- read.csv("data/processed/gsva_dmc1.csv")

# Create barplot for DMC1
ggplot(gsva_dmc1, aes(x=reorder(pathway, NES), y=NES, fill=NES > 0)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("FALSE"="blue", "TRUE"="red"), guide=FALSE) +
  labs(title="GSVA of DMC1",
       x="", 
       y="Normalized Enrichment Score") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5)) +
  coord_flip() +
  ggsave("figures/figure5B_gsva_dmc1.png", width=12, height=8, dpi=300)


#Figure 5C-D: GSEA Analysis
# Load required libraries
library(ggplot2)
library(enrichplot)

# Load GSEA results for SMARCA2
gsea_smarca2 <- read.csv("data/processed/gsea_smarca2.csv")

# Create dot plot for SMARCA2
ggplot(gsea_smarca2, aes(x=NES, y=reorder(pathway, NES), 
                         color=p.adjust, size=-log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  labs(title="GSEA of SMARCA2",
       x="Normalized Enrichment Score", 
       y="",
       color="Adjusted P-value",
       size="-log10(Adjusted P-value)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5)) +
  ggsave("figures/figure5C_gsea_smarca2.png", width=12, height=8, dpi=300)

# Load GSEA results for DMC1
gsea_dmc1 <- read.csv("data/processed/gsea_dmc1.csv")

# Create dot plot for DMC1
ggplot(gsea_dmc1, aes(x=NES, y=reorder(pathway, NES), 
                      color=p.adjust, size=-log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  labs(title="GSEA of DMC1",
       x="Normalized Enrichment Score", 
       y="",
       color="Adjusted P-value",
       size="-log10(Adjusted P-value)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5)) +
  ggsave("figures/figure5D_gsea_dmc1.png", width=12, height=8, dpi=300)



#Figure 5E: Correlation with Radiosensitization Genes
# Load required libraries
library(ggplot2)
library(reshape2)

# Load correlation data between NPC-RSS genes and radiosensitization genes
radio_correlation <- read.csv("data/processed/radiosensitization_correlation.csv")

# Reshape for plotting
corr_long <- melt(radio_correlation, 
                  id.vars="npc_rss_gene", 
                  variable.name="radio_gene", 
                  value.name="correlation")

# Load p-values
p_values <- read.csv("data/processed/radiosensitization_pvalues.csv")
p_long <- melt(p_values, 
               id.vars="npc_rss_gene", 
               variable.name="radio_gene", 
               value.name="pvalue")

# Combine correlation and p-values
plot_data <- merge(corr_long, p_long, by=c("npc_rss_gene", "radio_gene"))

# Filter for top 5 NPC-RSS genes
plot_data <- plot_data[plot_data$npc_rss_gene %in% 
                         c('SMARCA2', 'DMC1', 'CD9', 'PSG4', 'KNG1'), ]

# Add significance indicator
plot_data$significant <- plot_data$pvalue < 0.05
plot_data$point_size <- -log10(plot_data$pvalue)

# Create bubble plot
ggplot(plot_data, aes(x=radio_gene, y=npc_rss_gene, 
                      size=point_size, 
                      color=correlation)) +
  geom_point(alpha=0.7) +
  scale_color_gradient2(low="green", mid="white", high="orange", midpoint=0) +
  scale_size_continuous(range=c(0, 8)) +
  labs(title="Correlation with Radiosensitization-related Genes",
       x="Radiosensitization Gene", 
       y="NPC-RSS Gene",
       color="Correlation",
       size="-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggsave("figures/figure5E_radiosensitization_correlation.png", width=14, height=8, dpi=300)