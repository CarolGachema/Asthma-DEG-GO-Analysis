############################################## 
# Asthma Differential Gene Expression Analysis
# Dataset: GSE4302 (Affymetrix HG-UI33 Plus 2.0)
# Comparison: Placebo vs Flovent (Fluticasone Propionate)
# Method: Limma + GO Enrichment Analysis 
#################################################



# Load Required Libraries
library(GEOquery)       # For downloading and loading GEO datasets
library(limma)          # For differential expression analysis
library(hgu133plus2.db) # For Affymetrix probe to gene symbol annotation
library(AnnotationDbi)  # For gene ID mapping
library(ggplot2)        # For volcano plot visualization
library(ggrepel)        # For non-overlapping gene labels on volcano plot
library(clusterProfiler)# For GO enrichment analysis
library(org.Hs.eg.db)   # Human genome annotation database
library(enrichplot)     # For GO enrichment visualizations


# Set Working Directory
setwd("C:/Users/Caroline/Documents/MY PROJECTS/Asthma_Project")


# Load GEO Dataset

# Load the pre-downloaded series matrix file from NCBI GEO
gse <- getGEO(
  filename = "data/GSE4302_series_matrix.txt.gz",
  getGPL = FALSE   # Skip platform annotation download
)


# Extract Expression Data and Sample Metadata
expression_data <- exprs(gse)    # Gene expression matrix (rows = probes, cols = samples)
metadata <- pData(gse)           # Sample metadata (treatment groups, patient info, etc.)


# Preview available sample characteristics to identify group labels
unique(metadata$characteristics_ch1)


# Define Experimental Groups 

# Classify samples as either "Flovent" (treatment) or "Placebo" (control) based on the characteristics_ch1 metadata field
group <- ifelse(
  grepl("Flovent", metadata$characteristics_ch1, ignore.case = TRUE),
  "Flovent",
  ifelse(
    grepl("Placebo", metadata$characteristics_ch1, ignore.case = TRUE),
    "Placebo",
    NA   # Assign NA to any samples that don't match either group
  )
)


# Filter to Keep Only Flovent and Placebo Samples

# Remove any samples that were assigned NA (not Flovent or Placebo)
keep <- !is.na(group)
expression_data <- expression_data[, keep]
group <- group[keep]


# Convert group to a factor with Placebo as the reference level
group <- factor(group, levels = c("Placebo", "Flovent"))


# Verify sample counts per group
table(group)
message(paste("Total samples:", length(group)))
message(paste("Placebo samples:", sum(group == "Placebo")))
message(paste("Flovent samples:", sum(group == "Flovent")))


# Build the Design Matrix

# Create a no-intercept design matrix for two-group comparison
# Each column represents one group (Placebo or Flovent)
design <- model.matrix(~0 + group)
colnames(design) <- c("Placebo", "Flovent")

# Define the Contrast of Interest

# Define Placebo minus Flovent as our contrast
# Positive logFC = higher in Placebo | Negative logFC = higher in Flovent
contrast <- makeContrasts(Placebo - Flovent, levels = design)


# Fit Linear Model and Apply Empirical Bayes Statistics
fit  <- lmFit(expression_data, design)    # Fit linear model to expression data
fit2 <- contrasts.fit(fit, contrast)      # Apply the defined contrast
fit2 <- eBayes(fit2)                      # Apply empirical Bayes moderation for robust statistics


# Extract Differential Expression Results

# Retrieve statistics for all genes with FDR-adjusted p-values
# adjust.method = "fdr" applies Benjamini-Hochberg correction
results <- topTable(fit2, adjust.method = "fdr", number = Inf)


# Preview the top differentially expressed genes
head(results)


# Filter Statistically Significant Genes 

# Apply dual significance thresholds:
# FDR-adjusted p-value < 0.05 (statistically significant after correction)
sig_genes <- results[
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
]


# Report total number of significant genes identified
message(paste("Total significant DEGs:", nrow(sig_genes)))


# Annotate Probes with Gene Symbols

# Add a column containing probe IDs for mapping purposes
results$ID <- rownames(results)


# Translate Affymetrix probe IDs into human-readable HGNC gene symbols
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(results),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# Add gene symbols as a new column in the results dataframe
results$GeneSymbol <- gene_symbols


# Report how many probes were successfully annotated
num_matches <- sum(!is.na(results$GeneSymbol))
print(paste("Success! Number of genes annotated:", num_matches))


# Clean Annotated Results

# Remove probes with no matching gene symbol
# Unannotated probes cannot be biologically interpreted
annotated_results <- results[!is.na(results$GeneSymbol) & results$GeneSymbol != "", ]


# View the top 10 annotated differentially expressed genes
top10_annotated <- annotated_results[order(annotated_results$adj.P.Val), ][1:10, ]
print(top10_annotated[, c("GeneSymbol", "logFC", "adj.P.Val")])


# Save annotated results
write.csv(annotated_results, "results/all_genes_annotated.csv")
write.csv(top10_annotated,   "results/top10genes_annotated.csv")


# Generate ggplot2 Volcano Plot

# Classify all genes by expression status for colour coding
results$Expression <- "Not Significant"
results$Expression[results$adj.P.Val < 0.05 & results$logFC > 1]  <- "Higher in Placebo"
results$Expression[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Higher in Flovent"


# Select top 10 genes to label — top 5 from each direction
top5_placebo <- head(annotated_results[annotated_results$logFC > 1, ][
  order(annotated_results[annotated_results$logFC > 1, ]$adj.P.Val), ], 5)
top5_flovent <- head(annotated_results[annotated_results$logFC < -1, ][
  order(annotated_results[annotated_results$logFC < -1, ]$adj.P.Val), ], 5)
top10_label  <- rbind(top5_placebo, top5_flovent)


# Add labels to main results table
results$Label <- NA
results$Label[rownames(results) %in% rownames(top10_label)] <-
  top10_label$GeneSymbol[match(
    rownames(results)[rownames(results) %in% rownames(top10_label)],
    rownames(top10_label)
  )]


# Build the volcano plot
asthma_volcano <- ggplot(results, aes(
  x     = logFC,
  y     = -log10(adj.P.Val),
  color = Expression
)) +
  
  # Plot all gene points
  geom_point(alpha = 0.5, size = 1.8) +
  
  # Colour scheme: orange = higher in Placebo, green = higher in Flovent
  scale_color_manual(values = c(
    "Not Significant"  = "grey75",
    "Higher in Placebo"= "#E67E22",
    "Higher in Flovent"= "#27AE60"
  )) +
  
  # Significance threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Non overlapping gene labels
  geom_label_repel(
    aes(label = Label),
    na.rm         = TRUE,
    size          = 3.5,
    fontface      = "bold",
    box.padding   = 0.4,
    point.padding = 0.3,
    max.overlaps  = 20,
    color         = "black",
    fill          = "white",
    segment.color = "grey50"
  ) +
  
  # Titles and labels
  labs(
    title    = "Asthma Treatment Response: Placebo vs Flovent (Fluticasone Propionate)",
    subtitle = "Dataset: GSE4302 (Affymetrix HG-U133 Plus 2.0)  |  limma + FDR correction",
    caption  = "Positive logFC = higher in Placebo  |  Negative logFC = higher in Flovent",
    x        = "Log2 Fold Change",
    y        = "-log10 Adjusted P-value",
    color    = "Expression Status"
  ) +
  
  # Theme
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    plot.caption     = element_text(size = 9,  color = "grey40"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Save volcano plot 
ggsave(
  "figures/volcano_plot_upgraded.png",
  asthma_volcano,
  width  = 11,
  height = 8,
  dpi    = 300
)

message("Upgraded asthma volcano plot saved!!")


# Save All DEG Results to CSV
write.csv(results,    "results/all_genes.csv")             # All genes with full statistics
write.csv(sig_genes,  "results/significant_genes.csv")     # All significant DEGs

message("DEG results saved!!")



#############################################
# GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# Purpose: Identify biological processes enriched among significant DEGs
# Tool: clusterProfiler
###############################################


# Prepare Gene List for Enrichment 

# Extract gene symbols of significant DEGs for GO analysis
# Using the same thresholds: adj.P.Val < 0.05 and |logFC| > 1
genes_to_test <- annotated_results$GeneSymbol[
  annotated_results$adj.P.Val < 0.05 & abs(annotated_results$logFC) > 1
]

message(paste("Genes submitted for GO enrichment:", length(genes_to_test)))

# Run GO Enrichment Analysis (Biological Process)

# enrichGO tests whether specific GO Biological Process terms are significantly over-represented among our significant DEGs compared to the entire genome background
ego <- enrichGO(
  gene          = genes_to_test,
  OrgDb         = org.Hs.eg.db,   # Human annotation database
  keyType       = "SYMBOL",       # Input identifiers are gene symbols
  ont           = "BP",           # Biological Process ontology
  pAdjustMethod = "fdr",          # Multiple testing correction
  pvalueCutoff  = 0.05,           # P-value significance threshold
  qvalueCutoff  = 0.05            # Q-value (FDR) significance threshold
)


# Preview top 10 enriched biological processes
print("Top 10 GO Biological Processes:")
print(head(ego, 10))


# Visualize GO Enrichment Results

# Barplot: top 10 most significantly enriched GO terms
# Bar length = gene count | Colour = adjusted p-value significance
png("figures/go_enrichment_barplot.png", width = 900, height = 700)
barplot(
  ego,
  showCategory = 10,
  title        = "Top 10 GO Biological Processes — Asthma (Placebo vs Flovent)"
)
dev.off()


# Dotplot: visualization showing gene ratio and significance
# Dot size = gene ratio | Dot colour = adjusted p-value
png("figures/go_enrichment_dotplot.png", width = 900, height = 700)
dotplot(
  ego,
  showCategory = 10,
  title        = "GO Enrichment Dotplot — Asthma (Placebo vs Flovent)"
)
dev.off()

message("GO enrichment plots saved!!")


# Save GO Enrichment Results
write.csv(as.data.frame(ego), "results/go_enrichment_results.csv")

message("GO enrichment results saved!!")


#  Save Workspace
save.image("results/asthma_workspace.RData")
message("Workspace saved!!")
