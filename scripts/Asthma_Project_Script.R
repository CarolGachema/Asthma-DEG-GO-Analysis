############################################## 
# Asthma Differential Gene Expression Analysis
# Dataset: GSE4302 (Affymetrix HG-UI33 Plus 2.0)
# Comparison: Placebo vs Flovent (Fluticasone Propionate)
# Method: Limma + GO Enrichment Analysis (clusterProfiler)
#################################################



# Load Required Libraries
library(GEOquery)    # For downloading and loading GEO datasets
library(limma)       # For differential expression analysis 


# Set Working Directory
setwd("C:/Users/Caroline/Documents/MY PROJECTS/Asthma_Project")


# Load GEO Dataset
# Load the pre downloaded series matrix file from NCBI GEO
gse <- getGEO(
  filename = "data/GSE4302_series_matrix.txt.gz",
  getGPL = FALSE # Skip platform annotation download 
)


# Extract Expression Data and Sample Metadata
expression_data <- exprs(gse)   # Gene expression matrix (rows=probes, cols=samples)
metadata <- pData(gse)          # Sample metadata (treatment groups, patient info, etc)


# Preview available sample characteristics to identify group labels
unique(metadata$characteristics_ch1)


# Classify samples as either "Flovent" (treatment) or "Placebo" (control)
# based on the characteristics_ch1 metadata field
group <- ifelse(
  grepl("Flovent", metadata$characteristics_ch1, ignore.case = TRUE),
  "Flovent",
  ifelse(
    grepl("Placebo", metadata$characteristics_ch1, ignore.case = TRUE),
    "Placebo",
    NA  #Assign NA to any samples that don't match either group
  )
)


# Filter to keep only Flovent and Placebo samples
keep <- !is.na(group)
expression_data <- expression_data[, keep]
group <- group[keep]


# Convert group to a factor with Placebo as the reference level
group <- factor(group, levels = c("Placebo", "Flovent"))


# Verify sample counts per group
table(group)


# Build the design matrix

#Create a no intercept design matrix for twogroup comparison.
#Each column represents one group (Placebo or Flovent)
design <- model.matrix(~0 + group)
colnames(design) <- c("Placebo", "Flovent")


# Define contrast: Placebo minus flovent

# This identifies genes that are differentially expressed
# between untreated and treated asthma patients
contrast <- makeContrasts(Placebo - Flovent, levels = design)


# Fit Linear Model and apply Empirical Bayes Statistics
fit <- lmFit(expression_data, design)  # Fit linear model to expression data
fit2 <- contrasts.fit(fit, contrast)   # Apply the defined contrast.
fit2 <- eBayes(fit2)                   # Apply empirical Bayes moderation for robust statistics.


# Extract all genes with FDR-adjusted p-values

# adjust.method = "fdr" applies Benjamini-Hochberg correction
results <- topTable(fit2, adjust.method = "fdr", number = Inf)

# Preview the top differentially expressed genes
head(results)

# Apply significance thresholds

# Adjusted p-value < 0.05 (fdr corrected)
# |log2 fold change| >1 (at least 2 fold change in expression)
sig_genes <- results[
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
]


# Report total number of significant genes identified
nrow(sig_genes)


# Generate volcano plot

# visualizes effect side vs statistical significance
plot(
  results$logFC,
  -log10(results$adj.P.Val),
  pch = 20,
  col = "grey",
  main = "Placebo vs Flovent(Asthma Patients)",
  xlab = "Log2 Fold Change",
  ylab = "-log10 Adjusted P-value"
)

# Add threshold lines

# Horizontal blue line = significance threshold (adj.P.Val =0.05)
# Vertical blue lines = fold change thresholds (LogFC = ±1)
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)


# Save volcano plot to figures folder
png("figures/volcano_plot.png", width = 800, height = 600)


# Save Results to CSV
write.csv(results, "results/all_genes.csv")            #All genes with statistics
write.csv(sig_genes, "results/significant_genes.csv")  # Significant DEGs only


#Preview top 10 most significant genes
top10 <- head(results, 10)
top10


# Prepare for Gene Annotation

#Add a column containing the probe IDs (row names) for mapping purposes
results$ID <- rownames(results)


# Install and Load Annotation Package

# hgu133plus2.db contains probe to gene mappings for GPL570
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
library(AnnotationDbi)


# Map Probe IDs to Gene Symbols

# Translate Affymetrix probe IDs into human readable HGNC gene symbols
# multivials = "first" keeps only the first gene symbol if multipe map to one probe
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(results),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# Add gene symbols as a new column in the results dataframe
results$GeneSymbol <- gene_symbols


# Report how many probes were successfully annotated with gene symbols
num_matches <- sum(!is.na(results$GeneSymbol))
print(paste("Success! Number of genes annotated:", num_matches))


# Clean Annotated Results.

#Remove probes with no matching gene symbol
# (unannotated probes cannot be biologically interpreted)
annotated_results <- results[!is.na(results$GeneSymbol) & results$GeneSymbol != "", ]


# View the top 10 annotated differentially expressed genes
top10_annotated <- annotated_results[order(annotated_results$adj.P.Val), ][1:10, ]
print(top10_annotated[, c("GeneSymbol", "logFC", "adj.P.Val")])


# Save the final annotated results
write.csv(annotated_results, "results/all_genes_annotated.csv")
write.csv(top10_annotated, "results/top10genes_annotated.csv")



###################################
# GENE ONTOLOGY ENRICHMEN ANALYSIS
# Purpose: Identify biological processes significantly enriched among DEGs
# Tool: clusterProfiler
####################################


# Install nd Load GO Enrichment Libraries
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
library(clusterProfiler)   #Core enrichment analysis tool
library(org.Hs.eg.db)      # Human genome annotation database
library(enrichplot)        # Visualization tools for enrichment results


# Prepare Gene List for Enrichment

# Extract gene symbols of significant DEGs for GO analysis
# Using the same thresholds applied earlier: adj.P.Val < 0.05, |logFC|
genes_to_test <- annotated_results$GeneSymbol[annotated_results$adj.P.Val < 0.05 & abs(annotated_results$logFC) > 1]


# Run the GO Enrichment Analysis (Biological Process)
# enrichGO tests GO:BP terms
ego <- enrichGO(
  gene = genes_to_test,
  OrgDb = org.Hs.eg.db,    # Human annotation database
  keyType = "SYMBOL",      # Input gene identifiers are gene symbo;s
  ont = "BP",              # BP = Biological Process, CC = Cellular Component, MF = Molecular Function
  pAdjustMethod = "fdr",   # Multiple testing correction method
  pvalueCutoff = 0.05,     # Significance threshold for p-value
  qvalueCutoff = 0.05      # Significance threshold for q-value
)


# Preview the top 10 enriched biological processes
print("Top 10 GO Biological Processes:")
print(head(ego, 10))


# Visualize the GO Enrichment results

# Barplot shows the top 10 most significantly rencriched GO terms
# Bar length represents gene count; colour reps significance
png("figures/go_enrichment_barplot.png", width = 800, height = 600)
barplot(ego, showCategory = 10, title = "Top 10 GO Biological Processes in Asthma")
dev.off()


#Dot plot provides richer information including gene ratio and adjustedP value
png("figures/go_enrichment_dotplot.png", width = 800, height = 600)
dotplot(ego, showCategory = 10, title = "GO Enrichment Dotplot")
dev.off()


# Save GO Enrichment Results

# Export full GO enrichment results table for reporting and further analysis
write.csv(as.data.frame(ego), "results/go_enrichment_results.csv")
  
