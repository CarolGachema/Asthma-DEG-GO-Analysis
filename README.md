# Asthma-DEG-GO-Analysis
Differential gene expression and GO enrichment analysis of asthma patients (Placebo vs Flovent) using GSE4302 microarray data. Pipeline built in R with limma and clusterProfiler.


## Overview

This project performs **differential gene expression (DEG) analysis** and **Gene Ontology (GO) enrichment analysis** on publicly available asthma microarray data, comparing gene expression profiles between **Placebo** and **Flovent-treated** (fluticasone propionate) asthma patients. The goal is to identify genes and biological pathways significantly modulated by corticosteroid treatment,  providing molecular insight into the mechanisms of asthma therapy.

Asthma affects over 260 million people globally and remains a leading cause of chronic respiratory morbidity. Inhaled corticosteroids like fluticasone propionate are the cornerstone of asthma management, yet the precise transcriptomic mechanisms underlying their therapeutic effects remain an active area of investigation.


## Dataset

| Parameter | Details |
|-----------|---------|
| **GEO Accession** | [GSE4302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4302) |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **Comparison** | Placebo vs Flovent (fluticasone propionate) treated asthma patients |
| **Source** | NCBI Gene Expression Omnibus (GEO) |


## Methods

```
Raw GEO Data → Expression Matrix Extraction → Group Classification (Placebo/Flovent)
→ Design Matrix → Contrast Definition → limma Linear Model
→ Empirical Bayes Statistics → FDR Correction → DEG Filtering
→ Gene Annotation → GO Enrichment Analysis (clusterProfiler)
→ Visualisation → Results Export
```

### Statistical Framework
- **DEG Package:** `limma` (Linear Models for Microarray Analysis)
- **Moderation:** Empirical Bayes (eBayes)
- **Multiple Testing Correction:** Benjamini-Hochberg FDR
- **Significance Thresholds:** adj.P.Val < 0.05 AND |log2FC| > 1
- **GO Enrichment Package:** `clusterProfiler`
- **Ontology:** Biological Process (BP)


This project incorporates **Gene Ontology (GO) enrichment analysis** which is a computational approach that explains the biological processes that those genes are collectively involved in.

This two-layer approach provides:
- **Gene-level insight** - specific transcripts modulated by Flovent treatment
- **Pathway-level insight** - the broader biological processes disrupted in asthma and restored by corticosteroid therapy


## Key Results

| Metric | Value |
|--------|-------|
| **Comparison** | Placebo vs Flovent-treated asthma patients |
| **Statistical Method** | limma with FDR correction |
| **GO Enrichment** | Biological Process ontology |

Full results available in `results/` folder



## Visualisations

### Volcano Plot
Visualises effect size (log2 fold change) against statistical significance (−log10 adjusted p-value) for the Placebo vs Flovent comparison. Genes crossing both threshold lines represent significant DEGs.

 `figures/volcano_plot.png`

### GO Enrichment Barplot
Shows the top 10 most significantly enriched Gene Ontology Biological Process terms among the differentially expressed genes. Bar length represents gene count within each GO term.

 `figures/go_enrichment_barplot.png`

### GO Enrichment Dotplot
Provides richer visualisation of GO enrichment results — dot size represents gene ratio (proportion of DEGs in the GO term) while dot colour represents adjusted p-value significance level.

 `figures/go_enrichment_dotplot.png`



## Repository Structure

```
Asthma-DEG-GO-Analysis/
│
├── data/
│   └── GSE4302_series_matrix.txt.gz         # Raw GEO series matrix (download from NCBI)
│
├── scripts/
│   └── asthma_analysis.R                    # Full annotated R analysis pipeline
│
├── results/
│   ├── all_genes.csv                        # All genes with full statistics
│   ├── significant_genes.csv               # Filtered significant DEGs
│   ├── all_genes_annotated.csv             # All genes with HGNC gene symbols
│   ├── top10genes_annotated.csv            # Top 10 most significant annotated DEGs
│   └── go_enrichment_results.csv           # Full GO enrichment results table
│
└── figures/
    ├── volcano_plot.png                     # Volcano plot (Placebo vs Flovent)
    ├── go_enrichment_barplot.png           # Top 10 GO biological processes barplot
    └── go_enrichment_dotplot.png           # GO enrichment dotplot
```

---

## How to Reproduce This Analysis

### Prerequisites
Install the following R packages before running the analysis:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery",
  "limma",
  "hgu133plus2.db",
  "AnnotationDbi",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))
```

### Steps
1. Download the GSE4302 dataset from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4302) and place it in the `data/` folder
2. Open `scripts/asthma_analysis.R` in RStudio
3. Run the script end to end
4. Results will be saved to `results/` and `figures/`


## Biological Context

**Flovent (fluticasone propionate)** is a synthetic corticosteroid that binds to glucocorticoid receptors in airway cells, broadly suppressing inflammatory gene expression. It is one of the most widely prescribed inhaled corticosteroids for asthma management globally.

This analysis investigates the transcriptomic signature of Flovent treatment, identifying which genes are differentially expressed between treated and untreated asthma patients, and which biological pathways are significantly modulated by corticosteroid therapy.

Key biological explorations:
- The inflammatory pathways are suppressed by Flovent treatment.
- Whether immune response genes significantly downregulated in treated patients.
- What GO biological processes are enriched among differentially expressed genes?


## Limitations

- Exploratory analysis - findings require functional validation
- Sample size and patient heterogeneity may affect statistical power
- GO enrichment analysis is annotation-dependent and may miss novel pathways
- Future analyses could incorporate KEGG pathway analysis for additional mechanistic depth

## Author

**Caroline Gachema**


*Data sourced from NCBI GEO — publicly available for research use.*
