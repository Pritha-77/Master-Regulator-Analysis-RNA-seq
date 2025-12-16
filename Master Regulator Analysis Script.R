setwd("D:/PPARG-Adipogenesis/geo")

library(RegEnrich)
library(RTN)
library(DESeq2)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(RedeR)
library(igraph)
library(clusterProfiler)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(VennDiagram)
library(openxlsx)
library(grid)

# Read the raw count matrix
counts_raw <- read.xlsx("raw.genecount.xlsx")
                    

# Remove duplicate gene symbols (keep first occurrence)
counts_raw <- counts_raw[!duplicated(counts_raw[[1]]), ]

table(is.na(counts_raw[[1]]))

# Remove rows where SYMBOL (column 1) is NA
counts_raw <- counts_raw[!is.na(counts_raw[[1]]), ]

# Remove rows where SYMBOL is blank (just in case)
counts_raw <- counts_raw[counts_raw[[1]] != "" & counts_raw[[1]] != " ", ]

# Set gene symbols as row names
rownames(counts_raw) <- counts_raw[[1]] 

# Remove the gene symbol column
counts_raw <- counts_raw[, -1]


# Extract sample names (excluding SYMBOL column)
sample_ids <- colnames(counts_raw)

# Create conditions for each sample based on the naming pattern
conditions <- c(
  rep("Day0", 3),
  rep("Day4", 3),
  rep("Day8", 3)
)

# Build sample metadata dataframe
sample_info <- data.frame(
  sample_id = sample_ids,
  condition = conditions,
  stringsAsFactors = FALSE
)

# Print sample info
sample_info

# Save to file
write.table(
  sample_info,
  file = "sample_info.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Read sample information
sample_info <- read.table("sample_info.txt", 
                          header = TRUE,
                          stringsAsFactors = FALSE)

# Ensure sample order matches between count matrix and metadata
sample_info <- sample_info[match(colnames(counts_raw), sample_info$sample_id), ]

# Verify the match
all(colnames(counts_raw) == sample_info$sample_id)


#__________________Filtering Low-Expressed Genes___________________#
# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = counts_raw,
  colData = sample_info,
  design = ~ condition
)

# Filter genes: keep genes expressed in at least 80% of samples
perc_keep <- 0.8
gene_keep <- rowSums(counts_raw > 0) >= ceiling(perc_keep * ncol(counts_raw))
dds <- dds[gene_keep, ]

# Also create filtered count table for other analyses
count_tbl_filtered <- counts_raw[gene_keep, ]

#Normalizing Data for RTN Analysis
# Perform variance stabilizing transformation for RTN
vsd <- vst(dds, blind = FALSE)
normalized_counts <- assay(vsd)

#_______________Master Regulator Analysis with RegEnrich________________#
#Preparing Data for RegEnrich
#Load transcription factors from RegEnrich or a list of TFs specific to your dataset
#------------------- Load TFs -------------------#
data(TFs)
TFs_mouse_style <- stringr::str_to_title(TFs$TF_name)

#------------------- edgeR + voom -------------------#

dge <- DGEList(counts = count_tbl_filtered, samples = sample_info)
dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot = TRUE)

# Intersect TFs with expression matrix
mouse_tfs <- intersect(
  TFs_mouse_style,
  rownames(dge_v$E)
)

#------------------- Design matrix -------------------#
dge_v$targets$condition <- factor(
  dge_v$targets$condition,
  levels = c("Day0", "Day4", "Day8")
)

# RegEnrich expects 'group'
dge_v$targets$group <- dge_v$targets$condition

design_mtx <- model.matrix(~ 0 + condition, data = dge_v$targets)
colnames(design_mtx) <- levels(dge_v$targets$condition)

#------------------- Single contrast (REQUIRED) -------------------#
ctrst_Day8_vs_Day0 <- makeContrasts(
  Day8_vs_Day0 = Day8 - Day0,
  levels = design_mtx
)

#------------------- RegEnrich -------------------#

colData_regen <- data.frame(
  row.names = rownames(dge_v$targets),
  group = dge_v$targets$group,
  condition = dge_v$targets$condition
)

regen_obj <- RegenrichSet(
  expr = dge_v$E,
  colData = colData_regen,
  reg = mouse_tfs,
  method = "limma",
  design = design_mtx,
  contrast = ctrst_Day8_vs_Day0,
  networkConstruction = "COEN",
  enrichTest = "GSEA"
)

regen_obj <- regenrich_diffExpr(regen_obj)
regen_obj <- regenrich_network(regen_obj, RsquaredCut = 0.6)
regen_obj <- regenrich_enrich(regen_obj)
regen_obj <- regenrich_rankScore(regen_obj)

#------------------- Results -------------------#
res_regenrich <- na.omit(
  as.data.frame(results_score(regen_obj))
)

head(res_regenrich)

fwrite(res_regenrich, "res_mra_regenrich.csv")

#Network Reconstruction

# Create RTN object
rtni <- tni.constructor(
  expData = normalized_counts,
  regulatoryElements = mouse_tfs
)

# Compute mutual information between TFs and all genes
rtni <- tni.permutation(rtni, nPermutations = 100)

# Bootstrap analysis for network stability
rtni <- tni.bootstrap(rtni)

# Apply DPI filter to remove indirect interactions
rtni <- tni.dpi.filter(rtni)

# Check how many regulons were successfully constructed
regulon_summary <- tni.regulon.summary(rtni)



#Performing Differential Expression Analysis
# Set reference level for comparison
dds$condition <- factor(
  dds$condition,
  levels = c("Day0", "Day4", "Day8")
)

# Run DESeq2 differential expression analysis
dds_de <- DESeq(dds)

#Day4 vs Day0 (early differentiation)
res_day4_day0 <- results(
  dds_de,
  contrast = c("condition", "Day4", "Day0")
)
#Day8 vs Day0 (late differentiation / terminal state)
res_day8_day0 <- results(
  dds_de,
  contrast = c("condition", "Day8", "Day0")
)

# Sort by adjusted p-value
res_day8_day0 <- res_day8_day0[order(res_day8_day0$padj), ]
res_day8_day0_df <- as.data.frame(res_day8_day0)



#Create ranked gene list (for GSEA / pathway analysis)
logfc_ranked <- res_day8_day0_df$log2FoldChange
names(logfc_ranked) <- rownames(res_day8_day0_df)

# Remove NA values (important for GSEA)
logfc_ranked <- logfc_ranked[!is.na(logfc_ranked)]

# Rank genes
logfc_ranked <- res_day8_day0_df$log2FoldChange
names(logfc_ranked) <- rownames(res_day8_day0_df)

# Remove NA values (important for GSEA)
logfc_ranked <- logfc_ranked[!is.na(logfc_ranked)]

# Rank genes
logfc_ranked <- sort(logfc_ranked, decreasing = TRUE)


#Get top differentially expressed genes
#sig_genes <- rownames(res_day8_day0_df)[1:1000]

#OR
sig_genes <- rownames(
  subset(
    res_day8_day0_df,
    padj < 0.05 & abs(log2FoldChange) > 1
  )
)

#Method 1: Master Regulator Analysis with Fisher’s Exact Test (Day8 vs Day0)

# Prepare TNA object for Master Regulator Analysis
# phenotype: ranked log2FC values (Day8 vs Day0)
# hits: significantly DE genes driving differentiation

rtna <- tni2tna.preprocess(
  object    = rtni,
  phenotype = logfc_ranked,   # ranked DESeq2 log2FCs (Day8 vs Day0)
  hits      = sig_genes       # significant DE genes
)

# Perform MRA to identify transcription factors
# whose regulons are significantly enriched in DE genes

rtna_mra <- tna.mra(rtna)

# Retrieve all master regulator results
results_mra_fet <- tna.get(
  rtna_mra,
  what = "mra",
  ntop = -1
)

fwrite(
  results_mra_fet,
  file = "RTN_MRA_Day8_vs_Day0_FET.csv"
)

#Method 2: Master Regulator Analysis with GSEA


# Perform GSEA-based Master Regulator Analysis
# using ranked log2FC values (Day8 vs Day0)

rtna_gsea <- tna.gsea2(
  object           = rtna,
  pValueCutoff     = 0.05,
  pAdjustMethod    = "BH",
  minRegulonSize   = 15
)

# Retrieve all GSEA MRA results
results_mra_gsea <- tna.get(
  rtna_gsea,
  what = "gsea2",
  ntop = -1
)

fwrite(
  results_mra_gsea$differential,
  file = "RTN_MRA_GSEA_Day8_vs_Day0.csv"
)



#Comparing Results Between Methods

#------------------- COLLECT MASTER REGULATORS -------------------#

mra_list <- list(
  
  ## RegEnrich (GSEA-based enrichment score)
  RegEnrich = unique(
    res_regenrich$reg[
      res_regenrich$negLogPEnrich > 1.3   # ~ p < 0.05
    ]
  ),
  
  ## RTN – Fisher Exact Test MRA
  RTN_FET = unique(
    results_mra_fet$Regulon[
      results_mra_fet$Adjusted.Pvalue < 0.05
    ]
  ),
  
  ## RTN – GSEA-based MRA
  RTN_GSEA = unique(
    results_mra_gsea$differential$Regulon[
      results_mra_gsea$differential$Adjusted.Pvalue < 0.05
    ]
  )
)

# Inspect overlap sizes
lapply(mra_list, length)




venn.diagram(
  x = mra_list,
  filename = "Venn_Master_Regulators_Day8_vs_Day0.png",
  output = TRUE,
  
  lwd = 2,
  lty = "blank",
  
  fill = c("#440154FF", "#21908DFF", "#FDE725FF"),
  
  cex = 1.4,
  fontface = "bold",
  fontfamily = "sans",
  
  cat.cex = 1.6,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 135),
  cat.dist = c(0.06, 0.06, 0.08),
  
  rotation = 1
)


#high-confidence master regulators
common_MRs <- Reduce(intersect, mra_list)

common_MRs
length(common_MRs)

write.table(
  common_MRs,
  file = "Common_Master_Regulators_All3_Day8_vs_Day0.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



