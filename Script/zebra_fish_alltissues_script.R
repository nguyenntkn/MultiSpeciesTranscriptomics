# ============================== 0. Library ====================================
library('tidyverse')
library('DESeq2')
library('pheatmap')
library('msigdbr')
library('fgsea')
library('ggplot2')
library('dplyr')
library('umap')
library(ggrepel)
library(biomaRt)

# ============================= 1. Functions ===================================
Volcano_plot <- function(df, tissue){
  #1. Gene Annotation
  # Connect to Ensembl BioMart for zebrafish gene information
  ensembl <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl")
  
  # Extract Ensembl gene IDs from the row names of the input dataframe
  ensembl_ids <- rownames(df)
  
  # Query BioMart to get corresponding gene symbols
  gene_annotations <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Rename columns for clarity
  colnames(gene_annotations) <- c("ensembl_id", "gene_symbol")
  
  #2. Prepare Data for Volcano Plot 
  volcano_df <- as.data.frame(df) %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%     # Remove rows with NA padj or fold change
    mutate(ensembl_id = rownames(.)) %>%                  # Add Ensembl ID as a column for merging
    left_join(gene_annotations, by = "ensembl_id") %>%    # Add gene symbols
    mutate(
      neg_log10_padj = -log10(padj),                      # Calculate -log10 adjusted p-value
      significance = case_when(                           # Define significance category
        padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      label = ifelse(is.na(gene_symbol), ensembl_id, gene_symbol)  # Label with gene symbol or Ensembl ID
    )
  
  # 3. Select Genes to Label 
  # Identify top 10 most significant genes (lowest padj) to label on the plot
  top_genes <- volcano_df %>%
    filter(significance == "Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)
  
  # 4. Create Volcano Plot 
  ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.8, size = 1.5) +                                # Plot points
    geom_text_repel(data = top_genes, aes(label = label), size = 3, max.overlaps = Inf) +  # Add labels
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +    # Color coding
    coord_cartesian(xlim = c(-12, 12), ylim = c(0, 200)) +                    # Set fixed axis limits for consistency
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +  # Fold change cutoffs
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # padj cutoffs
    labs(
      title = paste("Volcano Plot -", tissue, "Tissue"),
      x = "Log2 Fold Change",
      y = "-log10 Adjusted P-value"
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),  # Add border around plot
      legend.position = "none",                                 # Remove legend
      plot.title = element_text(size = 18),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    )
}

# GSEA function for each tissue type
run_gsea <- function(res_df, tissue_name) {
  res_df <- res_df[!is.na(res_df$stat), ]  # remove NAs
  
  # Rank genes by stat
  genes_ranks <- res_df$stat
  names(genes_ranks) <- rownames(res_df)
  genes_ranks <- sort(genes_ranks, decreasing = TRUE)
  
  # Run fgsea
  set.seed(42)
  fgseaRes <- fgsea(pathways = gene_sets_list,
                    stats = genes_ranks,
                    minSize = 10,
                    maxSize = 500)
  
  # Save top results
  top_results <- fgseaRes %>% arrange(padj) %>% head(10)
  print(paste("Top GSEA results for", tissue_name))
  print(top_results[, c("pathway", "NES", "pval", "padj")])
  
  # Enrichment plot for a single pathway
  p <- plotEnrichment(gene_sets_list[["GOBP_RHYTHMIC_PROCESS"]], genes_ranks) + 
    labs(title="GOBP_RHYTHMIC_PROCESS")
  
  
  # Enrichment plot for several pathways (Top 10 pathways up vs Top 10 pathways down)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathways <- topPathwaysUp
  
  gsea_table_plot <- plotGseaTable(gene_sets_list[topPathways], genes_ranks, fgseaRes, gseaParam = 0.5)
  grid::grid.newpage()
  grid::grid.draw(gsea_table_plot)
  
  return(fgseaRes)
}
# ========================== 2. Work directory =================================
# Change the work directory appropriately
work_dir = '/Users/nguyennguyen/Documents/Clinical Bioinfo/Analytic and Storytelling/MultiSpeciesTranscriptomics'
setwd(work_dir)

# # ====================== 3. Making Count Data DF ===============================
# 
# raw_count_dir <- file.path(work_dir, "Data/raw_counts")
# 
# # Get count files for each tissue
# zebrafish_file_names <- list.files(raw_count_dir,
#                                    pattern = "dre_(brain|liver|skin)_.*\\.multiple\\.count$",
#                                    full.names = TRUE)
# 
# # ====================== 3.1. Read and Merge Count Files =====================
# 
# # Initialize empty list to store data frames
# temp_count_df_list <- list()
# 
# # Loop through all files and read them into the list
# for (i in seq_along(zebrafish_file_names)) {
#   file_path <- zebrafish_file_names[i]
#   sample_id <- gsub("\\.multiple\\.count$", "", basename(file_path))
#   
#   temp_count_df_list[[i]] <- read.table(file_path, col.names = c("gene", sample_id))
# }
# 
# # Merge all files by "gene"
# count_data_df <- temp_count_df_list %>%
#   purrr::reduce(full_join, by = "gene") %>%
#   column_to_rownames("gene")
# 
# # Export Merged Counts
# write.csv(count_data_df, file.path(work_dir, "Data/zebrafish_all_count_data.csv"), row.names = TRUE)
# 
# # Remove temp list to save space, may not be necessary if data is small
# rm(temp_count_df_list)
# 
# # ====================== 3.2. Making Meta Data DF ================================
# 
# # Prepare a list of meta data
# temp_meta_df_list <- list()
# 
# # Read all files within list and append into the list
# for (i in seq_along(zebrafish_file_names)) {
#   temp_file_name <- gsub("\\.multiple\\.count$", "", basename(zebrafish_file_names[i]))
#   file_components <- strsplit(temp_file_name, "_")[[1]]
#   
#   temp_meta_df_list[[i]] <- data.frame(
#     sample_ID = temp_file_name,
#     species = file_components[2],
#     tissue = file_components[3],
#     time = file_components[4],
#     time_unit = file_components[5],
#     stringsAsFactors = FALSE
#   )
# }
# 
# # Convert list to dataframe
# meta_data_df <- do.call(rbind, temp_meta_df_list)
# meta_data_df <- meta_data_df %>%
#   arrange(tissue, time)
# 
# rownames(meta_data_df) <- meta_data_df$sample_ID
# 
# # Export metadata to CSV
# write.csv(meta_data_df, file.path(work_dir, 'Data/zebrafish_all_meta_data.csv'), row.names = TRUE)
# 
# # Remove temp list to save space, may not be necessary if data is small
# rm(temp_meta_df_list)



# ======================== 4. Data Exploration =================================

count_data_df <- read.csv(file.path(work_dir, "Data/zebrafish_all_count_data.csv"), row.names = "X")
meta_data_df <- read.csv(file.path(work_dir, "Data/zebrafish_all_meta_data.csv"), row.names = "X")

# Organize count_data
ordered_samples <- meta_data_df$sample_ID
count_data_df <- count_data_df[, ordered_samples]

meta_data_df$time <- as.character(meta_data_df$time)


dim(count_data_df) # number of rows and columns
head(count_data_df) # to see first 5 rows
summary(count_data_df) # summary statistics
colSums(is.na(count_data_df)) # to confirm that this is no missing data
sum_zero=colSums(count_data_df == 0) # number of rows where expression = zero
sum_not_zero=colSums(count_data_df != 0) # number of rows per column where expression is over zero

# Spread of the non-zero counts

long_counts <- count_data_df %>%
  pivot_longer(
    cols = everything(),  # All sample columns (assuming no gene name column)
    names_to = "sample",
    values_to = "express"
  )%>%
  filter(express > 0 & express < 1000) 
  

ggplot(long_counts, aes(x = express)) +
  geom_histogram(bins = 1000) +
  labs(
    title = "Histogram of Raw Gene Expression (<1000)",
    x = "Expression Level",
    y = "Frequency"
  ) +
  theme_minimal()

# ======================== 5. DESeq2 analysis ==================================

# === For Brain Tissue ===

meta_brain <- filter(meta_data_df, tissue == "brain")
counts_brain <- count_data_df[, meta_brain$sample_ID]
dds_brain <- DESeqDataSetFromMatrix(countData = counts_brain,
                                    colData = meta_brain,
                                    design = ~ time)
dds_brain <- DESeq(dds_brain, test = "LRT", reduced = ~1)
res_brain <- results(dds_brain) %>% as.data.frame() %>% arrange(padj)

# === For Liver Tissue ===

meta_liver <- filter(meta_data_df, tissue == "liver")
counts_liver <- count_data_df[, meta_liver$sample_ID]
dds_liver <- DESeqDataSetFromMatrix(countData = counts_liver,
                                    colData = meta_liver,
                                    design = ~ time)
dds_liver <- DESeq(dds_liver, test = "LRT", reduced = ~1)
res_liver <- results(dds_liver) %>% as.data.frame() %>% arrange(padj)

# === For Skin Tissue ===

meta_skin <- filter(meta_data_df, tissue == "skin")
counts_skin <- count_data_df[, meta_skin$sample_ID]
dds_skin <- DESeqDataSetFromMatrix(countData = counts_skin,
                                   colData = meta_skin,
                                   design = ~ time)
dds_skin <- DESeq(dds_skin, test = "LRT", reduced = ~1)
res_skin <- results(dds_skin) %>% as.data.frame() %>% arrange(padj)

# === For all 3 Tissues ===
dds <- DESeqDataSetFromMatrix(countData = count_data_df,
                              colData = meta_data_df,
                              design = ~ tissue + time)
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds) %>% as.data.frame() %>% arrange(padj)

# ------------ Volcano Plot ---------
Volcano_plot(res_brain, tissue = "Brain")
Volcano_plot(res_liver, tissue = "Liver")
Volcano_plot(res_skin, tissue = "Skin")

# =============== 6. GSEA - Gene Set Enrichment Analysis =======================

### Load gene sets for GSEA
# MSigDB = Molecular Signatures Database
# Curated collection of gene sets used for GSEA

zebrafish_gene_sets <- msigdbr(species = "Danio rerio", category = "C5")
gene_sets_list <- split(zebrafish_gene_sets$ensembl_gene, zebrafish_gene_sets$gs_name)

fgsea_brain <- run_gsea(res_brain, "brain")
fgsea_liver <- run_gsea(res_liver, "liver")
fgsea_skin <- run_gsea(res_skin, "skin")

# ========================== 7. UMAP ===========================================

vsd <- vst(dds, blind = TRUE) # or use rlog(dds)

# Get vst-transformed expression matrix
vsd_mat <- assay(vsd)  # genes x samples
vsd_mat <- t(vsd_mat)  # transpose to samples x genes

# Run UMAP
set.seed(42)
umap_result <- umap(vsd_mat)  # Default settings work well for most RNA-seq

# Format UMAP result for ggplot
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sample_ID <- rownames(umap_df)

# Merge with metadata
umap_df <- left_join(umap_df, meta_data_df, by = "sample_ID")

# Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = tissue, shape = time)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(title = "UMAP of VST-transformed RNA-seq data",
       subtitle = "Colored by tissue, shaped by age (months)")

# =========================== 8. PCA ===========================================

# Run PCA on log-transformed counts for better scaling
vsd <- vst(dds, blind = TRUE) # or use rlog(dds)
vsd_mat <- assay(vsd)
log_counts_filtered <- vsd_mat[apply(vsd_mat, 1, var) > 0, ]
pca <- prcomp(t(log_counts_filtered), scale. = TRUE)

# Extract sample coordinates
pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df <- merge(pca_df, meta_data_df, by.x = "sample", by.y = "row.names")

# Plot with ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue, shape = time)) +
  geom_point(size = 5) +
  labs(title = "PCA of Gene Expression", x = "PC1", y = "PC2") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )


