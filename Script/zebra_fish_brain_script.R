# ============================== 1. Library ====================================
library('tidyverse')
library('DESeq2')
library('pheatmap')
library('msigdbr')
library('fgsea')
library('ggplot2')
library('dplyr')

# ========================== 2. Work directory =================================
# Change the work directory appropriately
work_dir = '/Users/nguyennguyen/Documents/Clinical Bioinfo/Analytic and Storytelling/MultiSpeciesTranscriptomics'
setwd(work_dir)

# ====================== 3. Making Count Data DF ===============================
# NOTE: Still figuring out how to make naming scheme more "universal" or easily editable
# in case we want to look at other species or tissue.

# Get zebrafish brain files
zebrafish_brain_file_names <- list.files(file.path(work_dir, "Data/raw_counts/"),
                                         pattern = "dre_brain_.*\\.multiple\\.count$",
                                         full.names = TRUE)

# Prepare a list of all data files
# Apparently it's more memory efficient to append to a list (then convert to df)
# rather than directly appending it into a df?
temp_count_df_list = list()

# Read all files within list and append into list
for (i in 1:length(zebrafish_brain_file_names)) {
  file_name = gsub(".multiple.count$", "", basename(zebrafish_brain_file_names[[i]]))
  temp_count_df_list[[i]] <- read.table(zebrafish_brain_file_names[[i]], col.names = c("gene", file_name))
}

# Merge all files within the list to make a count_data_df
count_data_df <- temp_count_df_list %>% 
  purrr::reduce(full_join, by = 'gene') %>% 
  column_to_rownames('gene')

# Export count data
write.csv(count_data_df, file.path(work_dir, 'Data/zebrafish_brain_count_data.csv'), row.names = FALSE) 

# Temp list is quite large, remove to save storage space
rm(temp_count_df_list)


# ====================== 4. Making Meta Data DF ================================
# Prepare a list of meta data
temp_meta_df_list = list()

# Read all files within list and append into the list
for (i in 1:length(zebrafish_brain_file_names)) {
  temp_file_name = gsub("\\.multiple\\.count$", "", basename(zebrafish_brain_file_names[[i]]))
  file_components = strsplit(temp_file_name, "_")[[1]]
  
  temp_meta_df_list[[i]] <- data.frame(
    sample_ID = temp_file_name,
    species = file_components[2],
    tissue = file_components[3],
    time = file_components[4],
    time_unit = file_components[5],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Convert list to dataframe
meta_data_df <- do.call(rbind, temp_meta_df_list)
rownames(meta_data_df) <- meta_data_df$sample_ID

# Export meta data
write.csv(meta_data_df, file.path(work_dir, 'Data/zebrafish_brain_meta_data.csv'))

# Remove temp list to save space, may not be necessary if data is small
rm(temp_meta_df_list)

# ======================== 4. Data Exploration =================================

dim(count_data_df) # number of rows and columns
head(count_data_df) # to see first 6 rows
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
  filter(express > 0 & express < 100) 
  

ggplot(long_counts, aes(x = express)) +
  geom_histogram(bins = 100) +
  labs(
    title = "Histogram of Raw Gene Expression (<100)",
    x = "Expression Level",
    y = "Frequency"
  ) +
  theme_minimal()

# ======================== 5. DESeq2 analysis ==================================

# Make DESeqDataSet object for DESeq
dds <- DESeqDataSetFromMatrix(
  countData = count_data_df,
  colData = meta_data_df,
  design = ~ time)

# Perform DEA
dds <- DESeq(dds)
res <- results(dds)

# Volcano plot of DESeq2 results
res_df <- as.data.frame(res) %>% # matrix out of results
  rownames_to_column("gene") %>%
  mutate(significant = padj < 0.03)

# volcano plot
# x - >0 upregulated, <0 downregulated
# y - statistical significance
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant, shape = significant)) +
  geom_point(alpha = 0.6) + 
  xlim(c(-3, 3)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot") +
  theme_minimal()

perform_DESeq2 <- function(dds, condition1, condition2) {
  res <- results(dds, contrast=c("time", condition1, condition2))
  res_df <- as.data.frame(res) %>% # converts to data frame
    rownames_to_column("gene") %>% # adds gene names
    # adds column with T/F if gene is statistically significant
    mutate(significant = padj < 0.05, 
           condition = ifelse(log2FoldChange > 0, # gene upregulated
                              condition1, # in condition 1 (a gene in one condition)
                              condition2)) # otherwise condition 2 (a gene in other condition)
  res_df <- na.omit(res_df) # removes incomplete genes
  # volcano plot
  ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant, shape=condition)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values = c("turquoise", "red")) +
    scale_shape_manual(values = c(16, 17)) +
    labs(title = paste("Volcano Plot:", condition1, "vs", condition2),
         x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "top")
}

# Plot volcano plots for each condition pair
plot1 <- perform_DESeq2(dds, "12", "24")
plot2 <- perform_DESeq2(dds, "12", "36")
plot3 <- perform_DESeq2(dds, "12", "42")

plot1
plot2
plot3


# =============== 6. GSEA - Gene Set Enrichment Analysis =======================

### Load gene sets for GSEA
# MSigDB = Molecular Signatures Database
# Curated collection of gene sets used for GSEA
zebrafish_gene_sets <- msigdbr(species = "Danio rerio", 
                               category = "C5") # C5 - gene ontology terms
# result: df from gene names and ensembl gene IDs

# converts df to named list to use in fgsea()
# Each list element = vector of genes for one gene set
gene_sets_list <- split(zebrafish_gene_sets$ensembl_gene, 
                        zebrafish_gene_sets$gs_name)

### Rank genes for GSEA analysis

genes_ranked <- -log10(res$pvalue) * sign(res$log2FoldChange)
names(genes_ranked) <- rownames(res)
genes_ranked <- genes_ranked[!is.na(genes_ranked) & is.finite(genes_ranked)]

### Perform GSEA
set.seed(42)
fgseaRes <- fgsea(pathways = gene_sets_list, 
                  stats = genes_ranked,
                  minSize = 15, # min no of genes a set should have to be tested
                  maxSize = 500)

head(fgseaRes[order(pval), ])

# Enrichment plot for a single pathway
plotEnrichment(gene_sets_list[["GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX"]],
               genes_ranked) + labs(title="GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX")

# Enrichment plot for several pathways (Top 10 pathways up vs Top 10 pathways down)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gene_sets_list[topPathways], genes_ranked, fgseaRes, 
              gseaParam=0.5)

