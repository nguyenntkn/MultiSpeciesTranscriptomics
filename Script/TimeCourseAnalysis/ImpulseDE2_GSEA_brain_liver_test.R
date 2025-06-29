# ========================== 1. Library & Setup ================================
# Bioconductor doesn't have ImpulseDE2 so we'll install it straight from the source
# First, install remotes if you don't have it
# install.packages("remotes") 
# Then try installing ImpulseDE2 from GitHub
# remotes::install_github("YosefLab/ImpulseDE2") 

library("tidyverse")
library("ImpulseDE2")
library("clusterProfiler")
library("org.Dr.eg.db") # Load the organism database for zebrafish



# Set your work directory (assuming you've defined work_dir previously)
work_dir = '/Users/nguyennguyen/Documents/Clinical Bioinfo/Analytic and Storytelling/MultiSpeciesTranscriptomics'
setwd(work_dir)

# ============================ 2. Data Loading =================================

# Get zebrafish files
file_names <- list.files(file.path(work_dir, "Data/raw_counts/"),
                         pattern = "dre.*\\.multiple\\.count$",
                         full.names = TRUE)

temp_count_df_list = list()
temp_meta_df_list = list()

# Read all files within list and append into list
for (i in 1:length(file_names)) {
  
  # Shorten file names
  temp_file_name = gsub("\\.multiple\\.count$", "", basename(file_names[[i]]))
  
  # Append count data into count data list
  # Use stringsAsFactors = FALSE to prevent automatic conversion to factors for 'gene' column
  temp_count_df_list[[i]] <- read.table(file_names[[i]], col.names = c("gene", temp_file_name),
                                        stringsAsFactors = FALSE)
  
  # Append meta data into list, formatted with appropriate information based on file names
  file_components = strsplit(temp_file_name, "_")[[1]]
  temp_meta_df_list[[i]] <- data.frame(
    Sample = temp_file_name,
    species = file_components[2],
    Condition = file_components[3],
    Time = file_components[4],
    TimeCateg = file_components[5],
    row.names = temp_file_name,
    stringsAsFactors = FALSE
  )
}

# Merge all files within the list to make a count_data_df
zebrafish_count_data_df <- temp_count_df_list %>%
  purrr::reduce(full_join, by = 'gene') %>%
  column_to_rownames('gene')

# Ensure the count data is strictly numeric and integer (as expected by DE methods)
zebrafish_count_data_df <- as.matrix(zebrafish_count_data_df)
mode(zebrafish_count_data_df) <- "numeric"
zebrafish_count_data_df <- round(zebrafish_count_data_df) # Round to nearest integer

# Handle any NA's that might have been introduced if raw data had non-numeric text
if (any(is.na(zebrafish_count_data_df))) {
  warning("NA values introduced in count data during numeric conversion. Please check raw files.")
  # Optionally, you might fill NAs with 0 or remove rows/columns with NAs
  # zebrafish_count_data_df[is.na(zebrafish_count_data_df)] <- 0
}

# Make medata df
zebrafish_meta_data_df <- do.call(rbind, temp_meta_df_list)

# Filter count data for only the samples in brain vs liver meta data
# This is crucial if your full count_data_df has more samples than your brain_liver_meta
brain_liver_samples <- rownames(zebrafish_meta_data_df %>% filter(Condition %in% c("brain", "liver")))
zebrafish_count_data_subset <- zebrafish_count_data_df[, brain_liver_samples]


# Subset brain vs liver and change to control vs case for ImpulseDE2
brain_liver_meta <- zebrafish_meta_data_df %>%
  filter(Condition %in% c("brain", "liver")) %>%
  mutate(
    Condition = case_when(
      Condition == "brain" ~ "case",
      Condition == "liver" ~ "control",
      TRUE ~ Condition # Keep other conditions as they are, if any
    ),
    # Ensure Time is numeric as ImpulseDE2 expects it
    Time = as.numeric(Time)
  )

# ========================== 3. ImpulseDE2 Run =================================

# To run ImpulseDE2, we need matCountData/count_data_df and dfAnnotation/meta_data_df
# Make sure matCountData dimensions match dfAnnotation rows
zebrafish_ImpulseDE2 <- runImpulseDE2(
  matCountData    = zebrafish_count_data_subset, # Use the subsetted count data
  dfAnnotation    = brain_liver_meta,            # Meta df
  boolCaseCtrl    = TRUE,                        # Case control analysis? (False for time series)
  vecConfounders  = NULL,                        # Confounding variables/ factor to correct for batch effects
  scaNProc        = 4                            # Number of process for parallelisation
)


# =============== 4. Gene Set Enrichment Analysis (GSEA) =======================

# 1. Prepare the ranked gene list from ImpulseDE2 results
# Get the full results data frame
results_df <- zebrafish_ImpulseDE2$dfImpulseDE2Results

# Ensure 'padj' column exists and handle potential NA or 0 values in q-value
if (!"padj" %in% colnames(results_df)) {
  stop("The 'padj' (adjusted p-value) column is missing from ImpulseDE2 results.")
}

# Replace any q-values of 0 with a very small positive number to avoid -log10(0) = Inf
min_padj_value_non_zero <- min(results_df$padj[results_df$padj > 0], na.rm = TRUE)
if (min_padj_value_non_zero == Inf) { # Handle case where all q-values might be 0 or NA
  warning("All q-values are 0 or NA. Cannot create a meaningful ranking for GSEA.")
  # You might want to stop or handle this specific case based on your data.
  # For now, we'll assign a very small default value.
  min_padj_value_non_zero = 1e-300
}
results_df$padj[results_df$padj <= 0] <- min_padj_value_non_zero / 10


# 2. Convert Gene Symbols to Entrez IDs
# GSEA functions in clusterProfiler generally prefer Entrez IDs
# 'Gene' column in dfImpulseDE2Results is assumed to be gene symbols (from your row names)
gene_map <- bitr(results_df$Gene,
                 fromType = "ENSEMBL", # Assuming your gene identifiers are Symbols
                 toType = "ENTREZID",
                 OrgDb = org.Dr.eg.db,
                 drop = TRUE) # drop = TRUE will remove genes that couldn't be mapped

# Merge Entrez IDs with your results and calculate the ranking score
ranked_genes_df <- results_df %>%
  inner_join(gene_map, by = c("Gene" = "ENSEMBL")) %>%
  # Calculate the ranking score: -log10(adj_p_value)
  mutate(rank_score = -log10(padj)) %>%
  # Group by ENTREZID and keep the row with the maximum rank_score
  group_by(ENTREZID) %>%
  summarize(rank_score = max(rank_score, na.rm = TRUE), .groups = 'drop') %>%
  dplyr::select(ENTREZID, rank_score) %>% # Use dplyr::select explicitly
  arrange(desc(rank_score)) # Sort in decreasing order for GSEA


# Create the named vector required by GSEA
geneList <- ranked_genes_df$rank_score
names(geneList) <- ranked_genes_df$ENTREZID

print(head(geneList))

# 3. Perform GSEA for Gene Ontology (Biological Process)
gse_go_bp <- gseGO(geneList      = geneList,
                   OrgDb         = org.Dr.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP", # Biological Process
                   pAdjustMethod = "BH", # Benjamini-Hochberg FDR correction
                   pvalueCutoff  = 0.05, # Cutoff for adjusted p-value of enriched terms
                   verbose       = FALSE)

gsea_results <- gse_go_bp@result
print(gsea_results$Description) # Get all gene sets description

# Optionally visualize
# dotplot(gse_go_bp, showCategory=15, split=".sign") + facet_grid(.~.sign)






# # For GO Biological Process results:
# if (!is.null(gse_go_bp) && nrow(as.data.frame(gse_go_bp)) > 0) {
#   go_results_df <- as.data.frame(gse_go_bp)
#   
#   # Define keywords for age-related and circadian rhythm
#   age_keywords <- c("age", "aging", "senescence", "longevity", "lifespan")
#   circadian_keywords <- c("circadian", "rhythm", "clock", "diurnal", "entrainment")
#   
#   # Combine all keywords into a single regex pattern (case-insensitive)
#   # The \\b ensures whole word matching for some keywords, but careful with "age" vs "storage"
#   # Let's start with flexible matching and refine if needed
#   all_keywords_pattern <- paste(c(age_keywords, circadian_keywords), collapse = "|")
#   all_keywords_pattern <- paste0("(?i)", all_keywords_pattern) # (?i) for case-insensitive
#   
#   # Filter out rows where the Description or ID contains any of these keywords
#   # We'll keep terms that *do not* match the pattern
#   filtered_go_results <- go_results_df %>%
#     dplyr::filter(!grepl(all_keywords_pattern, Description) &
#                     !grepl(all_keywords_pattern, ID)) # Check both Description and ID
#   
#   message("\nFiltered GO Biological Process GSEA results (first few rows):")
#   if (nrow(filtered_go_results) > 0) {
#     print(filtered_go_results$Description)
#   } else {
#     message("No GO Biological Process terms remain after filtering for age/circadian keywords.")
#   }
#   
#   # You can save the filtered results
#   # write.csv(filtered_go_results, file.path(work_dir, "GSEA_GO_BP_Filtered.csv"), row.names = FALSE)
# }






