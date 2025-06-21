# ============================== 1. Library ====================================
library('tidyverse')

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
temp_df_list = list()

# Read all files within list and append into list
for (i in 1:length(zebrafish_brain_file_names)) {
  file_name = gsub(".multiple.count$", "", basename(zebrafish_brain_file_names[[i]]))
  temp_df_list[[i]] <- read.table(zebrafish_brain_file_names[[i]], col.names = c("gene", file_name))
}

# Merge all files within the list to make a count_data_df
count_data_df <- temp_df_list %>% reduce(full_join, by = 'gene')

# Export count data
write.csv(count_data_df, file.path(work_dir, 'Data/zebrafish_brain_count_data.csv'), row.names = FALSE) 

# Remove the list to save space
rm(temp_df, temp_df_list)


# ====================== 4. Making Meta Data DF ================================
# Prepare a list of meta data
meta_df_list = list()

# Read all files within list and append into the list
for (i in 1:length(zebrafish_brain_file_names)) {
  temp_file_name = gsub(".multiple.count$", "", basename(zebrafish_brain_file_names[[i]]))
  file_components = strsplit(temp_file_name, "_")[[1]]
  
  meta_df_list[[i]] <- data.frame(
    sample_ID = basename(zebrafish_brain_file_names[[i]]),
    species = file_components[2],
    tissue = file_components[3],
    time_month = file_components[4],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Convert list to dataframe
meta_df_list <- do.call(rbind, meta_df_list)



# ============== Meta Data ===================
# file_names <- list.files("~/MultiSpeciesTranscriptomics/raw_counts")
# 
# metadata_list = list()
# 
# for (i in 1:length(file_names)) {
#   file_components <- strsplit(file_names[i], '_')[[1]]
#   
#   ID = file_components[1]
#   species = file_components[2]
#   tissue = file_components[3]
#   time = file_components[4]
#   time_unit = file_components[5]
#   
#   # Append as a data frame row
#   metadata_list[[i]] <- data.frame(
#     file_name = file_names[i],
#     ID = ID,
#     species = species,
#     tissue = tissue,
#     time = time,
#     time_unit = time_unit,
#     stringsAsFactors = FALSE
#   )
# }
# 
# # Combine list into a full data frame
# metadata_df <- do.call(rbind, metadata_list)
# print(metadata_df)

