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
# Apparently it's more memory efficient to append to a list (then convert to df)
# rather than directly appending it into a df?
temp_count_df_list = list()

# Read all files within list and append into list
for (i in 1:length(zebrafish_brain_file_names)) {
  file_name = gsub(".multiple.count$", "", basename(zebrafish_brain_file_names[[i]]))
  temp_count_df_list[[i]] <- read.table(zebrafish_brain_file_names[[i]], col.names = c("gene", file_name))
}

# Merge all files within the list to make a count_data_df
count_data_df <- temp_count_df_list %>% reduce(full_join, by = 'gene')

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

# Export meta data
write.csv(meta_data_df, file.path(work_dir, 'Data/zebrafish_brain_meta_data.csv'))

# Remove temp list to save space, may not be necessary if data is small
rm(temp_meta_df_list)

# ===================== 4. Normalizing Count Data ==============================
# Not sure what normalizing method to use yet :(














