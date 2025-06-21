# ================== Library =================
library('tidyverse')

# ============ Work directory ================

work_dir = '/Users/nguyennguyen/Documents/Clinical Bioinfo/Analytic and Storytelling/MultiSpeciesTranscriptomics/'
setwd(work_dir)

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

# ============= Making Count Data DF =================
# Get files
zebrafish_file_names <- list.files(file.path(work_dir, "Data/zebrafish_brain/"),
                                   pattern = "\\.multiple\\.count$",
                                   full.names = TRUE)

# Prepare a list of all data files
temp_df_list = list()

# Read all files within list and append into list
for (i in 1:length(zebrafish_file_names)) {
  
  file_name = gsub(".multiple.count$", "", basename(zebrafish_file_names[[i]]))
  
  temp_df_list[[i]] <- read.table(zebrafish_file_names[[i]], 
                                  col.names = c("gene", file_name))
}

# Merge all files within the list to make a count_data_df
count_data_df <- temp_df_list %>% reduce(full_join, by = 'gene')

# Remove the list to save space
rm(temp_df, temp_df_list)




