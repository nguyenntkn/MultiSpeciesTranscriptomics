# ================== Library =================
library('tidyverse')

# ============ Work directory ================
# setwd('')

# ============== Meta Data ===================
file_names <- list.files("/Users/nguyennguyen/Documents/Clinical Bioinfo/Analytic and Storytelling/MultiSpeciesTranscriptomics/raw_counts")

metadata_list = list()

for (i in 1:length(file_names)) {
  file_components <- strsplit(file_names[i], '_')[[1]]
  
  ID = file_components[1]
  species = file_components[2]
  tissue = file_components[3]
  time = file_components[4]
  time_unit = file_components[5]
  
  # Append as a data frame row
  metadata_list[[i]] <- data.frame(
    file_name = file_names[i],
    ID = ID,
    species = species,
    tissue = tissue,
    time = time,
    time_unit = time_unit,
    stringsAsFactors = FALSE
  )
}

# Combine list into a full data frame
metadata_df <- do.call(rbind, metadata_list)
print(metadata_df)






