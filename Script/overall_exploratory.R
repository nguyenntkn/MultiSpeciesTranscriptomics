
# ============================== 1. Library ====================================
library('tidyverse')


# ========================== 2. Work directory =================================
# Change the work directory appropriately
work_dir = '/Users/nguyennguyen/Documents/Clinical Bioinfo/Analytic and Storytelling/MultiSpeciesTranscriptomics'
setwd(work_dir)


# ================== 3. Making Count Data DF & Meta Data DF ====================
# NOTE: Still figuring out how to make naming scheme more "universal" or easily editable
# in case we want to look at other species or tissue.

# Get zebrafish files
file_names <- list.files(file.path(work_dir, "Data/raw_counts/"),
                                   pattern = "dre.*\\.multiple\\.count$",
                                   full.names = TRUE)

# Prepare a count data list and meta data list of all data files
# Apparently it's more memory efficient to append to a list (then convert to df)
# rather than directly appending it into a df?
temp_count_df_list = list()
temp_meta_df_list = list()


# Read all files within list and append into list
for (i in 1:length(file_names)) {
  
  # Shorten file names
  temp_file_name = gsub("\\.multiple\\.count$", "", basename(file_names[[i]]))

  # Append count data into count data list
  temp_count_df_list[[i]] <- read.table(file_names[[i]], col.names = c("gene", temp_file_name))
  
  # Append meta data into list, formatted with appropriate information based on file names
  file_components = strsplit(temp_file_name, "_")[[1]]
  temp_meta_df_list[[i]] <- data.frame(
    sample_ID = temp_file_name,
    species = file_components[2],
    tissue = file_components[3],
    time = paste0(file_components[4], file_components[5]),
    row.names = temp_file_name,
    stringsAsFactors = FALSE
  )
}


# Merge all files within the list to make a count_data_df
zebrafish_count_data_df <- temp_count_df_list %>% 
  purrr::reduce(full_join, by = 'gene') %>% 
  column_to_rownames('gene')
zebrafish_meta_data_df <- do.call(rbind, temp_meta_df_list)


# Export count data and meta data
write.csv(zebrafish_count_data_df, file.path(work_dir, 'Data/zebrafish_count_data.csv'), row.names = FALSE) 
write.csv(zebrafish_meta_data_df, file.path(work_dir, 'Data/zebrafish_meta_data.csv'))


# Temp list is quite large, remove to save storage space
rm(temp_count_df_list, temp_meta_df_list)


# ------------------------------------------------------------------------------
# Now do the same thing for killifish
# Get killifish files
file_names <- list.files(file.path(work_dir, "Data/raw_counts/"),
                                   pattern = "nfu.*\\.multiple\\.count$",
                                   full.names = TRUE)

# Prepare a count data list and meta data list of all data files
temp_count_df_list = list()
temp_meta_df_list = list()


# Read all files within list and append into list
for (i in 1:length(file_names)) {
  
  # Shorten file names
  temp_file_name = gsub("_mzm_tophat2_mapped\\.bam\\.multiple\\.count$", "", basename(file_names[[i]]))
  
  # Append count data into count data list
  temp_count_df_list[[i]] <- read.table(file_names[[i]], col.names = c("gene", temp_file_name))
  
  # Append meta data into list, formatted with appropriate information based on file names
  file_components = strsplit(temp_file_name, "_")[[1]]
  temp_meta_df_list[[i]] <- data.frame(
    sample_ID = temp_file_name,
    species = file_components[2],
    tissue = file_components[3],
    time = paste0(file_components[4], file_components[5]),
    row.names = temp_file_name,
    stringsAsFactors = FALSE
  )
}


# Merge all files within the list to make a count_data_df
killifish_count_data_df <- temp_count_df_list %>% 
  purrr::reduce(full_join, by = 'gene') %>% 
  column_to_rownames('gene')
killifish_meta_data_df <- do.call(rbind, temp_meta_df_list)


# Export count data and meta data
write.csv(killifish_count_data_df, file.path(work_dir, 'Data/killifish_count_data.csv'), row.names = FALSE) 
write.csv(killifish_meta_data_df, file.path(work_dir, 'Data/killifish_meta_data.csv'))


# Temp list is quite large, remove to save storage space
rm(temp_count_df_list, temp_meta_df_list)



# =========================== 4. Exploratory ===================================

any(is.na(zebrafish_count_data_df)) # FALSE: No missing data

dim(zebrafish_count_data_df) # number of rows and columns
head(zebrafish_count_data_df, 5) # to see first 5 rows
summary(zebrafish_count_data_df) # summary statistics
colSums(is.na(zebrafish_count_data_df)) # to confirm that this is no missing data
sum_zero=colSums(zebrafish_count_data_df == 0) # number of rows where expression = zero
sum_not_zero=colSums(zebrafish_count_data_df != 0) # number of rows per column where expression is over zero

print(mean(sum_zero))
print(mean(sum_not_zero))

# Getting the number of samples per tissue
table(zebrafish_meta_data_df$tissue)

# Spread of the non-zero counts and less than 100 counts 
long_counts <- zebrafish_count_data_df %>%
  pivot_longer(
    cols = everything(),  # All sample columns (assuming no gene name column)
    names_to = "sample",
    values_to = "express"
  ) %>%
 filter(express > 0 & express < 100) 

# Plotting histogram 
ggplot(long_counts, aes(x = express)) +
  geom_histogram(bins = 100) +
  labs(
    title = "Histogram of Raw Gene Expression (<100)",
    x = "Expression Level",
    y = "Frequency"
  ) +
  theme_minimal()



