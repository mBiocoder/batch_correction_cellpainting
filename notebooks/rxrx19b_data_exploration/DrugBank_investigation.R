############################################################################################################
# Goal of this script: Use misclassified_obs.csv to investigate the misclassified treatments with DrugBank.#
# Author: Mahima Arunkumar                                                                                 #
# ##########################################################################################################

library(dbdataset)
library(stringr)
library(dplyr)
library(ggplot2)


# General overview
names(dbdataset::drugbank)
head(dbdataset::drugbank[["salts"]], 5)
head(dbdataset::drugbank[["products"]], 5)
names(dbdataset::drugbank[["references"]])
names(dbdataset::drugbank[["cett"]])


# Read the CSV file
misclassified_obs <- read.csv("./misclassified_obs.csv")
misclassified_obs_harmony <- read.csv("./misclassified_harmony_obs.csv")

# Extract the treatment column
treatment <- misclassified_obs$treatment
treatment_harmony <- misclassified_obs_harmony$treatment

# Remove empty strings from the treatment column
treatment <- treatment[treatment != ""]
treatment_harmony <- treatment_harmony[treatment_harmony != ""]

# Convert treatment names to title case (capitalize first letter)
treatment <- str_to_title(treatment)
treatment_harmony <- str_to_title(treatment_harmony)


# Create dataframes
drugs <- dbdataset::drugbank[["drugs"]]
drug_names <- drugs[["general_information"]]$name

drug_df <- as.data.frame(drugs[["general_information"]])
drug_interactions_df <- as.data.frame(drugs[["drug_interactions"]])
category_df <- as.data.frame(drugs[["categories"]])
atc_codes_df <- as.data.frame(drugs[["atc_codes"]])
pharmacology_df <- as.data.frame(drugs[["pharmacology"]])

# Only contain information on misclassified treatments
subset_drug_df <- drug_df[drug_df$name %in% treatment, ] #50 out of 71 treatments found on DrungBank
subset_drug_df_harmony <- drug_df[drug_df$name %in% treatment_harmony, ] #59 out of 86 treatments found on DrungBank

drugbank_id_of_interest <- subset_drug_df$primary_key
drugbank_id_of_interest_harmony <- subset_drug_df_harmony$primary_key

#For subset_drug_df$primary_key search pharmacology_df and atc_codes_df
subset_pharmacology_df <- pharmacology_df[pharmacology_df$drugbank_id %in% drugbank_id_of_interest, ] 
subset_pharmacology_df_harmony <- pharmacology_df[pharmacology_df$drugbank_id %in% drugbank_id_of_interest_harmony, ] 

subset_atc_codes_df <- atc_codes_df[atc_codes_df$`drugbank-id` %in% drugbank_id_of_interest, ] 
subset_atc_codes_df_harmony <- atc_codes_df[atc_codes_df$`drugbank-id` %in% drugbank_id_of_interest_harmony, ] 

# Merge subset_drug_df, subset_pharmacology_df and subset_atc_codes_df
merged_df <- merge(subset_drug_df, subset_atc_codes_df, by.x = "primary_key", by.y = "drugbank-id")
merged_df <- merge(merged_df, subset_pharmacology_df, by.x = "primary_key", by.y = "drugbank_id")

merged_df_harmony <- merge(subset_drug_df_harmony, subset_atc_codes_df_harmony, by.x = "primary_key", by.y = "drugbank-id")
merged_df_harmony <- merge(merged_df_harmony, subset_pharmacology_df_harmony, by.x = "primary_key", by.y = "drugbank_id")

# Retain only columns of interest in merged_df
columns_of_interest <- c("primary_key", "name", "description", "monoisotopic_mass", "level_1", "level_2", "level_3", "level_4", "indication", "pharmacodynamics", "mechanism_of_action", "toxicity", "metabolism", "absorption", "half_life", "protein_binding", "route_of_elimination", "volume_of_distribution")
merged_df <- merged_df[, columns_of_interest]
merged_df_harmony <- merged_df_harmony[, columns_of_interest]


# Save dataframes to file
write.csv(merged_df, "./drugbank_infos.csv", row.names = FALSE)
write.csv(merged_df_harmony, "./drugbank_infos_harn´mony.csv", row.names = FALSE)

################################################ Visualize results ##################################################################

# Create a barplot of the frequency of highest "level"
ggplot(data = merged_df, aes(x = level_4)) +
  geom_bar(fill = "blue") +
  labs(x = "Level 4", y = "Frequency") +
  ggtitle("Before harmony integration: Misclassified treatments level 4 category") +
  coord_flip()

# Create a barplot of the frequency of highest "level"
ggplot(data = merged_df_harmony, aes(x = level_4)) +
  geom_bar(fill = "blue") +
  labs(x = "Level 4", y = "Frequency") +
  ggtitle("After harmony integration: Misclassified treatments level 4 category") +
  coord_flip()

