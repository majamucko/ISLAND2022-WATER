####Beta diversity comparisons bacteria vs. eukaryotes from ISLAND2022 WATER samples
##Maja Mucko
######Packages
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(fantaxtic)
library(factoextra)
library(Matrix)
library(reshape)
library(reshape2)
library(mdatools)
library(pairwiseAdonis)
library(magrittr)
library(dplyr)
library(microbiome)
library(file2meco)
library(vegan)
#install.packages("gplots")
library(gplots)
library(dendextend)
#devtools::install_github("zdk123/SpiecEasi", force = TRUE)
library(SpiecEasi)
#install.packages("BiocManager")
#devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports", "LinkingTo"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
library(NetCoMi)
#remotes::install_github("microbiome/OMA", dependencies = TRUE, upgrade = TRUE)
library(mia)


###Load 16S and 18S filtered OTU tables
Tang_x<-read.table("data/Filtered_B_ASV_table_7509.tsv") #Bacteria
Tang_y<-read.table("data/Filtered_E_ASV_table_5947.tsv") #Eukaryotes
# Load metadata
metadata <- read.table("data/metadata_water_ISLAND2022_envs.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Ensure row names match your OTU table sample IDs
all(rownames(Tang_x) %in% rownames(metadata)) # This should return TRUE
all(rownames(Tang_y) %in% rownames(metadata)) # This should return TRUE

###################
#####Tanglegram All
library(dendextend)

Tang_xd<-vegdist(Tang_x,method="bray")
Tang_xx<-as.dendrogram(hclust(Tang_xd,"average"))
Tang_yd<-vegdist(Tang_y,method="bray")
Tang_yy<-as.dendrogram(hclust(Tang_yd,"average"))

# Define custom colors for 'fraction'
fraction_colors <- c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F")
metadata$Color <- fraction_colors[metadata$fraction] # Map colors to samples

# Function to color labels based on metadata
color_labels <- function(dend, metadata) {
  labels <- labels(dend) # Get labels
  label_colors <- metadata[labels, "Color"] # Map labels to metadata colors
  labels_colors(dend) <- label_colors # Apply colors
  return(dend)
}

# Apply colors to both dendrograms
Tang_xx_colored <- color_labels(Tang_xx, metadata)
Tang_yy_colored <- color_labels(Tang_yy, metadata)

setdiff(labels(Tang_xx), rownames(metadata))
# Define symbols for 'ITW' categories
ITW_symbols <- c("yes" = 19, "no" = 1)  # Filled circle for 'yes', empty circle for 'no'
metadata$Symbol <- ITW_symbols[metadata$ITW]

# Function to add symbols to leaves based on metadata
add_symbols <- function(dend, metadata) {
  labels <- labels(dend)  # Get labels from the dendrogram
  label_symbols <- metadata[labels, "Symbol"]  # Map labels to metadata symbols
  
  # Assign symbols to leaves
  dend <- set(dend, "leaves_pch", label_symbols)
  return(dend)
}
# Apply symbols to both dendrograms
Tang_xx_symbol <- add_symbols(Tang_xx_colored, metadata)
Tang_yy_symbol <- add_symbols(Tang_yy_colored, metadata)

# Add titles during the tanglegram plotting
tanglegram(
  rank_branches(Tang_yy_symbol),
  rank_branches(Tang_xx_symbol),
  lab.cex = 0.9, edge.lwd = 1, lwd = 1,
  margin_inner = 5, center = TRUE,
  common_subtrees_color_lines = FALSE,
  highlight_distinct_edges = FALSE,
  columns_width = c(5, 3, 5), sort = TRUE,
  main_left = "Eukaryotes",  # Add title for left dendrogram
  main_right = "Bacteria"   # Add title for right dendrogram
)

###################
#####Mantel test 
Tang_Axd<-vegdist(Tang_xd,"bray")
Tang_Ayd<-vegdist(Tang_yd,"bray")
mantel_result <- mantel(Tang_xd,Tang_yd,permutations=999)

# Extract specific components
# Extract components
mantel_stat <- mantel_result$statistic
mantel_significance <- mantel_result$signif

# Calculate quantiles from permutation distribution
mantel_quantiles <- quantile(mantel_result$perm, probs = c(0.9, 0.95, 0.975, 0.99))

# Combine into a data frame for saving
mantel_details <- data.frame(
  Mantel_Statistic = mantel_stat,
  Significance = mantel_significance,
  `90% Quantile` = mantel_quantiles[1],
  `95% Quantile` = mantel_quantiles[2],
  `97.5% Quantile` = mantel_quantiles[3],
  `99% Quantile` = mantel_quantiles[4]
)

# Save to a CSV file
write.csv(mantel_details, file = "statistics/18S_16S_mantel_test_details.csv", row.names = FALSE)

#################
####Heatmap with agglomerated class level top 10 class bacteria and eukaryotes
# Load necessary libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(phyloseq)
B_ps_rare <- readRDS("data/phyloseq_16S_rare.rds")
E_ps_rare <- readRDS("data/phyloseq_18S_rare.rds")

# Step 1: Agglomerate to Class level for Bacteria
B_ps_class <- tax_glom(B_ps_rare, taxrank = "Class")

# Step 2: Agglomerate to Class level for Eukaryotes
E_ps_class <- tax_glom(E_ps_rare, taxrank = "Class")

# Step 3: Convert Bacteria Phyloseq Object to Data Frame
B_ps_df <- psmelt(B_ps_class) %>%
  group_by(Sample, Class) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))  # Calculate relative abundance

# Step 4: Convert Eukaryote Phyloseq Object to Data Frame
E_ps_df <- psmelt(E_ps_class) %>%
  group_by(Sample, Class) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))  # Calculate relative abundance

# Step 5: Get the top 10 Classes based on Relative Abundance for Bacteria
top_B_classes <- B_ps_df %>%
  group_by(Class) %>%
  summarise(Total_Relative_Abundance = sum(Relative_Abundance)) %>%
  arrange(desc(Total_Relative_Abundance)) %>%
  head(10)

# Step 6: Get the top 10 Classes based on Relative Abundance for Eukaryotes
top_E_classes <- E_ps_df %>%
  group_by(Class) %>%
  summarise(Total_Relative_Abundance = sum(Relative_Abundance)) %>%
  arrange(desc(Total_Relative_Abundance)) %>%
  head(10)

# Step 7: Filter original dataframes to include only top 10 classes for both
B_top_classes_df <- B_ps_df %>% filter(Class %in% top_B_classes$Class)
head(B_top_classes_df)
E_top_classes_df <- E_ps_df %>% filter(Class %in% top_E_classes$Class)
head(E_top_classes_df)

# Reshape the bacterial dataframe (B_top_classes_df)
B_top_classes_wide <- B_top_classes_df %>%
  pivot_wider(names_from = Sample, values_from = c(Relative_Abundance)) %>%
  arrange(Class) %>%
  rename_with(~paste("Bacterial", ., sep = "_"))  # Prefix for bacterial data
head(B_top_classes_wide)
# Load necessary library
library(dplyr)

# Group by 'Bacterial_Class' and summarize the abundance data
B_top_classes_summarized <- B_top_classes_wide %>%
  group_by(Bacterial_Class) %>%
  summarise(across(starts_with("Bacterial_W"), sum, na.rm = TRUE)) %>%
  ungroup()  # Ungroup after summarizing

# View the summarized data
head(B_top_classes_summarized)
# Rename columns by removing the 'Bacterial_' prefix
colnames(B_top_classes_summarized) <- gsub("Bacterial_", "", colnames(B_top_classes_summarized))

# Optionally, rename 'Bacterial_Class' to 'Class'
colnames(B_top_classes_summarized)[1] <- "Class"

# View the updated column names
colnames(B_top_classes_summarized)
# Export Bacterial data to .tsv
write.table(B_top_classes_summarized, file = "data/B_top_classes_summarized.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# Reshape the eukaryotic dataframe (E_top_classes_df)
E_top_classes_wide <- E_top_classes_df %>%
  pivot_wider(names_from = Sample, values_from = c(Relative_Abundance)) %>%
  arrange(Class) %>%
  rename_with(~paste("Eukaryotic", ., sep = "_"))  # Prefix for eukaryotic data
head(E_top_classes_wide)
# Group by 'Bacterial_Class' and summarize the abundance data
E_top_classes_summarized <- E_top_classes_wide %>%
  group_by(Eukaryotic_Class) %>%
  summarise(across(starts_with("Eukaryotic_W"), sum, na.rm = TRUE)) %>%
  ungroup()  # Ungroup after summarizing

# View the summarized data
head(E_top_classes_summarized)
# Rename columns by removing the 'Bacterial_' prefix
colnames(E_top_classes_summarized) <- gsub("Eukaryotic_", "", colnames(B_top_classes_summarized))

# Optionally, rename 'Bacterial_Class' to 'Class'
colnames(E_top_classes_summarized)[1] <- "Class"

# View the updated column names
colnames(E_top_classes_summarized)
# Export Eukaryotic data to .tsv
write.table(E_top_classes_summarized, file = "data/E_top_classes_summarized.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
##############I joined df by hand in excel and removed c__ prefix from Bacterial classes

# Load the merged data
heatmap_B_E_top_class <- read.table("data/merged_B_E_top_classes.txt", header = TRUE)
head(heatmap_B_E_top_class)
print(heatmap_B_E_top_class)
library(pheatmap)
# Calculate Bray-Curtis distance between samples
library(vegan)
numeric_data <- heatmap_B_E_top_class[, -1]  # Remove the first column (Class)
numeric_data <- apply(numeric_data, 2, as.numeric)
rownames(numeric_data) <- heatmap_B_E_top_class$Class  # Set class names as rownames
dist_matrix <- vegdist(t(numeric_data), method = "bray")
sum(is.na(heatmap_B_E_top_class))

# Plot heatmap with title indicating Bray-Curtis clustering
pheatmap(
  numeric_data,  # Use numeric-only data
  clustering_distance_cols = dist_matrix,  # Cluster columns based on Bray-Curtis distance
  clustering_method = "average",  # Use average linkage for clustering
  cluster_rows = FALSE,  # Disable row clustering
  color = colorRampPalette(c("skyblue1", "blue", "darkblue", "thistle1", "red"))(50),
  show_rownames = TRUE,  # Show class names
  show_colnames = TRUE,  # Show sample names
  fontsize = 12,
  scale = "none",  # Optional: scale data if needed
  main = "Heatmap with Clustering Based on Bray-Curtis Distance"
)
#### Heatmap with agglomerated phylum level top 10 phyla for bacteria and eukaryotes
# Load necessary libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(vegan)
library(pheatmap)

# Load phyloseq objects
B_ps_rare <- readRDS("data/phyloseq_16S_rare.rds")
E_ps_rare <- readRDS("data/phyloseq_18S_rare.rds")

# Step 1: Agglomerate to Phylum level for Bacteria
B_ps_phylum <- tax_glom(B_ps_rare, taxrank = "Phylum")

# Step 2: Agglomerate to Subdivision level for Eukaryotes (equivalent to Phylum)
E_ps_phylum <- tax_glom(E_ps_rare, taxrank = "Subdivision")

# Step 3: Convert Bacteria Phyloseq Object to Data Frame
B_ps_df <- psmelt(B_ps_phylum) %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))  # Calculate relative abundance

# Step 4: Convert Eukaryote Phyloseq Object to Data Frame
E_ps_df <- psmelt(E_ps_phylum) %>%
  group_by(Sample, Subdivision) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))  # Calculate relative abundance

# Rename "Subdivision" to "Phylum" for consistency
head(E_ps_df)
colnames(E_ps_df)
colnames(E_ps_df) <- trimws(colnames(E_ps_df))  # Remove extra spaces
library(dplyr)
colnames(E_ps_df)[colnames(E_ps_df) == "Subdivision"] <- "Phylum"

# Step 5: Get the top 10 Phyla based on Relative Abundance for Bacteria
top_B_phyla <- B_ps_df %>%
  group_by(Phylum) %>%
  summarise(Total_Relative_Abundance = sum(Relative_Abundance)) %>%
  arrange(desc(Total_Relative_Abundance)) %>%
  head(10)

# Step 6: Get the top 10 Phyla based on Relative Abundance for Eukaryotes
top_E_phyla <- E_ps_df %>%
  group_by(Phylum) %>%
  summarise(Total_Relative_Abundance = sum(Relative_Abundance)) %>%
  arrange(desc(Total_Relative_Abundance)) %>%
  head(10)

# Step 7: Filter original dataframes to include only top 10 phyla for both
B_top_phyla_df <- B_ps_df %>% filter(Phylum %in% top_B_phyla$Phylum)
E_top_phyla_df <- E_ps_df %>% filter(Phylum %in% top_E_phyla$Phylum)

# Step 8: Reshape the bacterial dataframe (B_top_phyla_df)
B_top_phyla_wide <- B_top_phyla_df %>%
  pivot_wider(names_from = Sample, values_from = c(Relative_Abundance)) %>%
  arrange(Phylum) %>%
  rename_with(~paste("Bacterial", ., sep = "_"))  # Prefix for bacterial data
head(B_top_phyla_wide)

# Step 9: Group by Phylum and summarize the abundance data for Bacteria
B_top_phyla_summarized <- B_top_phyla_wide %>%
  group_by(Bacterial_Phylum) %>%
  summarise(across(starts_with("Bacterial_W"), sum, na.rm = TRUE)) %>%
  ungroup()  # Ungroup after summarizing
# View the summarized data
head(B_top_phyla_summarized)
# Rename columns by removing the 'Bacterial_' prefix
colnames(B_top_phyla_summarized) <- gsub("Bacterial_", "", colnames(B_top_phyla_summarized))
# Optionally, rename 'Bacterial_Phylum' to 'Phylum'
colnames(B_top_phyla_summarized)[1] <- "Phylum"
# View the updated column names
colnames(B_top_phyla_summarized)
# Export Bacterial data to .tsv
write.table(B_top_phyla_summarized, file = "data/B_top_phyla_summarized.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Step 10: Reshape the eukaryotic dataframe (E_top_phyla_df)
E_top_phyla_wide <- E_top_phyla_df %>%
  pivot_wider(names_from = Sample, values_from = c(Relative_Abundance)) %>%
  arrange(Phylum) %>%
  rename_with(~paste("Eukaryotic", ., sep = "_"))  # Prefix for eukaryotic data
head(E_top_phyla_wide)

# Step 11: Group by Phylum and summarize the abundance data for Eukaryotes
E_top_phyla_summarized <- E_top_phyla_wide %>%
  group_by(Eukaryotic_Phylum) %>%
  summarise(across(starts_with("Eukaryotic_W"), sum, na.rm = TRUE)) %>%
  ungroup()  # Ungroup after summarizing
# View the summarized data
head(E_top_phyla_summarized)

# Rename columns by removing the 'Eukaryotic_' prefix
colnames(E_top_phyla_summarized) <- gsub("Eukaryotic_", "", colnames(E_top_phyla_summarized))
# Optionally, rename 'Eukaryotic_Phylum' to 'Phylum'
colnames(E_top_phyla_summarized)[1] <- "Phylum"
# View the updated column names
colnames(E_top_phyla_summarized)
# Export Eukaryotic data to .tsv
write.table(E_top_phyla_summarized, file = "data/E_top_phyla_summarized.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
##############I joined df by hand in excel and removed c__ prefix from Bacterial classes
# Load the merged data
heatmap_B_E_top_phylum <- read.table("data/merged_B_E_top_phylum.txt", header = TRUE)
head(heatmap_B_E_top_phylum)
library(pheatmap)
# Calculate Bray-Curtis distance between samples
library(vegan)
numeric_data <- heatmap_B_E_top_phylum[, -1]  # Remove the first column (Class)
numeric_data <- apply(numeric_data, 2, as.numeric)
rownames(numeric_data) <- heatmap_B_E_top_phylum$Phylum  # Set class names as rownames
dist_matrix <- vegdist(t(numeric_data), method = "bray")  # Transpose to get distance for columns
sum(is.na(heatmap_B_E_top_phylum))
# Check dimensions of the data
print(dim(numeric_data))  # Should be (20, 24) for 20 Phylum and 24 samples

# Plot heatmap with title indicating Bray-Curtis clustering
pheatmap(
  numeric_data,  # Use numeric-only data
  clustering_distance_cols = dist_matrix,  # Cluster columns based on Bray-Curtis distance
  clustering_method = "average",  # Use average linkage for clustering
  cluster_rows = FALSE,  # Disable row clustering
  color = colorRampPalette(c("skyblue1", "blue", "darkblue", "thistle1", "red"))(50),
  show_rownames = TRUE,  # Show class names
  show_colnames = TRUE,  # Show sample names
  fontsize = 12,
  scale = "none",  # Optional: scale data if needed
  main = "Phylum Heatmap with Clustering Based on Bray-Curtis Distance"
)
