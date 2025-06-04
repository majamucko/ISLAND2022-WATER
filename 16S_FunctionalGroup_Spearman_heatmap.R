#####Trying to test functional profiling with environmentals - bacterial dataset

# Load necessary libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)

# Load data
otu_table <- read.table("data/functionalGroup_otu_table.tsv", sep = "\t", header = TRUE, check.names = FALSE)
head(otu_table)
metadata <- read.table("data/metadata_water_ISLAND2022_envs.tsv", sep = "\t", header = TRUE, row.names = 1)
head(metadata)

# Transpose the OTU table (samples as rows, functions as columns)
otu_table_t <- t(otu_table[,-1])  # Exclude "FunctionalGroup" column
colnames(otu_table_t) <- otu_table$FunctionalGroup
otu_table_t <- as.data.frame(otu_table_t)
otu_table_t$Sample <- rownames(otu_table_t)

# Merge OTU table with metadata
merged_data <- merge(otu_table_t, metadata, by.x = "Sample", by.y = "row.names")
print(merged_data)
# Check the first few rows of the OTU table and the metadata
head(otu_table_t)
head(metadata)
# Check the column names of both dataframes
colnames(otu_table_t)
colnames(metadata)
# Check if row names of both dataframes match
all(rownames(otu_table_t) %in% rownames(metadata))

# Perform Spearman correlation for each functional group against environmental variables
results_spearman <- data.frame(FunctionalGroup = character(), Variable = character(), Correlation = numeric(), P_Value = numeric())

env_vars <- c("Temperature", "Salinity", "NO2", "NO3", "PO4", "Chl_a", "SiO4", "Chl_F")  # Adjust for your variables

for (func in colnames(otu_table_t)[-ncol(otu_table_t)]) {
  for (var in env_vars) {
    corr_test <- cor.test(merged_data[[func]], merged_data[[var]], method = "spearman", use = "complete.obs")
    results_spearman <- rbind(results_spearman, data.frame(
      FunctionalGroup = func,
      Variable = var,
      Correlation = corr_test$estimate,
      P_Value = corr_test$p.value
    ))
  }
}

# Adjust p-values for multiple comparisons
results_spearman$Adjusted_P <- p.adjust(results_spearman$P_Value, method = "fdr")

# View significant correlations
significant_spearman <- results_spearman %>% filter(Adjusted_P < 0.05)
print(significant_spearman)

# Reshape Spearman results for heatmap
library(pheatmap)

heatmap_data <- significant_spearman %>%
  pivot_wider(names_from = Variable, values_from = Correlation, values_fill = 0)

# Create heatmap
pheatmap(as.matrix(heatmap_data[,-1]), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Significant Spearman Correlations")

# Filter significant results for bubble plot
bubble_data <- otu_table %>%
  filter(FunctionalGroup %in% significant_spearman$FunctionalGroup)

ggplot(bubble_data_long, aes(x = Sample, y = FunctionalGroup, size = Percentage, fill = Percentage)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_gradient(low = "blue", high = "red") + # Color gradient
  scale_size_area(max_size = 20) +  # Increase bubble size by adjusting max_size
  ylab("Functional Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Load necessary libraries
library(tidyverse)
library(reshape2)
library(pheatmap)

# Load and preprocess data (assuming data is already loaded)
# Transpose OTU table
otu_table_t <- t(otu_table[,-1])
colnames(otu_table_t) <- otu_table$FunctionalGroup
otu_table_t <- as.data.frame(otu_table_t)
otu_table_t$Sample <- rownames(otu_table_t)

# Merge OTU table with metadata
merged_data <- merge(otu_table_t, metadata, by.x = "Sample", by.y = "row.names")

# Perform Spearman correlation
results_spearman <- data.frame(FunctionalGroup = character(), Variable = character(), Correlation = numeric(), P_Value = numeric())

env_vars <- c("Temperature", "Salinity", "NO2","NO3", "PO4", "SiO4", "Chl_a", "Chl_F")  # Adjust for your variables

for (func in colnames(otu_table_t)[-ncol(otu_table_t)]) {
  for (var in env_vars) {
    corr_test <- cor.test(merged_data[[func]], merged_data[[var]], method = "spearman", use = "complete.obs")
    results_spearman <- rbind(results_spearman, data.frame(
      FunctionalGroup = func,
      Variable = var,
      Correlation = corr_test$estimate,
      P_Value = corr_test$p.value
    ))
  }
}

# Adjust p-values for multiple testing
results_spearman$Adjusted_P <- p.adjust(results_spearman$P_Value, method = "fdr")

# Identify significant correlations
significant_spearman <- results_spearman %>% filter(Adjusted_P < 0.05)

# Reshape OTU data for heatmap
heatmap_data <- otu_table_t %>%
  pivot_longer(-Sample, names_to = "FunctionalGroup", values_to = "Percentage") %>%
  pivot_wider(names_from = Sample, values_from = Percentage, values_fill = 1)
print(heatmap_data)

# Ensure heatmap_matrix has row names corresponding to FunctionalGroups
heatmap_matrix <- heatmap_data %>%
  column_to_rownames("FunctionalGroup") %>%
  as.matrix()

# Create annotation matrix for significance (same dimensions as heatmap_matrix)
annotation_matrix <- matrix(
  "",  # Default: no annotation
  nrow = nrow(heatmap_matrix),
  ncol = ncol(heatmap_matrix),
  dimnames = list(rownames(heatmap_matrix), colnames(heatmap_matrix))
)

# Mark significant FunctionalGroups with stars
significant_groups <- significant_spearman$FunctionalGroup
annotation_matrix[rownames(heatmap_matrix) %in% significant_groups, ] <- "*"

# Plot the heatmap
pheatmap(
  mat = heatmap_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "Functional Group Percentages with Significant Correlations",
  display_numbers = annotation_matrix,  # Overlay annotation matrix
  number_color = "black",
  fontsize_number = 8,
  labels_row = rownames(heatmap_matrix),  # Show FunctionalGroups on rows
  color = colorRampPalette(c("white", "blue"))(100)
)
# Filter FunctionalGroups to include significant ones or those with max percentage >= 0.005, excluding "Other"
filtered_heatmap_data <- heatmap_data %>%
  rowwise() %>%
  mutate(MaxPercentage = max(c_across(-FunctionalGroup))) %>%  # Compute max percentage per group
  filter(
    (FunctionalGroup %in% significant_spearman$FunctionalGroup | MaxPercentage >= 0.005) & 
      FunctionalGroup != "other"  # Exclude "Other"
  ) %>%
  select(-MaxPercentage)  # Drop temporary column

# View the filtered heatmap data
print(filtered_heatmap_data)

# Ensure heatmap_matrix has row names corresponding to filtered FunctionalGroups
heatmap_matrix <- filtered_heatmap_data %>%
  column_to_rownames("FunctionalGroup") %>%
  as.matrix()

# Create annotation matrix for significance (same dimensions as heatmap_matrix)
annotation_matrix <- matrix(
  "",  # Default: no annotation
  nrow = nrow(heatmap_matrix),
  ncol = ncol(heatmap_matrix),
  dimnames = list(rownames(heatmap_matrix), colnames(heatmap_matrix))
)

# Mark significant FunctionalGroups with stars
significant_groups <- significant_spearman$FunctionalGroup
annotation_matrix[rownames(heatmap_matrix) %in% significant_groups, ] <- "*"

# Plot the heatmap
pheatmap(
  mat = heatmap_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "Functional Group Percentages with Significant Correlations",
  display_numbers = annotation_matrix,  # Overlay annotation matrix
  number_color = "black",
  fontsize_number = 10,
  labels_row = rownames(heatmap_matrix),  # Show FunctionalGroups on rows
  color = colorRampPalette(c("white","blue", "red"))(100)
)

# Filter for significant functional groups and environmental variables
significant_pairs <- significant_spearman %>%
  select(FunctionalGroup, Variable) %>%
  unique()
# Initialize a table for storing results
sample_significance <- data.frame(
  FunctionalGroup = character(),
  Variable = character(),
  Sample = character(),
  Value = numeric(),
  Environmental_Value = numeric()
)

# Loop through each significant FunctionalGroup and Variable
for (i in 1:nrow(significant_pairs)) {
  func <- significant_pairs$FunctionalGroup[i]
  variable <- significant_pairs$Variable[i]
  
  # Check each sample
  for (sample in rownames(merged_data)) {
    func_value <- merged_data[sample, func]
    env_value <- merged_data[sample, variable]
    
    if (!is.na(func_value) && !is.na(env_value)) {
      # Add to the table if environmental value is above/below threshold
      sample_significance <- rbind(sample_significance, data.frame(
        FunctionalGroup = func,
        Variable = variable,
        Sample = sample,
        Value = func_value,
        Environmental_Value = env_value
      ))
    }
  }
}

# View the resulting table
print(sample_significance)

# Optionally save the table
write.table(sample_significance, "sample_significance_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Define sample groups
sample_groups <- data.frame(
  Sample = colnames(heatmap_matrix),
  Group = ifelse(colnames(heatmap_matrix) %in% c("W1", "W2", "W3", "W4", "W5", "W6", 
                                                 "W7", "W8", "W9", "W10", "W11", "W12", 
                                                 "W13", "W14", "W15", "W22", "W23", "W24"),
                 "ITW", "noITW")  # ITW and noITW groups
)

# Set row names of sample_groups to match samples
rownames(sample_groups) <- sample_groups$Sample
sample_groups <- sample_groups %>% select(Group)  # Keep only Group column

# Validate sample_order
sample_order <- c(
  intersect(colnames(heatmap_matrix), c("W1", "W2", "W3", "W4", "W5", "W6",
                                        "W7", "W8", "W9", "W10", "W11", "W12",
                                        "W13", "W14", "W15", "W22", "W23", "W24")),
  intersect(colnames(heatmap_matrix), c("W16", "W17", "W18", "W19", "W20", "W21"))
)

# Debug: Check sample_order and heatmap_matrix columns
if (length(sample_order) != ncol(heatmap_matrix)) {
  stop("Mismatch between sample_order and heatmap_matrix columns.")
}

# Reorder heatmap_matrix based on sample_order
heatmap_matrix <- heatmap_matrix[, sample_order]

# Reorder annotation_matrix based on sample_order
annotation_matrix <- annotation_matrix[, sample_order]

# Replace NA/NaN/Inf in heatmap_matrix
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Plot heatmap with pre-defined clustering
pheatmap(
  mat = heatmap_matrix, 
  cluster_rows = TRUE,   # Cluster rows
  cluster_cols = FALSE,  # Disable automatic clustering for columns
  annotation_col = sample_groups,  # Add the sample grouping annotation
  display_numbers = annotation_matrix,  # Use significance stars
  number_color = "black",
  fontsize_number = 12,
  main = "Functional Group Percentages with ITW/NoITW Grouping",
  color = colorRampPalette(c("white", "blue", "red"))(100)
)

# Define custom colors for annotation_col
annotation_colors <- list(
  Group = c(ITW = "black", noITW = "grey")  # Customize colors for ITW and noITW
)
# Plot heatmap with annotation_col colors
pheatmap(
  mat = heatmap_matrix, 
  cluster_rows = TRUE,   # Cluster rows
  cluster_cols = FALSE,  # Disable automatic clustering for columns
  annotation_col = sample_groups,  # Add the sample grouping annotation
  annotation_colors = annotation_colors,  # Use custom colors for annotation
  display_numbers = annotation_matrix,  # Use significance stars
  number_color = "black",
  fontsize_number = 12,
  main = "Functional Group Percentages with ITW/NoITW Grouping",
  color = colorRampPalette(c("skyblue1", "blue", "darkblue", "thistle1", "red")) (100)
)
