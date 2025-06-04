## Code for processing 18S rRNA gene metabarcoding data (V4) from ISLAND-WATER072024
## Import of raw .qza files into phyloseq object and basic alpha / beta diversity index calculations
## Statistical analysis of indices
## Visualizations of jitter plots and nMDS
# Maja Mucko

# Load packages
library(tidyverse)
#install.packages("vegan")
library(vegan)
#install.packages("devtools")
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(devtools)
#devtools::install_github("gmteunisse/fantaxtic")
library(fantaxtic)
library(RColorBrewer)
#devtools::install_github("microbiome/microbiome")
library(microbiome)
#install.packages("factoextra")
library(factoextra)
library(Matrix)
#install.packages("microeco")
library(microeco)
#install.packages("file2meco")
library(file2meco)
#install.packages("mdatools")
library(mdatools)
#install.packages("ComplexUpset")
library(ComplexUpset)
#install.packages("reshape")
library(reshape)
library(reshape2)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(magrittr)
library(dplyr)
library(tibble)
library(tidyverse)

# Load 18S taxonomy file - assigned with new PR2 database release (V 5.0.1)
taxonomy <- read_qza(file="data/WATER_18s_taxonomy_pr2.qza")
tax_tab <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)

# Load 18S metadata file from BATS
metadata <- read.table("data/metadata_water_ISLAND2022_envs.tsv", header=TRUE, row.names=1, check.names=F, sep="\t") 

# Load rooted tree
rooted_tree <- read_qza(file = "data/18S_rooted_tree.qza")$data

# Load ASV count table 
table <- read_qza(file = "data/WATER_18s_table_trunc220.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 

# Merge into phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax_tab)), otu_table(count_tab, taxa_are_rows = T), sample_data(metadata), phy_tree(rooted_tree))

# Remove unwanted groups
ps_new = subset_taxa(ps, Division !="Streptophyta" |is.na(Division))
ps_new = subset_taxa(ps, Supergroup !="Eukaryota:nucl" |is.na(Supergroup))
ps_new = subset_taxa(ps, Division !="Eukaryota:nucl" |is.na(Division))
ps_new = subset_taxa(ps_new, Division !="Rhodophyta" |is.na(Division))
ps_new = subset_taxa(ps_new, Division !="Cryptophyta:nucl" |is.na(Division))
ps_new = subset_taxa(ps_new, Domain !="Bacteria" |is.na(Domain))
ps_new <- subset_taxa(ps_new, Subdivision!="Opisthokonta_X"|is.na(Subdivision))
ps_new <- subset_taxa(ps_new, Class!="Unassigned", Prune = T)
ps_new = subset_taxa(ps_new, Class!="Craniata"|is.na(Class))
ps_new <- subset_taxa(ps_new, Subdivision != "Unassigned" | is.na(Subdivision))
ps_new <- subset_taxa(ps_new, Class != "Chloroplast" | is.na(Class))
ps_new <- subset_taxa(ps_new, Order != "Mitochondria" | is.na(Order))
ps_new <- subset_taxa(ps_new, Class != "Embryophyceae" | is.na(Class))
ps_new <- subset_taxa(ps_new, Subdividion != "Ichthyophonus" | is.na(Subdivision))
ps_new <- subset_taxa(ps_new, Class != "Phaeophyceae" | is.na(Class))
ps_new <- subset_taxa(ps_new, Class != "Porifera" | is.na(Class))
ps_new <- subset_taxa(ps_new, Class != "Cryptophyceae:nucl" | is.na(Class))

ps_new = name_na_taxa(ps_new) # Adds an unassigned label to better identify lowest possible taxonomic assignment
#agglomerate to genus level and export that table - to have all singletons, just to annotate morphological data
taxa_names(ps_new) <- paste0("ASV", seq(ntaxa(ps_new)))
# Agglomerate ASVs to Genus level
#ps_genus <- tax_glom(ps_new, "Genus", NArm = TRUE)

# Extract the OTU table (counts) and taxonomic data
#otu_tab_genus <- as.data.frame(otu_table(ps_genus))
#tax_tab_genus <- as.data.frame(tax_table(ps_genus))

# Merge taxonomic info with count table
genus_table <- cbind(tax_tab_genus, otu_tab_genus)

# Export the table to CSV
#write.csv(genus_table, "genus_level_table.csv", row.names = TRUE)


# Remove singletons (ASVs present only once)
ps_filt = filter_taxa(ps_new, function (x) {sum(x) > 1}, prune=TRUE)
# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps_filt, sample.size = min(sample_sums(ps_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
# Rename ASVs in sequential order
taxa_names(ps_rare) <- paste0("ASV", seq(ntaxa(ps_rare)))
saveRDS(ps_rare, file = "data/phyloseq_18S_rare.rds")
asv_table_filtered_rarefied <- as.data.frame(t(otu_table(ps_rare)))
write.table(asv_table_filtered_rarefied, file = "data/Filtered_E_ASV_table_5947.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Rarefaction curves for 18S - Figure S1
ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

rare_18S <- ggrare(ps_filt, step = 100, plot = FALSE, parallel = FALSE, se = FALSE)
rare_18S + theme(legend.position = "none") + theme_q2r()+ theme(legend.position = "right") + 
  facet_wrap(~fraction, scales="free_y",nrow=4,ncol=3)
ggsave("plots/18S/18S_rarefaction_fractions_6501ASV.pdf", plot = last_plot(), height = 4, width = 6, device = "pdf")

# Estimate minimum, mean, and maximum 18S read counts
ps_min <- min(sample_sums(ps_rare))
ps_mean <- mean(sample_sums(ps_rare))
ps_max <- max(sample_sums(ps_rare))
#Create a data frame
summary_table <- data.frame(
  Metric = c("Min", "Mean", "Max"),
  Value = c(ps_min, ps_mean, ps_max)
)
# Export to CSV
write.csv(summary_table, file = "statistics/18S/summary_table_5947ASV.csv", row.names = FALSE)
# Export 18S ASV information - filtered, rarefied, and relabeled ASVs (used in data analysis for this study) - Table S4
OTU_filt = as(otu_table(ps_rare), "matrix")
TAX_filt = as(tax_table(ps_rare), "matrix")
merge_18S_filt <- cbind(OTU_filt,TAX_filt)
write.csv(merge_18S_filt, file="data/18SASVs_for_trophy_annotation.csv", row.names=T)

#######################################ALPHA DIVERSITY###########################################################
# Calculate alpha diversity metrics
alpha_diversity <- estimate_richness(ps_rare, measures = c("Shannon", "Simpson", "Observed"))

# Calculate Evenness as Shannon diversity divided by the log of Observed richness
alpha_diversity$Evenness <- alpha_diversity$Shannon / log(alpha_diversity$Observed)

# View the updated alpha diversity table
head(alpha_diversity)

# Check the sample data to see your grouping variables
sample_data(ps_rare)
# Create an empty summary table for pairwise comparisons
summary_pairwise_alpha_div <- data.frame()

# Function to run pairwise tests for a given alpha diversity metric and grouping variable
run_pairwise_test <- function(metric, grouping_variable) {
  # Perform pairwise comparisons using Wilcoxon rank sum test (or t-test if appropriate)
  pairwise_result <- pairwise.wilcox.test(
    alpha_diversity[[metric]], 
    sample_data(ps_rare)[[grouping_variable]], 
    p.adjust.method = "BH"
  )
  
  # Extract the p-values from the pairwise test results
  p_values <- pairwise_result$p.value
  # Convert the p-values to a tidy format
  p_values_tidy <- as.data.frame(as.table(p_values))
  # Filter out NA values
  p_values_tidy <- na.omit(p_values_tidy)
  # Add the metric and grouping variable to the results
  p_values_tidy$Metric <- metric
  p_values_tidy$GroupingVariable <- grouping_variable
  
  return(p_values_tidy)
}

# Loop through each alpha diversity metric and each grouping variable
for (metric in colnames(alpha_diversity)) {
  for (grouping_variable in c("fraction", "depth", "ITW", "heatwave")) {
    # Run the pairwise test and combine results
    pairwise_results <- run_pairwise_test(metric, grouping_variable)
    summary_pairwise_alpha_div <- rbind(summary_pairwise_alpha_div, pairwise_results)
  }
}

# Rename columns for clarity
colnames(summary_pairwise_alpha_div) <- c("Group1", "Group2", "p.value", "Metric", "GroupingVariable")

# View the summary of pairwise comparisons
head(summary_pairwise_alpha_div)
write.csv(summary_pairwise_alpha_div, "statistics/18S/summary_pairwise_Wilcox_alpha_diversity.csv", row.names = TRUE)

####NOW plot jitter plots for Shannon, Simpson, Evenness and Observed indices
library(tidyverse)

# Convert the alpha diversity table to a tidy format
alpha_diversity_tidy <- alpha_diversity %>%
  rownames_to_column("SampleID") %>%  # Convert row names to a column
  pivot_longer(cols = -SampleID,        # Pivot to longer format
               names_to = "Metric", 
               values_to = "Value")

# Check if SampleID overlaps with metadata
metadata_samples <- rownames(metadata)  # Get SampleIDs from metadata
common_samples <- intersect(metadata_samples, alpha_diversity_tidy$SampleID)

# Join alpha diversity metrics with metadata
metadata_1 <- metadata %>%
  rownames_to_column("SampleID") %>%
  filter(SampleID %in% common_samples) %>% # Filter metadata to include only common samples
  left_join(alpha_diversity_tidy, by = "SampleID") # Join with alpha diversity data

# View the first few rows of the merged data
head(metadata_1)

# Plot bar charts with jitter sample points for Shannon diversity
p1 <- metadata_1 %>%
  filter(Metric == "Shannon") %>%  # Filter for Shannon diversity
  ggplot(aes(x = fraction, y = Value, fill = fraction)) +
  stat_summary(geom = "bar", fun.data = mean_se, color = "black") + # Mean with standard error
  geom_jitter(shape = 21, width = 0.2, height = 0) + # Jitter points
  coord_cartesian(ylim = c(1, 6)) + # Adjust y-axis limits
  facet_grid(~ depth) + # Create panels for each depth
  xlab("Fraction") +
  ylab("Shannon Diversity") +
  theme_q2r() + # Change to a clean theme
  scale_fill_manual(values = c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F")) + # Custom colors
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Larger x-axis text
    axis.text.y = element_text(size = 10), # Larger y-axis text
    axis.title.x = element_text(size = 10), # Larger x-axis label
    axis.title.y = element_text(size = 10), # Larger y-axis label
    strip.text = element_text(size = 10) # Larger facet labels
  )
p1

# Save the plot as a PDF
ggsave("plots/18S/Shannon_by_fraction.pdf", height = 3, width = 4, device = "pdf") # Adjust width to 4 inches

# Plot for Simpson diversity index (p2)
p2 <- metadata_1 %>%
  filter(Metric == "Simpson") %>%  # Filter for Simpson diversity
  ggplot(aes(x = fraction, y = Value, fill = fraction)) +
  stat_summary(geom = "bar", fun.data = mean_se, color = "black") + # Mean with standard error
  geom_jitter(shape = 21, width = 0.2, height = 0) + # Jitter points
  coord_cartesian(ylim = c(0.3, 1)) + # Adjust y-axis limits (change based on your data)
  facet_grid(~ depth) + # Create panels for each depth
  xlab("Fraction") +
  ylab("Simpson Diversity") +
  theme_q2r() +
  scale_fill_manual(values = c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F")) + # Custom colors
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Larger x-axis text
    axis.text.y = element_text(size = 10), # Larger y-axis text
    axis.title.x = element_text(size = 10), # Larger x-axis label
    axis.title.y = element_text(size = 10), # Larger y-axis label
    strip.text = element_text(size = 10) # Larger facet labels
  )
p2
# Save p2 as a PDF
ggsave("plots/18S/Simpson_by_fraction.pdf", height = 3, width = 4, device = "pdf")

# Plot for Observed features index (p3)
p3 <- metadata_1 %>%
  filter(Metric == "Observed") %>%  # Filter for Observed features
  ggplot(aes(x = fraction, y = Value, fill = fraction)) +
  stat_summary(geom = "bar", fun.data = mean_se, color = "black") + # Mean with standard error
  geom_jitter(shape = 21, width = 0.2, height = 0) + # Jitter points
  facet_grid(~ depth) + # Create panels for each depth
  xlab("Fraction") +
  ylab("Observed Features") +
  theme_q2r() + # Clean theme
  scale_fill_manual(values = c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F")) + # Custom colors
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Larger x-axis text
    axis.text.y = element_text(size = 10), # Larger y-axis text
    axis.title.x = element_text(size = 10), # Larger x-axis label
    axis.title.y = element_text(size = 10), # Larger y-axis label
    strip.text = element_text(size = 10) # Larger facet labels
  )
p3
# Save p3 as a PDF
ggsave("plots/18S/Observed_by_fraction.pdf", height = 3, width = 4, device = "pdf")

# Plot for Evenness index (p4)
p4 <- metadata_1 %>%
  filter(Metric == "Evenness") %>%  # Filter for Evenness
  ggplot(aes(x = fraction, y = Value, fill = fraction)) +
  stat_summary(geom = "bar", fun.data = mean_se, color = "black") + # Mean with standard error
  geom_jitter(shape = 21, width = 0.2, height = 0) + # Jitter points
  facet_grid(~ depth) + # Create panels for each depth
  xlab("Fraction") +
  ylab("Evenness") +
  theme_q2r() + # Clean theme
  scale_fill_manual(values = c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F")) + # Custom colors
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Larger x-axis text
    axis.text.y = element_text(size = 10), # Larger y-axis text
    axis.title.x = element_text(size = 10), # Larger x-axis label
    axis.title.y = element_text(size = 10), # Larger y-axis label
    strip.text = element_text(size = 10) # Larger facet labels
  )
p4
# Save p4 as a PDF
ggsave("plots/18S/Evenness_by_fraction.pdf", height = 3, width = 4, device = "pdf")
###Put the plots in one image, Figure S2, panel Bacterial community
library(gridExtra)
grid_plot <- grid.arrange(p1, p2, p3, p4, ncol = 4)
ggsave("plots/18S/combined_plots_alpha_div.pdf", grid_plot, width = 12, height = 4)  # size of image with all 4 plots!

###################################BETA DIVERSITY#######################################################################
# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)

# Beta diversity indices
bray_curtis <- phyloseq::distance(ps_rare, method = "bray")
jaccard <- phyloseq::distance(ps_rare, method = "jaccard")
weighted_unifrac <- phyloseq::distance(ps_rare, method = "wunifrac")
unweighted_unifrac <- phyloseq::distance(ps_rare, method = "unifrac")

# Combine all distances into a list
beta_div_list <- list(
  Bray_Curtis = bray_curtis,
  Jaccard = jaccard,
  Weighted_UniFrac = weighted_unifrac,
  Unweighted_UniFrac = unweighted_unifrac
)
# Ensure grouping variables are factors
sample_data(ps_rare)$fraction <- as.factor(sample_data(ps_rare)$fraction)
sample_data(ps_rare)$depth <- as.factor(sample_data(ps_rare)$depth)
sample_data(ps_rare)$ITW <- as.factor(sample_data(ps_rare)$ITW)
sample_data(ps_rare)$heatwave <- as.factor(sample_data(ps_rare)$heatwave)

# Define a function to run PERMANOVA for a given distance matrix and grouping variable
run_permanova <- function(distance_matrix, grouping_variable) {
  # Create a data frame that contains the grouping variable and sample names
  sample_data_df <- data.frame(sample = rownames(sample_data(ps_rare)), 
                               group = grouping_variable)
  
  # Run PERMANOVA using the formula interface
  adonis_result <- adonis2(distance_matrix ~ group, data = sample_data_df)
  return(adonis_result)
}

# Run PERMANOVA for Bray-Curtis, Jaccard, WUniFrac and UniFrac distances using each grouping variable
permanova_results <- list(
  Bray_Curtis_Fraction = run_permanova(bray_curtis, sample_data(ps_rare)$fraction),
  Jaccard_Fraction = run_permanova(jaccard, sample_data(ps_rare)$fraction),
  
  Bray_Curtis_Depth = run_permanova(bray_curtis, sample_data(ps_rare)$depth),
  Jaccard_Depth = run_permanova(jaccard, sample_data(ps_rare)$depth),
  
  Bray_Curtis_ITW = run_permanova(bray_curtis, sample_data(ps_rare)$ITW),
  Jaccard_ITW = run_permanova(jaccard, sample_data(ps_rare)$ITW),
  
  Bray_Curtis_Heatwave = run_permanova(bray_curtis, sample_data(ps_rare)$heatwave),
  Jaccard_Heatwave = run_permanova(jaccard, sample_data(ps_rare)$heatwave),
  
  WUniFrac_Fraction = run_permanova(weighted_unifrac, sample_data(ps_rare)$fraction),
  UniFrac_Fraction = run_permanova(unweighted_unifrac, sample_data(ps_rare)$fraction),
  
  WUniFrac_Depth = run_permanova(weighted_unifrac, sample_data(ps_rare)$depth),
  UniFrac_Depth = run_permanova(unweighted_unifrac, sample_data(ps_rare)$depth),
  
  WUniFrac_ITW = run_permanova(weighted_unifrac, sample_data(ps_rare)$ITW),
  UniFrac_ITW = run_permanova(unweighted_unifrac, sample_data(ps_rare)$ITW),
  
  WUniFrac_Heatwave = run_permanova(weighted_unifrac, sample_data(ps_rare)$heatwave),
  UniFrac_Heatwave = run_permanova(unweighted_unifrac, sample_data(ps_rare)$heatwave)
)

# Print results
print(permanova_results)
extract_permanova_results <- function(permanova_list) {
  results <- data.frame(
    Metric = character(),
    Group = character(),
    Df = numeric(),
    SumOfSqs = numeric(),
    R2 = numeric(),
    F = numeric(),
    p_value = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE
  )
  
  significance_code <- function(p_value) {
    if (is.null(p_value) || is.na(p_value)) return("")
    if (p_value <= 0.001) return("***")
    if (p_value <= 0.01) return("**")
    if (p_value <= 0.05) return("*")
    if (p_value <= 0.1) return(".")
    return("")
  }
  
  for (metric in names(permanova_list)) {
    test_result <- permanova_list[[metric]]
    
    if (!is.null(test_result)) {
      # Extract group and method from the metric name
      group <- gsub("^(.*)_(.*)$", "\\2", metric)
      method <- gsub("^(.*)_(.*)$", "\\1", metric)
      
      # Extract values directly
      df <- test_result$Df[1]
      sum_of_sqs <- test_result$SumOfSqs[1]
      r2 <- test_result$R2[1]
      f_value <- test_result$F[1]
      p_value <- test_result$`Pr(>F)`[1]
      
      results <- rbind(
        results,
        data.frame(
          Metric = method,
          Group = group,
          Df = as.numeric(df),
          SumOfSqs = as.numeric(sum_of_sqs),
          R2 = as.numeric(r2),
          F = as.numeric(f_value),
          p_value = as.numeric(p_value),
          Significance = significance_code(p_value)
        )
      )
    } else {
      message(paste("No valid data for metric:", metric))
    }
  }
  
  return(results)
}
# Apply the function to permanova_results
table_data <- extract_permanova_results(permanova_results)

# Print the resulting table
print(table_data)

# Export the table to a CSV file for further formatting if needed
write.csv(table_data, "statistics/18S/18S_permanova_results_table.csv", row.names = FALSE)

# Define a function to run ANOSIM for a given distance matrix and grouping variable
run_anosim <- function(distance_matrix, grouping_variable) {
  # Create a data frame that contains the grouping variable and sample names
  sample_data_df <- data.frame(sample = rownames(sample_data(ps_rare)), 
                               group = grouping_variable)
  
  # Run ANOSIM
  anosim_result <- anosim(distance_matrix, sample_data_df$group)
  return(anosim_result)
}

# Run ANOSIM for Bray-Curtis and Jaccard distances
anosim_results <- list(
  Bray_Curtis_Fraction = run_anosim(bray_curtis, sample_data(ps_rare)$fraction),
  Jaccard_Fraction = run_anosim(jaccard, sample_data(ps_rare)$fraction),
  
  Bray_Curtis_Depth = run_anosim(bray_curtis, sample_data(ps_rare)$depth),
  Jaccard_Depth = run_anosim(jaccard, sample_data(ps_rare)$depth),
  
  Bray_Curtis_ITW = run_anosim(bray_curtis, sample_data(ps_rare)$ITW),
  Jaccard_ITW = run_anosim(jaccard, sample_data(ps_rare)$ITW),
  
  Bray_Curtis_Heatwave = run_anosim(bray_curtis, sample_data(ps_rare)$heatwave),
  Jaccard_Heatwave = run_anosim(jaccard, sample_data(ps_rare)$heatwave),
  
  WUniFrac_Fraction = run_anosim(weighted_unifrac, sample_data(ps_rare)$fraction),
  UniFrac_Fraction = run_anosim(unweighted_unifrac, sample_data(ps_rare)$fraction),
  
  WUniFrac_Depth = run_anosim(weighted_unifrac, sample_data(ps_rare)$depth),
  UniFrac_Depth = run_anosim(unweighted_unifrac, sample_data(ps_rare)$depth),
  
  WUniFrac_ITW = run_anosim(weighted_unifrac, sample_data(ps_rare)$ITW),
  UniFrac_ITW = run_anosim(unweighted_unifrac, sample_data(ps_rare)$ITW),
  
  WUniFrac_Heatwave = run_anosim(weighted_unifrac, sample_data(ps_rare)$heatwave),
  UniFrac_Heatwave = run_anosim(unweighted_unifrac, sample_data(ps_rare)$heatwave)
)

# Print results
print(anosim_results)
extract_anosim_results <- function(anosim_list) {
  results <- data.frame(
    Metric = character(),
    Group = character(),
    R_statistic = numeric(),
    p_value = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE
  )
  
  significance_code <- function(p_value) {
    if (is.null(p_value) || is.na(p_value)) return("")
    if (p_value <= 0.001) return("***")
    if (p_value <= 0.01) return("**")
    if (p_value <= 0.05) return("*")
    if (p_value <= 0.1) return(".")
    return("")
  }
  
  for (metric in names(anosim_list)) {
    test_result <- anosim_list[[metric]]
    
    if (!is.null(test_result)) {
      # Extract group and method from the metric name
      group <- gsub("^(.*)_(.*)$", "\\2", metric)
      method <- gsub("^(.*)_(.*)$", "\\1", metric)
      
      # Extract values from the ANOSIM result
      r_statistic <- test_result$statistic
      p_value <- test_result$signif
      
      results <- rbind(
        results,
        data.frame(
          Metric = method,
          Group = group,
          R_statistic = as.numeric(r_statistic),
          p_value = as.numeric(p_value),
          Significance = significance_code(p_value)
        )
      )
    } else {
      message(paste("No valid data for metric:", metric))
    }
  }
  
  return(results)
}
anosim_table <- extract_anosim_results(anosim_results)
print(anosim_table)
# Export the table to a CSV file for further formatting if needed
write.csv(anosim_table, "statistics/18S/18S_anosim_results_table.csv", row.names = FALSE)

# Define a common color palette and other aesthetic elements
common_color <- scale_color_manual(values = c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F"))
common_shape <- scale_shape_manual(values = c(16, 1), name = "Depth")

# Step 1: Compute nMDS for Bray-Curtis
set.seed(123)  # For reproducibility
nmds_bray <- metaMDS(bray_curtis, k = 2)

# Step 2: Prepare Data for Bray-Curtis Plot
nmds_bray_df <- as.data.frame(nmds_bray$points)
nmds_bray_df$fraction <- metadata$fraction
nmds_bray_df$depth <- metadata$depth

# Step 3: Create the nMDS Plot for Bray-Curtis
p5<- ggplot(nmds_bray_df, aes(x = MDS1, y = MDS2, color = fraction, shape = depth)) +
  geom_point(size = 3) +
  common_color +  # Add custom color scale for groups
  common_shape +  # Add custom shape scale for depth
  theme_q2r() +
  labs(title = "nMDS Plot - Bray-Curtis", x = "nMDS 1", y = "nMDS 2") +
  theme(legend.position = "right")
p5
# Step 4: Compute nMDS for Jaccard
set.seed(123)  # For reproducibility
nmds_jaccard <- metaMDS(jaccard, k = 2)

# Step 5: Prepare Data for Jaccard Plot
nmds_jaccard_df <- as.data.frame(nmds_jaccard$points)
nmds_jaccard_df$fraction <- metadata$fraction
nmds_jaccard_df$depth <- metadata$depth

# Step 6: Create the nMDS Plot for Jaccard
p6<- ggplot(nmds_jaccard_df, aes(x = MDS1, y = MDS2, color = fraction, shape = depth)) +
  geom_point(size = 3) +
  common_color +  # Add custom color scale for groups
  common_shape +  # Add custom shape scale for depth
  theme_q2r() +
  labs(title = "nMDS Plot - Jaccard", x = "nMDS 1", y = "nMDS 2") +
  theme(legend.position = "right")
p6
# Step 7: Compute nMDS for WUniFrac
set.seed(123)  # For reproducibility
nmds_wunifrac <- metaMDS(weighted_unifrac, k = 2)

# Step 8: Prepare Data for WUniFrac Plot
nmds_wunifrac_df <- as.data.frame(nmds_wunifrac$points)
nmds_wunifrac_df$fraction <- metadata$fraction
nmds_wunifrac_df$depth <- metadata$depth

# Step 9: Create the nMDS Plot for WUniFrac
p7<- ggplot(nmds_wunifrac_df, aes(x = MDS1, y = MDS2, color = fraction, shape = depth)) +
  geom_point(size = 3) +
  common_color +  # Add custom color scale for groups
  common_shape +  # Add custom shape scale for depth
  theme_q2r() +
  labs(title = "nMDS Plot - WUniFrac", x = "nMDS 1", y = "nMDS 2") +
  theme(legend.position = "right")
p7
# Step 10: Compute nMDS for UniFrac
set.seed(123)  # For reproducibility
nmds_unifrac <- metaMDS(unweighted_unifrac, k = 2)

# Step 11: Prepare Data for UniFrac Plot
nmds_unifrac_df <- as.data.frame(nmds_unifrac$points)
nmds_unifrac_df$fraction <- metadata$fraction
nmds_unifrac_df$depth <- metadata$depth

# Step 12: Create the nMDS Plot for UniFrac
p8<- ggplot(nmds_unifrac_df, aes(x = MDS1, y = MDS2, color = fraction, shape = depth)) +
  geom_point(size = 3) +
  common_color +  # Add custom color scale for groups
  common_shape +  # Add custom shape scale for depth
  theme_q2r() +
  labs(title = "nMDS Plot - UniFrac", x = "nMDS 1", y = "nMDS 2") +
  theme(legend.position = "right")
p8

# Extract individual legends
library(cowplot)
legend_color <- get_legend(
  p5 + theme(legend.position = "right")
)

legend_shape_size <- get_legend(
  p5 + 
    guides(
      size = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    theme(legend.position = "right")
)


# Combine the plots
combined_plots <- plot_grid(
  p5, p6, p7, p8,
  labels = c("A", "B", "C", "D"), 
  ncol = 2, align = "v"
)
# Combine the legends into a single legend
combined_legends <- plot_grid(
  legend_color, legend_shape_size, 
  ncol = 1, rel_heights = c(1, 1)
)
# Final combination of plots and legends
final_plot <- plot_grid(
  combined_plots, combined_legends, 
  ncol = 2, rel_widths = c(1, 0.1)
)
final_plot
# Save the final plot
ggsave("plots/18S/Combined_nMDS_Plots_Bray_Jaccard_WUnifrac_Unifrac.pdf", final_plot, height = 7, width = 8, device = "pdf")

# Beta diversity variation with fraction, ITW, depth
# Perform beta diversity variation analysis
tidy_taxonomy <- function(tax_table) {
  # Assuming you just want to clean the tax_table or convert it to a tibble
  return(as_tibble(tax_table))
}

dataset <- phyloseq2meco(ps_rare)
dataset$cal_betadiv(unifrac = FALSE)

# Calculate group distances and significance
t1 <- trans_beta$new(dataset = dataset, group = "fraction", measure = "bray")
t1$cal_group_distance(within_group = TRUE)
t1$cal_group_distance_diff(method = "anova") # Significance with ANOVA

# Prepare plot data
g1 <- t1$plot_group_distance(boxplot_add = "mean", color_values = "#DDAA33", xtext_keep = TRUE)
g1$data$fraction <- factor(g1$data$fraction, levels = c("MICRO", "NANO", "PICO"))

# Customize plot with custom colors and improved visuals
custom_colors <- c("MICRO" = "#FB9A99", "NANO" = "#A6CEE3", "PICO" = "#FDBF6F")
g1 <- g1 +
  geom_boxplot(aes(fill = fraction)) +  # Map fill to 'fraction' for custom colors
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(limits = c("MICRO", "NANO", "PICO")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  coord_flip() +
  ggtitle("18S 5947 ASV Bray-Curtis Dissimilarity by Fraction")
print(g1)
# Save the plot
ggsave("plots/18S/18S_Bray_box-whiskers_5947ASV_fraction.pdf", plot = g1, height = 4, width = 6, device = "pdf")

# Calculate Bray-Curtis distances for 'depth'
t1_depth <- trans_beta$new(dataset = dataset, group = "depth", measure = "bray")
t1_depth$cal_group_distance(within_group = TRUE)
t1_depth$cal_group_distance_diff(method = "anova")

g1_depth <- t1_depth$plot_group_distance(boxplot_add = "mean", color_values = "#DDAA33", xtext_keep = TRUE)
g1_depth$data$depth <- factor(g1_depth$data$depth, levels = c("SURFACE", "DCM"))

# Customize plot with custom colors for depth
custom_colors_depth <- c("SURFACE" = "red3", "DCM" = "blue3")

# Plot for Depth
g_depth <- g1_depth +
  geom_boxplot(aes(fill = depth)) +  # Map fill to 'depth' for custom colors
  scale_fill_manual(values = custom_colors_depth) +
  scale_x_discrete(limits = c("SURFACE", "DCM")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  coord_flip() +
  ggtitle("18S 5947 ASV Bray-Curtis Dissimilarity by Depth")
print(g_depth)
# Save the plot
ggsave("plots/18S/18S_Bray_box-whiskers_5947ASV_depth.pdf", plot = g_depth, height = 4, width = 6, device = "pdf")

#Calculate Bray-Curtis variation according to ITW
t1 <- trans_beta$new(dataset = dataset, group = "ITW", measure = "bray")
t1$cal_group_distance(within_group = T)
t1$cal_group_distance_diff(method = "anova") # Significance with ANOVA
g1 <- t1$plot_group_distance(boxplot_add = "mean",color_values = "#DDAA33",xtext_keep = TRUE)
g1$data$ITW <- factor(g1$data$ITW, levels = c("yes", "no"))

#Plot for ITW
g1 + 
  geom_boxplot(fill="#DDAA33") +
  scale_x_discrete(limits = c("yes", "no")) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12)) +
  coord_flip() +
  ggtitle("18S 5947ASV Bray-Curtis dissimilarity_ITW")
ggsave("plots/18S/Bray_box-whiskers_5947features_ITW.pdf", height=2, width=6, device="pdf")

####Community composition plots
####Preparing Class level for stacked barplots over depth & fraction
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
ps_rare <- readRDS("data/phyloseq_18S_rare.rds")
# 1. Transform OTU counts to relative abundances
ps_rel <- transform_sample_counts(ps_rare, function(x) x / sum(x))

# 2. Melt the phyloseq object into a long format data frame
ps_melted <- psmelt(ps_rel)

ps_rel_ITW <- ps_melted %>% filter(ITW == "yes")
ps_rel_noITW <- ps_melted %>% filter(ITW == "no")

# 3. Calculate top 15 Class for ITW == "yes" subset
top_Class_ITW <- ps_rel_ITW %>%
  group_by(Class) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  top_n(15, total_abundance) %>%
  pull(Class)

# Filter to top 15 Classes and group the rest as "Other"
ps_rel_ITW_top15 <- ps_rel_ITW %>%
  mutate(Class = ifelse(Class %in% top_Class_ITW, Class, "Other"))

# 4. Calculate total abundance for each fraction and depth combination
total_abundance_ITW <- ps_rel_ITW_top15 %>%
  group_by(fraction, depth) %>%
  summarise(Total = sum(Abundance), .groups = 'drop')

# 5. Group and summarize abundances for ITW == "yes"
ps_rel_ITW_top15_summary <- ps_rel_ITW_top15 %>%
  group_by(fraction, depth, Class) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  left_join(total_abundance_ITW, by = c("fraction", "depth")) %>%
  mutate(Relative_Abundance = (Abundance / Total) * 100)

# Remove rows where Relative_Abundance is NA or 0
ps_rel_ITW_top15_summary <- ps_rel_ITW_top15_summary %>%
  filter(!is.na(Relative_Abundance) & Relative_Abundance > 0)

# 6. Calculate top 15 Class for ITW == "no" subset
top_Class_noITW <- ps_rel_noITW %>%
  group_by(Class) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  top_n(15, total_abundance) %>%
  pull(Class)

# Filter to top 15 Classes and group the rest as "Other"
ps_rel_noITW_top15 <- ps_rel_noITW %>%
  mutate(Class = ifelse(Class %in% top_Class_noITW, Class, "Other"))

# 7. Calculate total abundance for each fraction and depth combination
total_abundance_noITW <- ps_rel_noITW_top15 %>%
  group_by(fraction, depth) %>%
  summarise(Total = sum(Abundance), .groups = 'drop')

# 8. Group and summarize abundances for ITW == "no"
ps_rel_noITW_top15_summary <- ps_rel_noITW_top15 %>%
  group_by(fraction, depth, Class) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  left_join(total_abundance_noITW, by = c("fraction", "depth")) %>%
  mutate(Relative_Abundance = (Abundance / Total) * 100)

# Remove rows where Relative_Abundance is NA or 0
ps_rel_noITW_top15_summary <- ps_rel_noITW_top15_summary %>%
  filter(!is.na(Relative_Abundance) & Relative_Abundance > 0)

# 9. Define a color palette for top 15 Classes + "Other"
class_colors_top15 <- c(
  "#A6CEE3", "red4", "#CAB2D6", "#6A3D9A" ,"mediumseagreen", 
  "#1F78B4", "#FDBF6F", "#FF7F00", "#FFFF99", "gray50", 
  "orchid", "khaki4", "peachpuff","#FB9A99" , "yellow", 
  "lawngreen"  # "Other" will be in gray50
)

# 10. Plot for ITW == "yes"
g_ITW <- ggplot(ps_rel_ITW_top15_summary, 
                aes(x = fraction, y = Relative_Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, scales = "free_y") +
  coord_flip() +  # Horizontal bars
  labs(title = "Relative Abundance of Top Class for ITW = 'Yes'",
       x = "Fraction",
       y = "Relative Abundance (%)",
       fill = "Class") +
  scale_fill_manual(values = class_colors_top15) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))

# 11. Plot for ITW == "no"
g_noITW <- ggplot(ps_rel_noITW_top15_summary, 
                  aes(x = fraction, y = Relative_Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, scales = "free_y") +
  coord_flip() +  # Horizontal bars
  labs(title = "Relative Abundance of Top Class for ITW = 'No'",
       x = "Fraction",
       y = "Relative Abundance (%)",
       fill = "Class") +
  scale_fill_manual(values = class_colors_top15) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))

# 12. Print the plots
print(g_ITW)
print(g_noITW)

# Load the patchwork package for combining plots
library(patchwork)

# Adjust the bar width, theme, and font sizes
g_ITW <- g_ITW +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +  # Adjust bar width
  labs(tag = "A") +  # Add panel label
  theme_q2r() +  # Apply theme_q2r
  theme(
    legend.key.size = unit(0.3, "cm"),         # Reduce legend size
    legend.text = element_text(size = 10),    # Adjust legend text size
    legend.title = element_text(size = 12),   # Adjust legend title size
    strip.text = element_text(size = 14),     # Increase font size of facet wraps
    axis.text = element_text(size = 12),      # Increase font size of axes
    axis.title = element_text(size = 14)      # Increase font size of axis titles
  )

g_noITW <- g_noITW +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +  # Adjust bar width
  labs(tag = "B") +  # Add panel label
  theme_q2r() +  # Apply theme_q2r
  theme(
    legend.key.size = unit(0.3, "cm"),         # Reduce legend size
    legend.text = element_text(size = 10),    # Adjust legend text size
    legend.title = element_text(size = 12),   # Adjust legend title size
    strip.text = element_text(size = 14),     # Increase font size of facet wraps
    axis.text = element_text(size = 12),      # Increase font size of axes
    axis.title = element_text(size = 14)      # Increase font size of axis titles
  )

# Combine the two plots vertically
combined_plot <- g_ITW - g_noITW

# Display the combined plot
print(combined_plot)

# Step 1: List unique Classes for each condition
classes_ITW <- unique(ps_rel_ITW_top15_summary$Class)
classes_noITW <- unique(ps_rel_noITW_top15_summary$Class)

# Step 2: Identify matching and mismatched Classes
matching_classes <- intersect(classes_ITW, classes_noITW)
mismatched_ITW <- setdiff(classes_ITW, matching_classes)
mismatched_noITW <- setdiff(classes_noITW, matching_classes)

# Step 3: Assign colors using `class_colors_top15`
# Matching Classes get the first colors in order
matching_colors <- setNames(class_colors_top15[1:length(matching_classes)], matching_classes)

# Remaining colors are used for mismatched Classes
# Recycle colors if mismatched Classes exceed the remaining palette
remaining_colors <- class_colors_top15[(length(matching_classes) + 1):length(class_colors_top15)]
if (length(remaining_colors) < length(mismatched_ITW) + length(mismatched_noITW)) {
  remaining_colors <- rep(remaining_colors, length.out = length(mismatched_ITW) + length(mismatched_noITW))
}

# Split remaining colors for mismatched Classes in ITW and noITW
mismatched_colors_ITW <- setNames(remaining_colors[1:length(mismatched_ITW)], mismatched_ITW)
mismatched_colors_noITW <- setNames(remaining_colors[(length(mismatched_ITW) + 1):(length(mismatched_ITW) + length(mismatched_noITW))], mismatched_noITW)

# Combine all colors into palettes for each condition
palette_ITW <- c(matching_colors, mismatched_colors_ITW)
palette_noITW <- c(matching_colors, mismatched_colors_noITW)

# Step 4: Plot for ITW == "yes" using its palette
g_ITW <- ggplot(ps_rel_ITW_top15_summary, 
                aes(x = fraction, y = Relative_Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, scales = "free_y") +
  coord_flip() +  # Horizontal bars
  labs(title = "Relative Abundance of Top Class for ITW = 'Yes'",
       x = "Fraction",
       y = "Relative Abundance (%)",
       fill = "Class") +
  scale_fill_manual(values = palette_ITW) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))

# Step 5: Plot for ITW == "no" using its palette
g_noITW <- ggplot(ps_rel_noITW_top15_summary, 
                  aes(x = fraction, y = Relative_Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, scales = "free_y") +
  coord_flip() +  # Horizontal bars
  labs(title = "Relative Abundance of Top Class for ITW = 'No'",
       x = "Fraction",
       y = "Relative Abundance (%)",
       fill = "Class") +
  scale_fill_manual(values = palette_noITW) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))

# Step 6: Combine the two plots vertically
combined_plot <- g_ITW / g_noITW

# Display the combined plot
print(combined_plot)
