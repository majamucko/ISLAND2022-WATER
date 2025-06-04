###Genus eukaryotic networks in SURFACE and DCM layers

###Data for these networks were derived from rarefied phyloseq object: 
#phyloseq_18S <- readRDS("data/phyloseq_18S_rare.rds")

###Phyloseq object was agglomerated to genus level with (agglomerateByRank) function
###from tse package
###Table was extracted from R, manually checked for terrigenus and benthic taxa and annotated with trophic_mode, body size and size class according to literature
###Final dataset for SURFACE encompasses 204 genus, for DCM 225 genus divided to: 
###"autotrophic protist", "heterotrophic protist", "mixotrophic protist", "parasitic protist" and "zooplankton"

###Network calculation
###SURFACE
#Load necessery packages
library(SpiecEasi)
library(mia)
library(TreeSummarizedExperiment)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)

#functions
# Arguments:
# - assoMat: association matrix
# - threshold: associations below the threshold are set to zero
# - dissTrans: dissimilarity transformation ("signed" or "unsigned")

transform_asso <- function(assoMat, thresh = NULL, dissTrans = "signed") {
  # Sparsification
  if (!is.null(thresh)) {
    assoMat[abs(assoMat) < thresh] <- 0
  }
  
  # Compute dissimilarity matrix
  if (dissTrans == "signed") {
    dissMat <- sqrt(0.5 * (1 - assoMat))
  } else {
    dissMat <- sqrt(1 - assoMat^2)
  }
  
  # Dissimilarity between nodes with zero correlation is set to 1
  # (these nodes are unconnected and thus should have maximum dissimilarity)
  dissMat[assoMat == 0] <- 1
  
  # Compute similarity matrix
  simMat <- 1 - dissMat
  
  # Turn into igraph object
  graphObj <- SpiecEasi::adj2igraph(simMat)
  
  return(list(graph = graphObj, adja = simMat, asso = assoMat, diss = dissMat))
}

#Load data
SURFACE_net <- read.table("data/SURFACE_18S_204genus_matrix_annotated.tsv", header = TRUE, sep = "\t")
print(SURFACE_net)
head(SURFACE_net)
# Set ASV IDs as row names and remove the ASV_ID column
rownames(SURFACE_net) <- SURFACE_net$Taxa
SURFACE_net$Taxa <- NULL
print(SURFACE_net)

# Extract OTU count columns (numeric part only)
otu_counts_SURFACE <- SURFACE_net[, 1:12]

# Convert OTU counts to numeric matrix
otu_counts_SURFACE <- as.matrix(otu_counts_SURFACE)

# Normalize the OTU data by converting counts to relative abundances
otu_counts_SURFACE_norm <- t(otu_counts_SURFACE)  # Transpose for samples as rows
otu_counts_SURFACE_norm <- sweep(otu_counts_SURFACE_norm, 1, rowSums(otu_counts_SURFACE_norm), "/")  # Normalize by row sums

# Ensure the normalized data is numeric
otu_counts_SURFACE_norm <- as.matrix(otu_counts_SURFACE_norm)
head(otu_counts_SURFACE_norm)
# Check for NAs or Infs
if (any(is.na(otu_counts_SURFACE_norm)) || any(is.infinite(otu_counts_SURFACE_norm))) {
  stop("Data contains NA or Inf values, please clean data.")
}

# Log transformation (optional, prevents log(0) errors)
otu_counts_SURFACE_log <- log1p(otu_counts_SURFACE_norm)  # log1p handles log(1 + x), safe for 0 values

# Ensure data is numeric matrix
otu_counts_SURFACE_log <- as.matrix(otu_counts_SURFACE_log)

# Build the network using SPIEC-EASI with pulsar parameters
set.seed(13075)  # For reproducibility
se_mb_SURFACE <- spiec.easi(
  otu_counts_SURFACE_log, 
  method = 'mb', 
  nlambda = 50, 
  lambda.min.ratio = 1e-3, 
  pulsar.params = list(rep.num = 50)
)

# Inspect the result
print(se_mb_SURFACE)# Inspect the output
str(se_mb_SURFACE)

#Network analysis
# Get optimal matrix with partial correlations
se_mb_SURFACE_cor <- as.matrix(getOptBeta(se_mb_SURFACE))
se_mb_SURFACE_cor <- as.matrix(symBeta(se_mb_SURFACE_cor))
rownames(se_mb_SURFACE_cor) <- colnames(se_mb_SURFACE_cor) <- rownames(otu_counts_SURFACE)
diag(se_mb_SURFACE_cor) <- 1

# Extract edges from the correlation matrix
edges_SURFACE <- as.data.frame(as.table(se_mb_SURFACE_cor))

# Remove self-loops (correlations of taxa with themselves)
edges_SURFACE <- edges_SURFACE[edges_SURFACE$Var1 != edges_SURFACE$Var2, ]

# Rename columns for clarity
colnames(edges_SURFACE) <- c("Source", "Target", "Weight")

# Separate positive and negative links
positive_links <- subset(edges_SURFACE, Weight > 0)
negative_links <- subset(edges_SURFACE, Weight < 0)

# Save results to files
write.table(positive_links, "SPIEC_EASI_positive_links_SURFACE_Genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(negative_links, "SPIEC_EASI_negative_links_SURFACE_Genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Print summaries
cat("Number of positive links:", nrow(positive_links), "\n")
cat("Number of negative links:", nrow(negative_links), "\n")

################################################################################################
################################################################################################

### DCM Network Calculation
# Load data
DCM_net <- read.table("data/DCM_18S_225genus_matrix_annotated.tsv", header = TRUE, sep = "\t")
print(DCM_net)
head(DCM_net)

# Set ASV IDs as row names and remove the ASV_ID column
rownames(DCM_net) <- DCM_net$Taxa
DCM_net$Taxa <- NULL
print(DCM_net)

# Extract OTU count columns (numeric part only)
otu_counts_DCM <- DCM_net[, 1:12]

# Convert OTU counts to numeric matrix
otu_counts_DCM <- as.matrix(otu_counts_DCM)

# Normalize the OTU data by converting counts to relative abundances
otu_counts_DCM_norm <- t(otu_counts_DCM)  # Transpose for samples as rows
otu_counts_DCM_norm <- sweep(otu_counts_DCM_norm, 1, rowSums(otu_counts_DCM_norm), "/")  # Normalize by row sums

# Ensure the normalized data is numeric
otu_counts_DCM_norm <- as.matrix(otu_counts_DCM_norm)
head(otu_counts_DCM_norm)

# Check for NAs or Infs
if (any(is.na(otu_counts_DCM_norm)) || any(is.infinite(otu_counts_DCM_norm))) {
  stop("Data contains NA or Inf values, please clean data.")
}

# Log transformation (optional, prevents log(0) errors)
otu_counts_DCM_log <- log1p(otu_counts_DCM_norm)  # log1p handles log(1 + x), safe for 0 values

# Ensure data is numeric matrix
otu_counts_DCM_log <- as.matrix(otu_counts_DCM_log)

# Build the network using SPIEC-EASI with pulsar parameters
set.seed(13075)  # For reproducibility
se_mb_DCM <- spiec.easi(
  otu_counts_DCM_log, 
  method = 'mb', 
  nlambda = 50, 
  lambda.min.ratio = 1e-3, 
  pulsar.params = list(rep.num = 50)
)

# Inspect the result
print(se_mb_DCM)  # Inspect the output
str(se_mb_DCM)

# Network analysis
# Get optimal matrix with partial correlations
se_mb_DCM_cor <- as.matrix(getOptBeta(se_mb_DCM))
se_mb_DCM_cor <- as.matrix(symBeta(se_mb_DCM_cor))
rownames(se_mb_DCM_cor) <- colnames(se_mb_DCM_cor) <- rownames(otu_counts_DCM)
diag(se_mb_DCM_cor) <- 1

# Extract edges from the correlation matrix
edges_DCM <- as.data.frame(as.table(se_mb_DCM_cor))

# Remove self-loops (correlations of taxa with themselves)
edges_DCM <- edges_DCM[edges_DCM$Var1 != edges_DCM$Var2, ]

# Rename columns for clarity
colnames(edges_DCM) <- c("Source", "Target", "Weight")

# Separate positive and negative links
positive_links_DCM <- subset(edges_DCM, Weight > 0)
negative_links_DCM <- subset(edges_DCM, Weight < 0)

# Save results to files
write.table(positive_links_DCM, "SPIEC_EASI_positive_links_DCM_Genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(negative_links_DCM, "SPIEC_EASI_negative_links_DCM_Genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Print summaries
cat("Number of positive links:", nrow(positive_links_DCM), "\n")
cat("Number of negative links:", nrow(negative_links_DCM), "\n")

library(igraph)
library(ggraph)
library(scales)

# Create an igraph object from the corr matrix
# Get the signed correlation matrix
se_mb_SURFACE_cor <- as.matrix(getOptBeta(se_mb_SURFACE))

# Convert to a symmetric format (ensuring all edges are bidirectional)
se_mb_SURFACE_cor <- as.matrix(symBeta(se_mb_SURFACE_cor))

# Set diagonal values to 1 (optional, depends on visualization needs)
diag(se_mb_SURFACE_cor) <- 1

# Check range of values
range(se_mb_SURFACE_cor)  # Now you should see negative and positive values!
graph_SURFACE <- graph_from_adjacency_matrix(se_mb_SURFACE_cor, mode = "directed", weighted = TRUE)
V(graph_SURFACE)$name <- SURFACE_net$Taxonomy
V(graph_SURFACE)$Trophy <- SURFACE_net$Trophy
summary(graph_SURFACE)

# Assign edge colors based on weight sign
E(graph_SURFACE)$color <- ifelse(E(graph_SURFACE)$weight > 0, "blue", "red")

# Set a weight threshold (e.g., 0.5)
#weight_threshold <- 0.4

# Filter edges based on the threshold
#graph_SURFACE <- delete_edges(graph_SURFACE, E(graph_SURFACE)[abs(weight) < weight_threshold])

# Calculate node degree
V(graph_SURFACE)$degree <- degree(graph_SURFACE)

# Define a custom color map for Trophy levels
color_map <- c(
  "zooplankton" = "red3",
  "heterotrophic protist" = "purple",
  "autotrophic protist" = "green2",
  "mixotrophic protist" = "blue2",
  "parasitic protist" = "black"
)

# Scale node size by degree
p<- ggraph(graph_SURFACE, layout = "graphopt") +  
  geom_edge_link(aes(edge_width = abs(weight), color = factor(sign(weight))), alpha = 0.4) +
  scale_edge_color_manual(values = c("-1" = "red", "1" = "blue"), name = "Link Sign") +  
  scale_edge_width_continuous(range = c(0.2, 2)) +  
  geom_node_point(aes(size = degree, color = Trophy), alpha = 0.4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, force = 10) +
  scale_size(range = c(1, 5), name = "Node Degree") +  # Explicit legend title
  scale_color_manual(values = color_map, name = "Trophic Mode") +
  theme_graph() +
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.margin = margin(0.2, 0.2, 0.2, 0.2),
    legend.box.spacing = unit(0.2, "lines")
  ) +
  ggtitle("SPIEC-EASI SURFACE Genus Network")
# Save the plot as an SVG file
install.packages("svglite")
library(svglite)
ggsave("SPIEC-EASI SURFACE Genus_Network.svg", plot = p, width = 12, height = 8)

library(igraph)
library(ggraph)
library(scales)

# Create an igraph object from the corr matrix
# Get the signed correlation matrix
se_mb_DCM_cor <- as.matrix(getOptBeta(se_mb_DCM))

# Convert to a symmetric format (ensuring all edges are bidirectional)
se_mb_DCM_cor <- as.matrix(symBeta(se_mb_DCM_cor))

# Set diagonal values to 1 (optional, depends on visualization needs)
diag(se_mb_DCM_cor) <- 1

# Check range of values
range(se_mb_DCM_cor)  # Now you should see negative and positive values!
graph_DCM <- graph_from_adjacency_matrix(se_mb_DCM_cor, mode = "directed", weighted = TRUE)
V(graph_DCM)$name <- DCM_net$Taxonomy
V(graph_DCM)$Trophy <- DCM_net$Trophy
summary(graph_DCM)

# Assign edge colors based on weight sign
E(graph_DCM)$color <- ifelse(E(graph_DCM)$weight > 0, "blue", "red")

# Set a weight threshold (e.g., 0.5)
#weight_threshold <- 0.4

# Filter edges based on the threshold
#graph_SURFACE <- delete_edges(graph_SURFACE, E(graph_SURFACE)[abs(weight) < weight_threshold])

# Calculate node degree
V(graph_DCM)$degree <- degree(graph_DCM)

# Define a custom color map for Trophy levels
color_map <- c(
  "zooplankton" = "red3",
  "heterotrophic protist" = "purple",
  "autotrophic protist" = "green2",
  "mixotrophic protist" = "blue2",
  "parasitic protist" = "black"
)

# Scale node size by degree
p1<- ggraph(graph_DCM, layout = "graphopt") +  
  geom_edge_link(aes(edge_width = abs(weight), color = factor(sign(weight))), alpha = 0.4) +
  scale_edge_color_manual(values = c("-1" = "red", "1" = "blue"), name = "Link Sign") +  
  scale_edge_width_continuous(range = c(0.2, 2)) +  
  geom_node_point(aes(size = degree, color = Trophy), alpha = 0.4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, force = 10) +
  scale_size(range = c(1, 5), name = "Node Degree") +  # Explicit legend title
  scale_color_manual(values = color_map, name = "Trophic Mode") +
  theme_graph() +
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.margin = margin(0.2, 0.2, 0.2, 0.2),
    legend.box.spacing = unit(0.2, "lines")
  ) +
  ggtitle("SPIEC-EASI DCM Genus Network")
# Save the plot as an SVG file
install.packages("svglite")
library(svglite)
ggsave("SPIEC-EASI DCM Genus_Network.svg", plot = p1, width = 12, height = 8)