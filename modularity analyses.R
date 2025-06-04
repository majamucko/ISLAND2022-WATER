library(igraph)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

# DCM GENUS
dcm_genus_edges <- read.csv2("SPIEC_EASI_links_DCM_genus.csv")
dcm_genus_nodes <- read.csv2("SPIEC_EASI_nodes_DCM_genus.csv")

network <- graph_from_data_frame(dcm_genus_edges, dcm_genus_nodes, directed = F) 

set.seed(111)
# Calculation of modularity and modularity roles
degree=degree(network)
Nnodes=length(degree)
Z=degree
Z[]=0
P=Z
Wtc=cluster_louvain(network, resolution = 0.5)
Membership=membership(Wtc)
Seq=seq(1:Nnodes)
for(i in 1:Nnodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(network,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(network,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/degree[i])^2}
  P[i]=1-SUMP
}


dcm_attribute_node <- dcm_genus_nodes %>%
  mutate(modularity = Membership,
         Pi = P,
         Zi = Z,
         role = case_when(
           Zi > 2.5 & Pi <= 0.62 ~ "module hubs",
           Zi <= 2.5 & Pi > 0.62 ~ "connectors",
           Zi > 2.5 & Pi > 0.62 ~ "network hubs",
           Zi <= 2.5 & Pi <= 0.62 ~ "peripheral"
         ))

# Zi-Pi plot
ggplot() + 
  geom_point(mapping = aes(x = as.numeric(dcm_attribute_node$Pi), y = as.numeric(dcm_attribute_node$Zi), color = as.factor(dcm_attribute_node$Trophy)), size = 3.5)  + 
  theme_classic()  + 
  geom_hline(yintercept=2.5, linetype = "dashed") + 
  geom_vline(xintercept=0.62, linetype = "dashed") +  
  annotate("text", x = 0.1, y = 3, label = "Module hubs", size = 5) +
  annotate("text", x = 0.7, y = 3, label = "Network hubs", size = 5) +
  annotate("text", x = 0.1, y = -1.5, label = "Peripherals", size = 5) +
  annotate("text", x = 0.7, y = -1.5, label = "Connectors", size = 5) +
  theme_minimal() +  
  labs(x = "Among-module connectivity (Pi)",
       y = "Within-module connectivity (Zi)",
       color = "Functional groups") #+ # To visualize node labels
  geom_text(mapping = aes(x = as.numeric(dcm_attribute_node$Pi), y = as.numeric(dcm_attribute_node$Zi), label = dcm_attribute_node$Label)) 

# Calculate percentage of size classes and trophic modes within each module
size_percentages_dcm <- dcm_attribute_node %>%
    group_by(modularity, size_class) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(modularity) %>%
    mutate(percentage = round(100 * count / sum(count), 2)) %>%
    select(-count) %>%
    pivot_wider(names_from = size_class, values_from = percentage, values_fill = 0)
  
trophism_percentages_dcm <- dcm_attribute_node %>%
    group_by(modularity, Trophy) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(modularity) %>%
    mutate(percentage = round(100 * count / sum(count), 2)) %>%
    select(-count) %>%
    pivot_wider(names_from = Trophy, values_from = percentage, values_fill = 0)
  
aggregated_dcm <- size_percentages_dcm %>%
    left_join(trophism_percentages_dcm, by = "modularity")

# Convert size data to long format
size_long_dcm <- aggregated_dcm  %>%
    pivot_longer(cols = meso:mega, names_to = "Size_Category", values_to = "Percentage")
  
# Convert trophic mode data to long format
trophic_long_dcm  <- aggregated_dcm  %>%
    pivot_longer(cols = "autotrophic protist":zooplankton, names_to = "Trophic_mode", values_to = "Percentage")
  
# Size category plot
plot_size_dcm <- ggplot(size_long_dcm , aes(x = factor(modularity), y = Percentage, fill = Size_Category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Modularity", y = "Percentage", fill = "Size Class")  +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  
# Trophic mode plot
plot_trophic_dcm  <- ggplot(trophic_long_dcm , aes(x = factor(modularity), y = Percentage, fill = Trophic_mode)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Modularity", y = "Percentage", fill = "Functional Group")  +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")
  
# Combine plots with cowplot
plot_grid(plot_size_dcm , plot_trophic_dcm , ncol = 1, align = "v")



# SURFACE GENUS
surface_genus_edges <- read.csv2("SPIEC_EASI_links_SURFACE_genus.csv")
surface_genus_nodes <- read.csv2("SPIEC_EASI_nodes_SURFACE_genus.csv")

network <- graph_from_data_frame(surface_genus_edges, surface_genus_nodes, directed = F) 

set.seed(111)
# Calculation of modularity and modularity roles
degree=degree(network)
Nnodes=length(degree)
Z=degree
Z[]=0
P=Z
Wtc=cluster_louvain(network, resolution = 0.5)
Membership=membership(Wtc)
Seq=seq(1:Nnodes)
for(i in 1:Nnodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(network,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(network,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/degree[i])^2}
  P[i]=1-SUMP
}


surface_attribute_node <- surface_genus_nodes %>%
  mutate(modularity = Membership,
         Pi = P,
         Zi = Z,
         role = case_when(
           Zi > 2.5 & Pi <= 0.62 ~ "module hubs",
           Zi <= 2.5 & Pi > 0.62 ~ "connectors",
           Zi > 2.5 & Pi > 0.62 ~ "network hubs",
           Zi <= 2.5 & Pi <= 0.62 ~ "peripheral"
         ))

# Zi-Pi plot
ggplot() + 
  geom_point(mapping = aes(x = as.numeric(surface_attribute_node$Pi), y = as.numeric(surface_attribute_node$Zi), color = surface_attribute_node$Trophy), size = 3.5)  + 
  theme_classic()  + 
  geom_hline(yintercept=2.5, linetype = "dashed") + 
  geom_vline(xintercept=0.62, linetype = "dashed") +  
  annotate("text", x = 0.1, y = 3, label = "Module hubs", size = 5) +
  annotate("text", x = 0.7, y = 3, label = "Network hubs", size = 5) +
  annotate("text", x = 0.1, y = -1.5, label = "Peripherals", size = 5) +
  annotate("text", x = 0.7, y = -1.5, label = "Connectors", size = 5) +
  theme_minimal() +  
  labs(x = "Among-module connectivity (Pi)",
       y = "Within-module connectivity (Zi)",
       color = "Functional groups") #+ # To visualize node labels
  geom_text(mapping = aes(x = surface_attribute_node$Pi, y = surface_attribute_node$Zi, label = surface_attribute_node$Label)) 


# Calculate percentage of size classes and trophic modes within each module
size_percentages_surface <- surface_attribute_node %>%
    group_by(modularity, size_class) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(modularity) %>%
    mutate(percentage = round(100 * count / sum(count), 2)) %>%
    select(-count) %>%
    pivot_wider(names_from = size_class, values_from = percentage, values_fill = 0)
  
trophism_percentages_surface <- surface_attribute_node %>%
    group_by(modularity, Trophy) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(modularity) %>%
    mutate(percentage = round(100 * count / sum(count), 2)) %>%
    select(-count) %>%
    pivot_wider(names_from = Trophy, values_from = percentage, values_fill = 0)
  
aggregated_surface <- size_percentages_surface %>%
    left_join(trophism_percentages_surface, by = "modularity")

# Convert size data to long format
size_long_surf <- aggregated_surface %>%
  pivot_longer(cols = meso:mega, names_to = "Size_Category", values_to = "Percentage")

# Convert trophic mode data to long format
trophic_long_surf <- aggregated_surface %>%
  pivot_longer(cols = "autotrophic protist":zooplankton, names_to = "Trophic_mode", values_to = "Percentage")

# Size category plot
plot_size_surf <- ggplot(size_long_surf, aes(x = factor(modularity), y = Percentage, fill = Size_Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Modularity", y = "Percentage", fill = "Size Class")  +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# Trophic mode plot
plot_trophic_surf <- ggplot(trophic_long_surf, aes(x = factor(modularity), y = Percentage, fill = Trophic_mode)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Modularity", y = "Percentage", fill = "Functional Group")  +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Combine plots with cowplot
plot_grid(plot_size_surf, plot_trophic_surf, ncol = 1, align = "v")













