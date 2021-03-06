### Extract all PPIs that contain proteins that are associated with carbonhydrate transport
### and then check network in different environments (SD and Raffinose)
setwd("~/Desktop/PPiSeq_additional_data/") # GO:0008643 carbonhydrate transport
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
GO_slim = as.matrix(read.table("Outsourced_datasets/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Gene_Carbon = unique(GO_slim[which(GO_slim[,6] == "GO:0008643"), 1])
PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_pos_SD_merge_filter.csv")
# The fitness values of negative PPIs are zero.
PPI_count = csvReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv")
check_specific_protein = function(PPI, Gene_Carbon){
        PPI_chosen = "0"
        protein_pair = split_string_vector(PPI[,1])
        for(i in 1:nrow(protein_pair)){
                if(protein_pair[i,1] %in% Gene_Carbon | protein_pair[i,2] %in% Gene_Carbon){
                        PPI_chosen = c(PPI_chosen, PPI[i,1])
                }
        }
        PPI_chosen = PPI_chosen[2:length(PPI_chosen)]
        return(PPI_chosen)
}
PPI_carbon = check_specific_protein(PPI, Gene_Carbon)
PPI_carbon_fitness = PPI[which(PPI[,1] %in% PPI_carbon),]
PPI_carbon_count = PPI_count[match(PPI_carbon_fitness[,1], PPI_count[,1]),]
name_exchange = csvReader_T("Outsourced_datasets/Systematic_standard_protein.csv")
PPI_split = split_string_vector(PPI_carbon_fitness[,1])
protein_1 = name_exchange[match(PPI_split[,1], name_exchange[,1]),2]
protein_2 = name_exchange[match(PPI_split[,2], name_exchange[,1]),2]

weight = rep(0, nrow(PPI_carbon_fitness))
for(i in 1:nrow(PPI_carbon_fitness)){
        if(as.numeric(PPI_carbon_count[i,3]) == 1){
                weight[i] = as.numeric(PPI_carbon_fitness[i,10]) 
                # change this number to get network for different environments 
                #DMSO:4, H2O2:5, HU:6, Dox:7, Forskolin:8, Raffinose:9, NaCl:10, 16C:11, FK506:12
        }
}

PPI_net = data.frame(protein_1, protein_2, weight)
group = rep(6, nrow(PPI_net))
col_nodes = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#1f78b4", "#CECED2")
library(scales)
col_nodes_transparent = alpha(col_nodes, 0.4)
col_nodes_transparent[3] = alpha(col_nodes[3], 0.8)
col_nodes_transparent[5] = alpha(col_nodes[5], 0.8)
for(i in 1:nrow(PPI_net)){
        if(PPI_net[i,1] == "HXT1" | PPI_net[i,2] == "HXT1"){
                group[i] = 1
        }else if(PPI_net[i,1] == "HXT3" | PPI_net[i,2] == "HXT3"){
                group[i] = 2
        }else if(PPI_net[i,1] == "HXT5" | PPI_net[i,2] == "HXT5"){
                group[i] = 3
        }else if(PPI_net[i,1] == "HXT7" | PPI_net[i,2] == "HXT7"){
                group[i] = 4
        }else if(PPI_net[i,1] == "HXT2" | PPI_net[i,2] == "HXT2"){
                group[i] = 5
        }
}
PPI_net = data.frame(PPI_net, group)

colnames(PPI_net) = c("Protein01", "Protein02", "weight", "group")
library("igraph")
net = graph_from_data_frame(d = PPI_net, directed = F)
#net = simplify(net,  remove.loops = TRUE)
deg = igraph::degree(net)
deg_order = deg[order(deg, decreasing = T)]
V(net)$size = log2(deg)# Node size stands for the degree of protein
E(net)$color[E(net)$group == 1 ] = col_nodes_transparent[1]
E(net)$color[E(net)$group == 2 ] = col_nodes_transparent[2]
E(net)$color[E(net)$group == 3 ] = col_nodes_transparent[3]
E(net)$color[E(net)$group == 4 ] = col_nodes_transparent[4]
E(net)$color[E(net)$group == 5 ] = col_nodes_transparent[5]
#E(net)$color = apple_colors[10]
V(net)$color = col_nodes[6]
V(net)["HXT1"]$color = col_nodes[1]
V(net)["HXT3"]$color = col_nodes[2]
V(net)["HXT5"]$color = col_nodes[3]
V(net)["HXT7"]$color = col_nodes[4]
V(net)["HXT2"]$color = col_nodes[5]
V(net)$label = NA
E(net)$width <- 2*exp(E(net)$weight)

order_random = sample(order(degree(net)), length(degree(net)), replace = F)
a = rownames(data.frame(degree(net)))[order(degree(net))]
protein = a[order_random]
protein_order = cbind(protein, order_random)
colnames(protein_order) = c("protein", "random_order")
#l = layout_in_circle(net, order=order_random) # Once choose a good layout do not change it. And make network for another environment
net_clean <- delete.edges(net, which(E(net)$weight == 0))

pdf("Figures/Figure5/Figure5C_PPI_carbonhydrate_transport_network_NaCl.pdf", height =5, width = 5)
plot(net_clean, layout = l,  vertex.frame.color=apple_colors[8], vertex.label.color = apple_colors[11],
     vertex.label.cex = 0.35, margin = c(0,0,0,0))
dev.off()


