###########################
# Source some basic functions froma function.R in Github repository
source_https <- function(u, unlink.tmp.certs = FALSE) {
        # load package
        require(RCurl)
        # read script lines from website using a security certificate
        if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
        script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
        if(unlink.tmp.certs) unlink("cacert.pem")
        
        # parase lines and evealuate in the global environement
        eval(parse(text = script), envir= .GlobalEnv)
}
source_https("https://raw.githubusercontent.com/sashaflevy/PPiSeq/master/working_code/function.R", unlink.tmp.certs = TRUE)

#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

### Extract all PPIs that contain proteins that are associated with carbonhydrate transport
### and then check network in different environments (SD and Raffinose)
setwd("~/Dropbox/PPiSeq_02/") # GO:0008643 carbonhydrate transport
GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Gene_Carbon = unique(GO_slim[which(GO_slim[,6] == "GO:0008643"), 1])
PPI = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_pos.csv")

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
PPI_carbon_fitness = PPI_carbon_fitness[which(PPI_carbon_fitness[,1] != "YHR096C_YHR096C"),]
name_exchange = csvReader_T("Working_data/Systematic_standard_protein.csv")
PPI_split = split_string_vector(PPI_carbon_fitness[,1])
protein_1 = name_exchange[match(PPI_split[,1], name_exchange[,1]),2]
protein_2 = name_exchange[match(PPI_split[,2], name_exchange[,1]),2]
weight = as.numeric(PPI_carbon_fitness[,5])
# DMSO:4, H2O2:5, HU:6, Dox:7, Forskolin:8, Raffinose:9, NaCl:10, 16C:11, FK506:12
#PPI_net = data.frame(protein_1, protein_2, weight, label)
#PPI_net = PPI_net[which(PPI_net$weight!= 0),]
PPI_net = data.frame(protein_1, protein_2, weight)
group = rep(6, nrow(PPI_net))
#col_nodes = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#CECED2")

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
#E(net)$color[E(net)$reported == 1 ] = apple_colors[1]
#E(net)$color[E(net)$reported == 2 ] = apple_colors[4]


#order_random = sample(order(degree(net)), 144, replace = F)
#a = rownames(data.frame(degree(net)))[order(degree(net))]
#protein = a[order_random]
#protein_order = cbind(protein, order_random)
#colnames(protein_order) = c("protein", "random_order")
#csvWriter(protein_order, "Working_figure/Figure4/Figure4B_dynamic_PPI_network/carbohydrate_transport/Index_order_network_circle.csv")
#l = layout_in_circle(net, order=order_random)
net_clean <- delete.edges(net, which(E(net)$weight == 0))

pdf("Working_figure/Figure4/Figure4B_dynamic_PPI_network/carbohydrate_transport/PPI_carbonhydrate_transport_H2O2.pdf", height =5, width = 5)
#pdf("~/Desktop/Figure4/PPI_carbonhydrate_transport_Fk506_circle.pdf", height =5, width = 5)
plot(net_clean, layout = l,  vertex.frame.color=apple_colors[8], vertex.label.color = apple_colors[11],
     vertex.label.cex = 0.35, margin = c(0,0,0,0))
dev.off()


