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

###### Make a heatmap for PPIs detected in all environments
## Make a heatmap to show the network density for each pair of GO terms
setwd("~/Dropbox/PPiSeq_02")
### Write a function to produce matrix for pos_PPI_count and pos_PPI_density
transform_list_matrix = function(all_count_name, pos_count_name, output_pos_count, output_density){
        all_count= as.matrix(read.table(all_count_name,sep = "\t", header = T))
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        pos_match_all_count = pos_count[match(all_count_name, pos_count_name),3]
        pos_match_all_count[is.na(pos_match_all_count)] = 0
        density_all = as.numeric(pos_match_all_count)/as.numeric(all_count[,3])
        network_density = cbind(all_count, pos_match_all_count, density_all)
        colnames(network_density) = c("GO_1", "GO_2", "all_count", "pos_count", "density")
        GO_unique = unique(network_density[,1])
        
        matrix_pos_count = matrix(0, length(GO_unique), length(GO_unique))
        matrix_density = matrix(0, length(GO_unique), length(GO_unique))
        for(i in 1:length(GO_unique)){
                a = GO_unique[i]
                for(j in 1:length(GO_unique)){
                        b = GO_unique[j]
                        if(length(which(network_density[,1] == a & network_density[,2] == b)) != 0){
                                matrix_pos_count[i,j] = network_density[which(network_density[,1] == a & network_density[,2] == b),4]
                                matrix_density[i,j] = network_density[which(network_density[,1] == a & network_density[,2] == b),5]  
                        }
                        else{
                                matrix_pos_count[i,j] = 0
                                matrix_density[i,j] = 0
                        }
                        
                }
        }
        colnames(matrix_pos_count) = GO_unique
        colnames(matrix_density) = GO_unique
        write.table(matrix_pos_count, output_pos_count, sep = "\t", quote = F)
        write.table(matrix_density, output_density, sep = "\t", quote=F)
}
# Cellular compartment
all_count_name = "Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Working_data_2/PPI_pair_GO/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Working_data_2/PPI_pair_GO/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Working_data_2/PPI_pair_GO/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Working_data_2/PPI_pair_GO/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Working_data_2/PPI_pair_GO/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Working_data_2/PPI_pair_GO/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Working_data_2/PPI_pair_GO/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

# Molecular function
all_count_name = "Working_data_2/PPI_pair_GO/Network_all_count_PPI_MF_new.txt"
pos_count_name = "Working_data_2/PPI_pair_GO/Network_pos_count_PPI_MF_new.txt"
output_pos_count = "Working_data_2/PPI_pair_GO/Network_pos_count_PPI_MF_matrix.txt"
output_density = "Working_data_2/PPI_pair_GO/Network_density_PPI_MF_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


########### Make a plot to show the enrichment of GO_GO pairs 
setwd("~/Dropbox/PPiSeq_02")
# Cellular compartment
network_density = as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_density_PPI_CC_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column

GO_CC = colnames(network_density)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]

extract_random_density = function(order_all_count, pos_count_name, GO_GO_name_order){
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        order_pos_count = pos_count[match(GO_GO_name_order, pos_count_name),3]
        order_all_count[is.na(order_all_count)] = 0
        order_pos_count[is.na(order_pos_count)] = 0
        density_order = rep(0, length(GO_GO_name_order))
        for(i in 1:length(density_order)){
                if(as.numeric(order_pos_count[i]) != 0 & as.numeric(order_all_count[i]) != 0){
                        density_order[i] = as.numeric(order_pos_count[i])/as.numeric(order_all_count[i])
                }
                else{
                        density_order[i] = 0
                }
        }
        return(density_order)
}

random_name_front = "Working_data_2/PPI_pair_GO/random_network/random_network_density/random_network_CC_pos_count_"
density_CC_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_CC_random = cbind(density_CC_random, random_density)
}
density_CC_random = density_CC_random[,2:ncol(density_CC_random)]

### Calculate the p-value based on 1000 random networks
### Permutation test
library(coin)
p_value = rep(0, length(network_density_vector))
for(i in 1:length(network_density_vector)){
        m_real = network_density_vector[i]
        m_random = as.numeric(density_CC_random[i,])
        density_all = c(m_real, m_random)
        density_name = c("real", rep("random", 1000))
        matrix_compare = cbind.data.frame(density_name, density_all)
        permutation=oneway_test(density_all~density_name, matrix_compare, alternative = "less") # permutation test
        p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
}

GO_order = read.table("Working_data_2/PPI_pair_GO/environment/GO_CC_order.txt",header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) 
max(Network_density) #  0.05169792
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0,  0.02, 0.04,  0.06, 0.08, 0.1), 
                   label = c("0", "2%", "4%", "6%", "8%", "10%"), range = c(0,1.396))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black",angle = -30), axis.title = element_blank())
ggsave("Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_density_CC_PPI_network.pdf", width = 8, height =4.5)

# Only plot the heatmap so that I can cut diagonally
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0,  0.02, 0.04,  0.06, 0.08, 0.1), 
                   label = c("0", "2%", "4%", "6%", "8%", "10%"),
                   range = c(0, 1.396))+ ## to tune the size of circles
        theme(legend.position= "none", plot.margin = margin(0.2, 0.2, 2, 2, "cm")) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_blank(), axis.title = element_blank())
ggsave("Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_density_CC_PPI_network_no_label.pdf", 
       width = 4, height = 4)


#### Biological process
setwd("~/Dropbox/PPiSeq_02/")
# (1) Generate a full picture for all BP GO terms
network_density = as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 99)
columnv = rep(GO_BP, 99)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_all_count_PPI_BP_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

#random_name_front = "Working_data/Positive_PPI_environment/PPI_pair_GO/environment/DMSO/random_network/random_network_density/random_network_BP_pos_count_"
random_name_front = "Working_data_2/PPI_pair_GO/random_network/random_network_density/random_network_BP_pos_count_"
density_BP_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_BP_random = cbind(density_BP_random, random_density)
}
density_BP_random = density_BP_random[,2:ncol(density_BP_random)]

### Calculate the p-value based on 1000 random networks
library(coin)
p_value = rep(0, length(network_density_vector)) 
for(i in 1:length(network_density_vector)){
        m_real = as.numeric(network_density_vector[i])
        m_random = as.numeric(density_BP_random[i,])
        if(m_real == 0){
                p_value[i] = 1
        }
        else if(sum(m_random) == 0){
                p_value[i] = 0
        }
        else{
                density_all = c(m_real, m_random)
                density_name = c("real", rep("random", 1000))
                matrix_compare = cbind.data.frame(density_name, density_all)
                permutation=oneway_test(density_all~density_name, matrix_compare, alternative = "less") # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
        
}

GO_order = read.table("Working_data_2/PPI_pair_GO/environment/GO_BP_order_primary.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) 
max(Network_density) # 0.1111111
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)

library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0,  0.02, 0.04,  0.06, 0.08, 0.1), 
                   label =c("0", "2%" , "4%", "6%", "8%", "10%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_text(size = 9, color = "black"), axis.title = element_blank())

ggsave("Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_density_BP_PPI_network_primary.pdf", 
       width = 16, height = 12)

## No label so that I can cut the figure
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0,  0.02, 0.04,  0.06, 0.08, 0.1), 
                   label =c("0", "2%" , "4%", "6%", "8%", "10%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "none") +
  
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_blank(), axis.title = element_blank())
ggsave("Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_density_BP_PPI_network_primary_no_label.pdf", 
       width = 12, height = 12)

## (2) Choose only 59 GO terms to make figures
setwd("~/Dropbox/PPiseq_02/")
network_density = as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
chosen_BP = as.matrix(read.table("Working_data_2/PPI_pair_GO/chosen_BP_GO_59.txt",
                                 header = F, sep = "\t"))
all_BP = colnames(network_density)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = network_density[row_column_good,row_column_good] # 59

network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 59)
columnv = rep(GO_BP, 59)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_all_count_PPI_BP_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

random_name_front = "Working_data_2/PPI_pair_GO/random_network/random_network_density/random_network_BP_pos_count_"
density_BP_random = rep(0, length(GO_GO_name_order))
extract_random_density = function(order_all_count, pos_count_name, GO_GO_name_order){
        
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        order_pos_count = pos_count[match(GO_GO_name_order, pos_count_name),3]
        order_all_count[is.na(order_all_count)] = 0
        order_pos_count[is.na(order_pos_count)] = 0
        density_order = rep(0, length(GO_GO_name_order))
        for(i in 1:length(density_order)){
                if(as.numeric(order_pos_count[i]) != 0 & as.numeric(order_all_count[i]) != 0){
                        density_order[i] = as.numeric(order_pos_count[i])/as.numeric(order_all_count[i])
                }
                else{
                        density_order[i] = 0
                }
        }
        return(density_order)
}

for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_BP_random = cbind(density_BP_random, random_density)
}
density_BP_random = density_BP_random[,2:ncol(density_BP_random)]

### Calculate the p-value based on 1000 random networks
library(coin)
p_value = rep(0, length(network_density_vector)) 
for(i in 1:length(network_density_vector)){
        m_real = as.numeric(network_density_vector[i])
        m_random = as.numeric(density_BP_random[i,])
        if(m_real == 0){
                p_value[i] = 1
        }
        else if(sum(m_random) == 0){
                p_value[i] = 0
        }
        else{
                density_all = c(m_real, m_random)
                density_name = c("real", rep("random", 1000))
                matrix_compare = cbind.data.frame(density_name, density_all)
                permutation=oneway_test(density_all~density_name, matrix_compare, alternative = "less") # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
        
}

GO_order = read.table("Working_data_2/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) # 0.11111111 (2.149237)folds of network density for cellular compartment)
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)

library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x, position = "right") +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0,  0.02, 0.04,  0.06, 0.08, 0.1), 
                   label =c("0", "2%" , "4%", "6%", "8%", "10%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "left",legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
     
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.right = element_text(size = 9, color = "black"), axis.title = element_blank())

ggsave("Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_density_BP_PPI_network_chosen.pdf", 
       width = 12, height = 7.5)

## No label so that I can cut the figure
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0,  0.02, 0.04,  0.06, 0.08, 0.1), 
                   label =c("0", "2%" , "4%", "6%", "8%", "10%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "none") +
  
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_blank(), axis.title = element_blank())

ggsave("Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_density_BP_PPI_network_chosen_no_label.pdf", 
       width = 7.5, height = 7.5)


###############################################################################################


###### Make a heatmap for SD environment
## Make a heatmap to show the network density for each pair of GO terms
setwd("~/Dropbox/PPiSeq_02")
### Write a function to produce matrix for pos_PPI_count and pos_PPI_density
transform_list_matrix = function(all_count_name, pos_count_name, output_pos_count, output_density){
        all_count= as.matrix(read.table(all_count_name,sep = "\t", header = T))
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        pos_match_all_count = pos_count[match(all_count_name, pos_count_name),3]
        pos_match_all_count[is.na(pos_match_all_count)] = 0
        density_all = as.numeric(pos_match_all_count)/as.numeric(all_count[,3])
        network_density = cbind(all_count, pos_match_all_count, density_all)
        colnames(network_density) = c("GO_1", "GO_2", "all_count", "pos_count", "density")
        GO_unique = unique(network_density[,1])
        
        matrix_pos_count = matrix(0, length(GO_unique), length(GO_unique))
        matrix_density = matrix(0, length(GO_unique), length(GO_unique))
        for(i in 1:length(GO_unique)){
                a = GO_unique[i]
                for(j in 1:length(GO_unique)){
                        b = GO_unique[j]
                        if(length(which(network_density[,1] == a & network_density[,2] == b)) != 0){
                                matrix_pos_count[i,j] = network_density[which(network_density[,1] == a & network_density[,2] == b),4]
                                matrix_density[i,j] = network_density[which(network_density[,1] == a & network_density[,2] == b),5]  
                        }
                        else{
                                matrix_pos_count[i,j] = 0
                                matrix_density[i,j] = 0
                        }
                        
                }
        }
        colnames(matrix_pos_count) = GO_unique
        colnames(matrix_density) = GO_unique
        write.table(matrix_pos_count, output_pos_count, sep = "\t", quote = F)
        write.table(matrix_density, output_density, sep = "\t", quote=F)
}
# Cellular compartment
all_count_name = "Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Working_data_2/PPI_pair_GO/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

# Molecular function
all_count_name = "Working_data_2/PPI_pair_GO/Network_all_count_PPI_MF_new.txt"
pos_count_name = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_MF_new.txt"
output_pos_count = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_MF_matrix.txt"
output_density = "Working_data_2/PPI_pair_GO/environment/DMSO/Network_density_PPI_MF_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


########### Make a plot to show the enrichment of GO_GO pairs 
setwd("~/Dropbox/PPiSeq_02")
# Cellular compartment
network_density = as.matrix(read.table("Working_data_2/PPI_pair_GO/environment/DMSO/Network_density_PPI_CC_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column

GO_CC = colnames(network_density)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]

extract_random_density = function(order_all_count, pos_count_name, GO_GO_name_order){
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        order_pos_count = pos_count[match(GO_GO_name_order, pos_count_name),3]
        order_all_count[is.na(order_all_count)] = 0
        order_pos_count[is.na(order_pos_count)] = 0
        density_order = rep(0, length(GO_GO_name_order))
        for(i in 1:length(density_order)){
                if(as.numeric(order_pos_count[i]) != 0 & as.numeric(order_all_count[i]) != 0){
                        density_order[i] = as.numeric(order_pos_count[i])/as.numeric(order_all_count[i])
                }
                else{
                        density_order[i] = 0
                }
        }
        return(density_order)
}

random_name_front = "Working_data_2/PPI_pair_GO/environment/DMSO/random_network/random_network_density/random_network_CC_pos_count_"
density_CC_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_CC_random = cbind(density_CC_random, random_density)
}
density_CC_random = density_CC_random[,2:ncol(density_CC_random)]

### Calculate the p-value based on 1000 random networks
### Permutation test
library(coin)
p_value = rep(0, length(network_density_vector))
for(i in 1:length(network_density_vector)){
        m_real = network_density_vector[i]
        m_random = as.numeric(density_CC_random[i,])
        density_all = c(m_real, m_random)
        density_name = c("real", rep("random", 1000))
        matrix_compare = cbind.data.frame(density_name, density_all)
        permutation=oneway_test(density_all~density_name, matrix_compare, alternative = "less") # permutation test
        p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
}

GO_order = read.table("Working_data_2/PPI_pair_GO/environment/GO_CC_order.txt",header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) 
max(Network_density) # 0.02793417
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), 
                   label = c("0", "1%" , "2%", "3%", "4%", "5%", "6%", "7%"), range = c(0,1.257))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black",angle = -30), axis.title = element_blank())
ggsave("Working_figure/Figure2_PPI_enrichment_GO/heatmap_density_CC_PPI_network.pdf", width = 8, height =4.5)

# Only plot the heatmap so that I can cut diagonally
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), 
                   label = c("0", "1%" , "2%", "3%", "4%", "5%", "6%", "7%"),
                   range = c(0, 1.257))+ ## to tune the size of circles
        theme(legend.position= "none", plot.margin = margin(0.2, 0.2, 2, 2, "cm")) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_blank(), axis.title = element_blank())
ggsave("Working_figure/Figure2_PPI_enrichment_GO/heatmap_density_CC_PPI_network_no_label.pdf", 
       width = 4, height = 4)


#### Biological process
setwd("~/Dropbox/PPiSeq_02/")
# (1) Generate a full picture for all BP GO terms
network_density = as.matrix(read.table("Working_data_2/PPI_pair_GO/environment/DMSO/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 99)
columnv = rep(GO_BP, 99)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data_2/PPI_pair_GO/Network_all_count_PPI_BP_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

#random_name_front = "Working_data/Positive_PPI_environment/PPI_pair_GO/environment/DMSO/random_network/random_network_density/random_network_BP_pos_count_"
random_name_front = "Working_data_2/PPI_pair_GO/environment/DMSO/random_network/random_network_density/random_network_BP_pos_count_"
density_BP_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_BP_random = cbind(density_BP_random, random_density)
}
density_BP_random = density_BP_random[,2:ncol(density_BP_random)]

### Calculate the p-value based on 1000 random networks
library(coin)
p_value = rep(0, length(network_density_vector)) 
for(i in 1:length(network_density_vector)){
        m_real = as.numeric(network_density_vector[i])
        m_random = as.numeric(density_BP_random[i,])
        if(m_real == 0){
                p_value[i] = 1
        }
        else if(sum(m_random) == 0){
                p_value[i] = 0
        }
        else{
                density_all = c(m_real, m_random)
                density_name = c("real", rep("random", 1000))
                matrix_compare = cbind.data.frame(density_name, density_all)
                permutation=oneway_test(density_all~density_name, matrix_compare, alternative = "less") # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
        
}

GO_order = read.table("Working_data_2/PPI_pair_GO/environment/GO_BP_order_primary.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) # 0.06666667
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)

#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), 
                   label =c("0", "1%" , "2%", "3%", "4%", "5%", "6%", "7%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_text(size = 9, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("Working_figure/Figure2_PPI_enrichment_GO/heatmap_density_BP_PPI_network_primary.pdf", 
       width = 16, height = 12)

## No label so that I can cut the figure
ggplot() + 
geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), 
                   label =c("0", "1%" , "2%", "3%", "4%", "5%", "6%", "7%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "none") +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_blank(), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("Working_figure/Figure2_PPI_enrichment_GO/heatmap_density_BP_PPI_network_primary_no_label.pdf", 
       width = 12, height = 12)

## (2) Choose only 59 GO terms to make figures
setwd("~/Dropbox/PPiseq_02/")
network_density = as.matrix(read.table("Working_data_2/PPI_pair_GO/environment/DMSO/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
chosen_BP = as.matrix(read.table("Working_data_2/PPI_pair_GO/chosen_BP_GO_59.txt",
                       header = F, sep = "\t"))
all_BP = colnames(network_density)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = network_density[row_column_good,row_column_good] # 59

network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 59)
columnv = rep(GO_BP, 59)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

random_name_front = "Working_data_2/PPI_pair_GO/environment/DMSO/random_network/random_network_density/random_network_BP_pos_count_"
density_BP_random = rep(0, length(GO_GO_name_order))
extract_random_density = function(order_all_count, pos_count_name, GO_GO_name_order){
     
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        order_pos_count = pos_count[match(GO_GO_name_order, pos_count_name),3]
        order_all_count[is.na(order_all_count)] = 0
        order_pos_count[is.na(order_pos_count)] = 0
        density_order = rep(0, length(GO_GO_name_order))
        for(i in 1:length(density_order)){
                if(as.numeric(order_pos_count[i]) != 0 & as.numeric(order_all_count[i]) != 0){
                        density_order[i] = as.numeric(order_pos_count[i])/as.numeric(order_all_count[i])
                }
                else{
                        density_order[i] = 0
                }
        }
        return(density_order)
}

for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_BP_random = cbind(density_BP_random, random_density)
}
density_BP_random = density_BP_random[,2:ncol(density_BP_random)]

### Calculate the p-value based on 1000 random networks
library(coin)
p_value = rep(0, length(network_density_vector)) 
for(i in 1:length(network_density_vector)){
        m_real = as.numeric(network_density_vector[i])
        m_random = as.numeric(density_BP_random[i,])
        if(m_real == 0){
                p_value[i] = 1
        }
        else if(sum(m_random) == 0){
                p_value[i] = 0
        }
        else{
                density_all = c(m_real, m_random)
                density_name = c("real", rep("random", 1000))
                matrix_compare = cbind.data.frame(density_name, density_all)
                permutation=oneway_test(density_all~density_name, matrix_compare, alternative = "less") # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
        
}

GO_order = read.table("Working_data_2/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) # 0.0666667 (2.387 folds of network density for cellular compartment)
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)

#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x, position = "right") +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), 
                   label =c("0", "1%" , "2%", "3%", "4%", "5%", "6%", "7%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "left",legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.right = element_text(size = 9, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("Working_figure/Figure2_PPI_enrichment_GO/heatmap_density_BP_PPI_network_chosen.pdf", 
       width = 12, height = 7.5)

## No label so that I can cut the figure
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), 
                   label =c("0", "1%" , "2%", "3%", "4%", "5%", "6%", "7%"), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "none") +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_blank(), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("Working_figure/Figure2_PPI_enrichment_GO/heatmap_density_BP_PPI_network_chosen_no_label.pdf", 
       width = 7.5, height = 7.5)


