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
all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Molecular function
all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_MF_new.txt"
pos_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_MF_new.txt"
output_pos_count = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_MF_matrix.txt"
output_density = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_MF_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


########### Make a plot to show the enrichment of GO_GO pairs 
# Cellular compartment
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_CC_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column

GO_CC = colnames(network_density)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_CC_new.txt",sep = "\t", header = T))
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

random_name_front = "Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_CC_pos_count_"
density_CC_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_CC_random = cbind(density_CC_random, random_density)
}
density_CC_random = density_CC_random[,2:ncol(density_CC_random)]

### Calculate the z_score based on 1000 random networks
z_score = rep(0, length(network_density_vector))
for(i in 1:length(network_density_vector)){
        mean = mean(as.numeric(na.omit(density_CC_random[i,])))
        sd = sd(as.numeric(na.omit(density_CC_random[i,])))
        z_score[i] = (as.numeric(network_density_vector[i]) - mean)/sd
}

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
        permutation=oneway_test(density_all~density_name, matrix_compare) # permutation test
        p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
}
# dendrogram
library(ggplot2)
library(ggdendro)
rownames(network_density) = colnames(network_density)
hc = hclust(dist(network_density), method="ward.D")
dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 3) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_CC_GO_pair_cluster.pdf", width = 5, height = 5)
label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),]

GO_order = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_CC_order.txt",
                      header = T, sep = "\t")
Network_density = as.numeric(network_density_vector)
dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)
#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.015, 0.03, 0.044), 
                   label =as.character(c(0, 0.015, 0.03, 0.044)), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
       
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = apple_colors[11], angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_CC_PPI_network_p.value.pdf", 
       width = 7, height = 5)

############# Make scatter plot to compare the diagonal densities (the same GO term) with other densities
dimer_index = rep(0, nrow(dataf))
for(i in 1:nrow(dataf)){
        if(dataf[i,1] == dataf[i,2]){
                dimer_index[i] = 1
        }
}
dimer_density = dataf[which(dimer_index == 1),3]
non_dimer_density = dataf[which(dimer_index != 1),3]

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Density_CC_same_different.pdf", width = 5, height = 5)
plot(density(non_dimer_density), col = apple_colors[5], xlab = "Network density", 
     main = NA)
lines(density(dimer_density), col = apple_colors[7], add = T)
legend("topright", c("Different CC", "Same CC"), 
       col = apple_colors[c(5,7)], lty = c(1,1), bty = "n")
dev.off()

density_all = as.numeric(c(dimer_density, non_dimer_density))
name = c(rep("Same cellular compartment", length(dimer_density)), 
         rep("Different cellular compartments", length(non_dimer_density)))
compare_matrix = data.frame(name, density_all)
colnames(compare_matrix) = c("GO_pair", "density")
t.test(dimer_density, non_dimer_density) # 4.756e-05
library(ggplot2)

ggplot(compare_matrix, aes(x = GO_pair, y = density, group = GO_pair, col = GO_pair))+
        geom_violin( draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
        #geom_boxplot(width = 0.1, col ="black", alpha = 0.3)+
        #geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
        stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
        scale_color_manual(name = "", breaks = c('Same cellular compartment', "Different cellular compartments"),
                           values  = apple_colors[c(5,7)])  +
        scale_y_continuous(name = "Network density", 
                           limits=c(0, 0.05),
                           breaks = seq(0,0.05, by =0.01),
                           labels = seq(0,0.05, by= 0.01))+
        theme(legend.position ="none") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"), axis.title.x= element_blank(),
              axis.text.y.left = element_text(size = 10, color = "black")) 
        
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Violin_CC_same_different.pdf", width= 5, height =5)


###### Make a heatmap for positive count
network_pos_count = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_CC_matrix.txt",
                                         sep = "\t", header = T))
colnames(network_pos_count) = gsub("\\.", " ", colnames(network_pos_count))
GO = colnames(network_pos_count)
network_count_vector = as.vector(network_pos_count)
Network_count = as.numeric(network_count_vector)
dataC = data.frame(rowv = rep(GO, each = 22),
                   columnv = rep(GO, 22),
                   Network_count = as.numeric(network_count_vector))
Network_count = as.numeric(network_count_vector)
count_label = c("0", "200", "800", "2000", "4000", "8214")
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Number of PPIs", breaks = c(10, 200, 800, 2000, 4000, 8214), 
                   label = count_label, 
                   range = c(0, 4))+
        theme(legend.position ="right", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_pos_count_CC_PPI_network.pdf", width = 7, height = 5)

########################################
## Biological process
setwd("~/Dropbox/PPiSeq_02")
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 25)
columnv = rep(GO_BP, 25)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_BP_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

random_name_front = "Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_BP_pos_count_"
density_BP_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_BP_random = cbind(density_BP_random, random_density)
}
density_BP_random = density_BP_random[,2:ncol(density_BP_random)]

### Calculate the z_score based on 1000 random networks
z_score = rep(0, length(network_density_vector))
for(i in 1:length(network_density_vector)){
        mean = mean(as.numeric(na.omit(density_BP_random[i,])))
        sd = sd(as.numeric(na.omit(density_BP_random[i,])))
        z_score[i] = (as.numeric(network_density_vector[i]) - mean)/sd
}

### Calculate the p-value based on 1000 random networks
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
                permutation=oneway_test(density_all~density_name, matrix_compare) # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
       
}

# dendrogram
library(ggdendro)
rownames(network_density) = colnames(network_density)
hc = hclust(dist(network_density), method="ward.D")

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 6) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
#ggsave("~/Desktop/Dendrogram_BP_GO_pair_cluster.pdf", width = 6, height = 10)
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_BP_GO_pair_cluster.pdf", width = 8, height = 6)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),] # order the BP GO term by the cluster analysis

Network_density = as.numeric(network_density_vector) # 0.11764
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
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.03, 0.05, 0.07), 
                   label =as.character(c(0, 0.01, 0.03, 0.05, 0.07)), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 6, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_BP_PPI_network.pdf", 
       width = 8.5, height = 5)

###### Make a heatmap for positive count
network_pos_count = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_BP_matrix.txt", header=T, sep = "\t"))
colnames(network_pos_count) = gsub("\\.", " ", colnames(network_pos_count))
GO = colnames(network_pos_count) # 25
network_count_vector = as.vector(network_pos_count)

dataC = data.frame(rowv = rep(GO, each = 25),
                   columnv = rep(GO, 25),
                   Network_count = as.numeric(network_count_vector)
)
Network_count = as.numeric(network_count_vector) # 2970
count_label = c("0", "500", "1000", "1500", "2000", "2500")
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Number of PPIs", breaks = c(0,500,1000,1500,2000,2500), 
                   label = count_label, 
                   range = c(0, 4))+
        theme(legend.position ="right", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 6, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_pos_count_BP_PPI_network.pdf", 
       width = 8.5, height = 5)

########################################
## Molecular function
setwd("~/Dropbox/PPiSeq_02")
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_MF_matrix.txt",
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column
GO_MF = colnames(network_density)
rowv = rep(GO_MF, each = 41)
columnv = rep(GO_MF, 41)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_MF_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
#all_count_name[which(!all_count_name %in% GO_GO_name_order)][1]
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

random_name_front = "Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_MF_pos_count_"
density_MF_random = rep(0, length(GO_GO_name_order))
for (i in 1:1000){
        random_name = paste(random_name_front, as.character(i), ".txt", sep = "")
        random_density = extract_random_density(order_all_count, random_name, GO_GO_name_order)
        density_MF_random = cbind(density_MF_random, random_density)
}
density_MF_random = density_MF_random[,2:ncol(density_MF_random)]

### Calculate the z_score based on 1000 random networks
#z_score = rep(0, length(network_density_vector))
#for(i in 1:length(network_density_vector)){
        #mean = mean(as.numeric(na.omit(density_BP_random[i,])))
        #sd = sd(as.numeric(na.omit(density_BP_random[i,])))
        #z_score[i] = (as.numeric(network_density_vector[i]) - mean)/sd
#}

### Calculate the p-value based on 1000 random networks
p_value = rep(0, length(network_density_vector)) 
for(i in 1:length(network_density_vector)){
        m_real = as.numeric(network_density_vector[i])
        m_random = as.numeric(density_MF_random[i,])
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
                permutation=oneway_test(density_all~density_name, matrix_compare) # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
        
}


rownames(network_density) = colnames(network_density)
# dendrogram
library(ggdendro)
#rownames(network_density) = as.character(1:98)
hc = hclust(dist(network_density), method="ward.D")

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 5) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
#ggsave("~/Desktop/Dendrogram_MF_GO_pair_cluster.pdf", width = 6, height = 6)
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_MF_GO_pair_cluster.pdf", width = 5, height = 5)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),] # order the BP GO term by the cluster analysis

GO_order = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_MF_order.txt",
                      header = T, sep = "\t")

dataf = data.frame(rowv = rep(GO_MF, each = 41),
                   columnv = rep(GO_MF, 41),
                   Network_density = as.numeric(network_density_vector),
                   label = p_value)
Network_density = as.numeric(network_density_vector)
max(Network_density) # 0.1222807
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12), 
                   label =as.character(c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12)), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = apple_colors[11], angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_MF_PPI_network.pdf", 
       width = 11, height = 7)

############# Make scatter plot to compare the diagonal densities (the same GO term) with other densities
dimer_index = rep(0, nrow(dataf))
for(i in 1:nrow(dataf)){
        if(dataf[i,1] == dataf[i,2]){
                dimer_index[i] = 1
        }
}
dimer_density = dataf[which(dimer_index == 1),3]
non_dimer_density = dataf[which(dimer_index != 1),3]

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Density_MF_same_different.pdf", width = 5, height = 5)
plot(density(non_dimer_density), col = apple_colors[5], xlab = "Network density", 
     main = NA)
lines(density(dimer_density), col = apple_colors[7], add = T)
legend("topright", c("Different molecular functions", "Same molecular funciton"), 
       col = apple_colors[c(5,7)], lty = c(1,1), bty = "n")
dev.off()

density_all = as.numeric(c(dimer_density, non_dimer_density))
name = c(rep("Same molecular funciton", length(dimer_density)), 
         rep("Different molecular functions", length(non_dimer_density)))
compare_matrix = data.frame(name, density_all)
colnames(compare_matrix) = c("GO_pair", "density")
t.test(dimer_density, non_dimer_density) # 2.602e-05
library(ggplot2)

ggplot(compare_matrix, aes(x = GO_pair, y = density, group = GO_pair, col = GO_pair))+
        geom_violin( draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
        #geom_boxplot(width = 0.1, col ="black", alpha = 0.3)+
        #geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
        stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
        scale_color_manual(name = "", breaks = c('Same biological process', "Different biological processes"),
                           values  = apple_colors[c(5,7)])  +
        scale_y_continuous(name = "Network density", 
                           limits=c(0, 0.12),
                           breaks = seq(0,0.12, by =0.03),
                           labels = seq(0,0.12, by= 0.03))+
        theme(legend.position ="none") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 11, color = "black"), axis.title.x= element_blank(),
              axis.text.y.left = element_text(size = 11, color = "black")) +
        
        ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Violin_MF_same_different.pdf", width= 5, height =5)

###### Make a heatmap for positive count
network_pos_count = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_pos_count_PPI_MF_matrix.txt", header=T, sep = "\t"))
colnames(network_pos_count) = gsub("\\.", " ", colnames(network_pos_count))
GO = colnames(network_pos_count) # 41
network_count_vector = as.vector(network_pos_count)

dataC = data.frame(rowv = rep(GO, each = 41),
                   columnv = rep(GO, 41),
                   Network_count = as.numeric(network_count_vector)
)
Network_count = as.numeric(network_count_vector) # 393
count_label = c("0", "80", "160", "240", "320", "393")
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Number of PPIs", breaks = c(0,80,160,240,320,393), 
                   label = count_label, 
                   range = c(0, 4))+
        theme(legend.position ="right", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_pos_count_MF_PPI_network.pdf", width = 11, height = 7)


################################ Make plots with the original BP GO terms
setwd("~/Dropbox/PPiSeq_02/")
# Biological process
all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 99)
columnv = rep(GO_BP, 99)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network

#### extract density from each random network and put them as a column into a matrix
all_count= as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt",sep = "\t", header = T))
all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
all_count_name = gsub("-", " ", all_count_name)
all_count_name = gsub(",", " ", all_count_name)
order_all_count = all_count[match(GO_GO_name_order, all_count_name), 3]
#overlap = intersect(all_count_name, GO_GO_name_order)

random_name_front = "Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/BP_primary/random_network_BP_pos_count_"
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
                permutation=oneway_test(density_all~density_name, matrix_compare) # permutation test
                p_value[i] =  permutation@distribution@pvalue(permutation@statistic@teststatistic)
        }
        
}

# dendrogram
library(ggdendro)
rownames(network_density) = colnames(network_density)
hc = hclust(dist(network_density), method="ward.D")

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 10) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
#ggsave("~/Desktop/Dendrogram_BP_GO_pair_cluster.pdf", width = 6, height = 10)
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_BP_GO_pair_cluster_primary.pdf", width = 6, height = 10)

GO_order = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_BP_order.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) # 0.11764
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
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.03, 0.05, 0.07, 0.09, 0.117), 
                   label =as.character(c(0, 0.01, 0.03, 0.05, 0.07, 0.09, 0.117)), 
                   range = c(0, 3))+ ## to tune the size of circles
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = apple_colors[11], angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_BP_PPI_network_primary.pdf", 
       width = 16, height = 13)

############# Make scatter plot to compare the diagonal densities (the same GO term) with other densities
dimer_index = rep(0, nrow(dataf))
for(i in 1:nrow(dataf)){
        if(dataf[i,1] == dataf[i,2]){
                dimer_index[i] = 1
        }
}
dimer_density = dataf[which(dimer_index == 1),3]
non_dimer_density = dataf[which(dimer_index != 1),3]

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Density_BP_same_different_primary.pdf", width = 5, height = 5)
plot(density(non_dimer_density), col = apple_colors[5], xlab = "Network density", 
     main = NA)
lines(density(dimer_density), col = apple_colors[7], add = T)
legend("topright", c("Different biological processes", "Same biological process"), 
       col = apple_colors[c(5,7)], lty = c(1,1), bty = "n")
dev.off()

density_all = as.numeric(c(dimer_density, non_dimer_density))
name = c(rep("Same biological process", length(dimer_density)), 
         rep("Different biological processes", length(non_dimer_density)))
compare_matrix = data.frame(name, density_all)
colnames(compare_matrix) = c("GO_pair", "density")
t.test(dimer_density, non_dimer_density) # 8.364e-13
library(ggplot2)

ggplot(compare_matrix, aes(x = GO_pair, y = density, group = GO_pair, col = GO_pair))+
        geom_violin( draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
        #geom_boxplot(width = 0.1, col ="black", alpha = 0.3)+
        #geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
        stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
        scale_color_manual(name = "", breaks = c('Same biological process', "Different biological processes"),
                           values  = apple_colors[c(5,7)])  +
        scale_y_continuous(name = "Network density", 
                           limits=c(0, 0.12),
                           breaks = seq(0,0.12, by =0.03),
                           labels = seq(0,0.12, by= 0.03))+
        theme(legend.position ="none") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 11, color = "black"), axis.title.x= element_blank(),
               axis.text.y.left = element_text(size = 11, color = "black")) +

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Violin_BP_same_different_primary.pdf", width= 5, height =5)


###### Make a heatmap for positive count
network_pos_count = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_pos_count_PPI_BP_matrix.txt", header=T, sep = "\t"))
colnames(network_pos_count) = gsub("\\.", " ", colnames(network_pos_count))
GO = colnames(network_pos_count) # 99
network_count_vector = as.vector(network_pos_count)

dataC = data.frame(rowv = rep(GO, each = 99),
                   columnv = rep(GO, 99),
                   Network_count = as.numeric(network_count_vector)
)
Network_count = as.numeric(network_count_vector) # 293
count_label = c("0", "50", "100", "150", "200", "290")
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Number of PPIs", breaks = c(0,50,100,150,200,290), 
                   label = count_label, 
                   range = c(0, 4))+
        theme(legend.position ="right", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_pos_count_BP_PPI_network_primary.pdf", width = 16, height = 13)

