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

##### Check the variation score for each GO-GO pair across different environments
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/")
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
## 16C
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "16C/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "16C/Network_pos_count_PPI_CC_matrix.txt"
output_density = "16C/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "16C/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "16C/Network_pos_count_PPI_BP_matrix.txt"
output_density = "16C/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

##Dox
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Dox/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Dox/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Dox/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Dox/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Dox/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Dox/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## FK506
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "FK506/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "FK506/Network_pos_count_PPI_CC_matrix.txt"
output_density = "FK506/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "FK506/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "FK506/Network_pos_count_PPI_BP_matrix.txt"
output_density = "FK506/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## Forskolin
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Forskolin/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Forskolin/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Forskolin/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Forskolin/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Forskolin/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Forskolin/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## H2O2
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "H2O2/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "H2O2/Network_pos_count_PPI_CC_matrix.txt"
output_density = "H2O2/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "H2O2/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "H2O2/Network_pos_count_PPI_BP_matrix.txt"
output_density = "H2O2/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## HU
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "HU/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "HU/Network_pos_count_PPI_CC_matrix.txt"
output_density = "HU/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "HU/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "HU/Network_pos_count_PPI_BP_matrix.txt"
output_density = "HU/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## NaCl
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "NaCl/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "NaCl/Network_pos_count_PPI_CC_matrix.txt"
output_density = "NaCl/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "NaCl/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "NaCl/Network_pos_count_PPI_BP_matrix.txt"
output_density = "NaCl/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## Raffinose
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Raffinose/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Raffinose/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Raffinose/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Raffinose/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Raffinose/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Raffinose/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


##############################################################################
## Calculate the variation score of the network density for each GO_GO pair across different environments
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/")

# Cellular compartment
sd_CC_density = as.matrix(read.table("DMSO/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
cold_CC_density = as.matrix(read.table("16C/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
Dox_CC_density = as.matrix(read.table("Dox/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
FK506_CC_density = as.matrix(read.table("FK506/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
Forskolin_CC_density = as.matrix(read.table("Forskolin/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
H2O2_CC_density = as.matrix(read.table("H2O2/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
HU_CC_density = as.matrix(read.table("HU/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
NaCl_CC_density = as.matrix(read.table("NaCl/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
Raffinose_CC_density = as.matrix(read.table("Raffinose/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))

cc_list = list(sd_CC_density,cold_CC_density,Dox_CC_density,FK506_CC_density, Forskolin_CC_density,
               H2O2_CC_density,HU_CC_density,NaCl_CC_density,Raffinose_CC_density)

cc_list_mean = apply(simplify2array(cc_list), 1:2, mean)
cc_list_sd = apply(simplify2array(cc_list), 1:2, sd)

cc_list_CV = cc_list_sd/cc_list_mean

colnames(cc_list_CV) = gsub("\\.", " ", colnames(cc_list_CV))
GO_CC = colnames(cc_list_CV)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)

network_density_vector = as.vector(cc_list_CV) # by column
CV = as.numeric(network_density_vector)
hist(CV)
max(as.numeric(network_density_vector), na.rm = T) # 3
min(as.numeric(network_density_vector), na.rm = T) # 0.180
GO_order = read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_CC_order.txt",header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector))
csvWriter(dataf, "Variation_CC_all_environments.csv")

library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = Network_density), dataf)  + 
        scale_color_gradientn(name = "Coefficient variation", colors = apple_colors[c(10,5)], limits = c(0, 3),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Coefficient variation", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3), 
                   label = c("0", "0.1" , "0.2", "0.3", "0.4", "0.5", "1", "2", "3"), range = c(0,3))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS4_PPI_enrichment_GO/heatmap_density_CC_PPI_network_variation.pdf", width = 8, height =4.5)

# Biological process
sd_BP_density = as.matrix(read.table("DMSO/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
cold_BP_density = as.matrix(read.table("16C/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
Dox_BP_density = as.matrix(read.table("Dox/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
FK506_BP_density = as.matrix(read.table("FK506/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
Forskolin_BP_density = as.matrix(read.table("Forskolin/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
H2O2_BP_density = as.matrix(read.table("H2O2/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
HU_BP_density = as.matrix(read.table("HU/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
NaCl_BP_density = as.matrix(read.table("NaCl/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
Raffinose_BP_density = as.matrix(read.table("Raffinose/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))

BP_list = list(sd_BP_density,cold_BP_density,Dox_BP_density,FK506_BP_density, Forskolin_BP_density,
               H2O2_BP_density,HU_BP_density,NaCl_BP_density,Raffinose_BP_density)

BP_list_mean = apply(simplify2array(BP_list), 1:2, mean)
BP_list_sd = apply(simplify2array(BP_list), 1:2, sd)

BP_list_CV = BP_list_sd/BP_list_mean

colnames(BP_list_CV) = gsub("\\.", " ", colnames(BP_list_CV))
chosen_BP = as.matrix(read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/chosen_BP_GO_59.txt",
                                 header = F, sep = "\t"))
all_BP = colnames(BP_list_CV)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = BP_list_CV[row_column_good,row_column_good] # 59


network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 59)
columnv = rep(GO_BP, 59)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

max(as.numeric(network_density_vector), na.rm = T) # 3
min(as.numeric(network_density_vector), na.rm = T) # 0

GO_order = read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                      header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Network_density = as.numeric(network_density_vector))
csvWriter(dataf, "Variation_BP_all_environments.csv")
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = Network_density), dataf)  + 
        scale_color_gradientn(name = "Coefficient variation", colors = apple_colors[c(10,5)], limits = c(0, 3),
                              oob = scales::squish)+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Coefficient variation", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3), 
                   label = c("0", "0.1" , "0.2", "0.3", "0.4", "0.5", "1", "2", "3"), range = c(0,3))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS4_PPI_enrichment_GO/heatmap_density_BP_PPI_network_variation.pdf", width = 12, height =7.5)

#### Take the mean CV for each GO term across all GO terms
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/")
CC_variation = dataFrameReader_T("Variation_CC_all_environments.csv")
BP_variation = dataFrameReader_T("Variation_BP_all_environments.csv")
CC_unique = unique(CC_variation$rowv)
BP_unique = unique(BP_variation$rowv)
CC_CV = rep(0, length(CC_unique))
BP_CV = rep(0, length(BP_unique))
for(i in 1:length(CC_unique)){
        index = which(CC_variation$rowv == CC_unique[i])
        CC_CV[i] = mean(CC_variation[index, 3])
}
for(i in 1:length(BP_unique)){
        index = which(BP_variation$rowv == BP_unique[i])
        BP_CV[i] = mean(BP_variation[index, 3])
}

matrix_CC = cbind(as.character(CC_unique), CC_CV)
matrix_BP = cbind(as.character(BP_unique), BP_CV)
colnames(matrix_CC) = c("CC", "CV")
colnames(matrix_BP) = c("BP", "CV")
matrix_CC_order = matrix_CC[order(matrix_CC[,2], decreasing = T),]
matrix_BP_order = matrix_BP[order(matrix_BP[,2], decreasing = T),]
csvWriter(matrix_CC_order, "Variation_single_CC_order.csv")
csvWriter(matrix_BP_order, "Variation_single_BP_order.csv")

### plot the top 10 most dynamic GO terms
pdf("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS4_PPI_enrichment_GO/Top10_dynamic_CCs.pdf", width= 5.5, height=5)
par(mar = c(7,4,1,1))
barCenter = barplot(as.numeric(matrix_CC_order[1:10,2]), horiz=F, beside=F, ylim=c(0,1.5), 
                    ylab="Coefficient variration",
                    #space= c(0.15, 0.15,  0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
                    axisnames=F, border=NA, col = apple_colors[5], cex.axis=0.8)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x = barCenter, y = -0.05, labels = as.character(matrix_CC_order[1:10,1]), 
     adj = 1, srt =45, xpd = TRUE)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

pdf("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS4_PPI_enrichment_GO/Top10_dynamic_BPs.pdf", width= 6, height=6)
par(mar = c(10,5,1,1))
barCenter = barplot(as.numeric(matrix_BP_order[1:10,2]), horiz=F, beside=F, ylim=c(0,1.5), 
                    ylab="Coefficient variration",
                    #space= c(0.15, 0.15,  0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
                    axisnames=F, border=NA, col = apple_colors[5], cex.axis=0.7)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x = barCenter, y = -0.05, labels = as.character(matrix_BP_order[1:10,1]), 
     adj = 1, srt =45, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()