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

### Check all gene features enrichment between different groups of PPIs
setwd("~/Dropbox/PPiSeq_02/")
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
PPI_env_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_env_list = vector("list", 9)
for(i in 1:9){
        PPI_chosen = PPI_env_count[which(PPI_env_count[,2] == as.character(i)),1]
        PP_pair = split_string_vector(PPI_chosen)
        PPI_env_list[[i]] = unique(c(PP_pair[,1], PP_pair[,2]))
}
protein_count = rep(0, 9)
for(i in 1:9){
        protein_count[i] = length(PPI_env_list[[i]])
}
protein_count # 2089, 754, 462, 395, 394, 459, 390, 396, 235

features = colnames(gene_feature)
features_chosen = features[2:length(features)]

library(ggplot2)
feature_sort = function(PPI_env_list, gene_feature, specific_feature){
        feature_column = which(features == specific_feature)
        feature_gene = rep(0, 2)
        for(i in 1:9){
                feature = gene_feature[which(gene_feature[,1] %in% PPI_env_list[[i]]),feature_column]
                environment_feature = cbind(rep(i, length(feature)), as.numeric(feature))
                feature_gene = rbind(feature_gene, environment_feature)   
        }
        feature_gene = feature_gene[2:nrow(feature_gene),]
        Environment = as.character(feature_gene[,1])
        Feature = as.numeric(feature_gene[,2])
        feature_plot = data.frame(Environment, Feature)
        plot = ggplot(feature_plot, aes(x = Environment, y = Feature, group = Environment))+
                #geom_boxplot(col = apple_colors[5])+
                geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), col= apple_colors[5])+
                #geom_point(col = apple_colors[3])
                stat_summary(aes(x = Environment, y = Feature, group = Environment), feature_plot,
                             fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
                scale_x_discrete(limits = c("1","2", "3","4", "5", "6", "7", "8", "9"))+
                ylab(specific_feature)+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_text(size = 10, color = "black"),
                      axis.text.y.left = element_text(size = 10, color = "black"))+
                xlab("Number of environments in which a PPI is observed") + 
                ggtitle(specific_feature)
        return(plot)
}
specific_feature = features_chosen[1]
require(gridExtra)
pdf("Working_figure/Figure2/Figure2F_gene_features/Gene_features_across_different_PPIs_1.pdf", width = 15, height = 10)
plot_1= feature_sort(PPI_env_list, gene_feature, features_chosen[1])
plot_2= feature_sort(PPI_env_list, gene_feature, features_chosen[2])
plot_3= feature_sort(PPI_env_list, gene_feature, features_chosen[3])
plot_4= feature_sort(PPI_env_list, gene_feature, features_chosen[4])
plot_5= feature_sort(PPI_env_list, gene_feature, features_chosen[5])
plot_6= feature_sort(PPI_env_list, gene_feature, features_chosen[6])
grid.arrange(plot_1,plot_2, plot_3,plot_4, plot_5,plot_6,nrow= 2)
dev.off()
pdf("Working_figure/Figure2/Figure2F_gene_features/Gene_features_across_different_PPIs_2.pdf", width = 15, height = 10)
plot_7= feature_sort(PPI_env_list, gene_feature, features_chosen[7])
plot_8= feature_sort(PPI_env_list, gene_feature, features_chosen[8])
plot_9= feature_sort(PPI_env_list, gene_feature, features_chosen[9])
plot_10= feature_sort(PPI_env_list, gene_feature, features_chosen[10])
plot_11= feature_sort(PPI_env_list, gene_feature, features_chosen[11])
plot_12= feature_sort(PPI_env_list, gene_feature, features_chosen[12])
grid.arrange(plot_7,plot_8, plot_9,plot_10, plot_11,plot_12,nrow= 2)
dev.off()

pdf("Working_figure/Figure2/Figure2F_gene_features/Gene_features_across_different_PPIs_3.pdf", width = 15, height = 10)
plot_13= feature_sort(PPI_env_list, gene_feature, features_chosen[13])
#plot_14= feature_sort(PPI_env_list, gene_feature, features_chosen[14])
plot_15= feature_sort(PPI_env_list, gene_feature, features_chosen[15])
plot_16= feature_sort(PPI_env_list, gene_feature, features_chosen[16])
plot_17= feature_sort(PPI_env_list, gene_feature, features_chosen[17])
plot_18= feature_sort(PPI_env_list, gene_feature, features_chosen[18])
plot_19= feature_sort(PPI_env_list, gene_feature, features_chosen[19])
grid.arrange(plot_13,  plot_15, plot_16, plot_17,plot_18, plot_19,nrow= 2)
dev.off()

pdf("Working_figure/Figure2/Figure2F_gene_features/Gene_features_across_different_PPIs_4.pdf", width = 15, height = 10)
plot_20= feature_sort(PPI_env_list, gene_feature, features_chosen[20])
plot_21= feature_sort(PPI_env_list, gene_feature, features_chosen[21])
plot_22= feature_sort(PPI_env_list, gene_feature, features_chosen[22])
plot_23= feature_sort(PPI_env_list, gene_feature, features_chosen[23])
plot_24= feature_sort(PPI_env_list, gene_feature, features_chosen[24])
plot_25= feature_sort(PPI_env_list, gene_feature, features_chosen[25])
grid.arrange(plot_20,  plot_21, plot_22, plot_23,plot_24, plot_25,nrow= 2)
dev.off()

pdf("Working_figure/Figure2/Figure2F_gene_features/Gene_features_across_different_PPIs_5.pdf", width = 15, height = 10)
plot_26= feature_sort(PPI_env_list, gene_feature, features_chosen[26])
plot_27= feature_sort(PPI_env_list, gene_feature, features_chosen[27])
plot_28= feature_sort(PPI_env_list, gene_feature, features_chosen[28])
plot_29= feature_sort(PPI_env_list, gene_feature, features_chosen[29])
plot_30= feature_sort(PPI_env_list, gene_feature, features_chosen[30])
plot_31= feature_sort(PPI_env_list, gene_feature, features_chosen[31])
grid.arrange(plot_26,  plot_27, plot_28, plot_29,plot_30, plot_31,nrow= 2)
dev.off()

pdf("Working_figure/Figure2/Figure2F_gene_features/Gene_features_across_different_PPIs_6.pdf", width = 15, height = 15)
plot_32= feature_sort(PPI_env_list, gene_feature, features_chosen[32])
plot_33= feature_sort(PPI_env_list, gene_feature, features_chosen[33])
plot_34= feature_sort(PPI_env_list, gene_feature, features_chosen[34])
plot_35= feature_sort(PPI_env_list, gene_feature, features_chosen[35])
plot_36= feature_sort(PPI_env_list, gene_feature, features_chosen[36])
plot_37= feature_sort(PPI_env_list, gene_feature, features_chosen[37])
plot_38 = feature_sort(PPI_env_list, gene_feature, features_chosen[38])
grid.arrange(plot_32,  plot_33, plot_34, plot_35, plot_36, plot_37, plot_38, nrow= 3)
dev.off()