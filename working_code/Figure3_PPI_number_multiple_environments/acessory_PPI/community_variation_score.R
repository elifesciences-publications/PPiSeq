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

#### Detect communities and check their distribution of mean variation score
setwd("~/Dropbox/PPiSeq_02/")
vScore_PPI = csvReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")

##### Calculate the vScore for each protein
PPI_pair = split_string_vector(vScore_PPI[,1])
protein_unique = unique(as.vector(PPI_pair))
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
        index_protein = unique(c(which(PPI_pair[,1] == protein_unique[i]),
                                 which(PPI_pair[,2] == protein_unique[i])))
        vScore_protein[i,2] = mean(as.numeric(vScore_PPI[index_protein,3]))
}

## Input a network and find its neighbours for a specific node
PPI_network = data.frame(protein_A = PPI_pair[,1], protein_B = PPI_pair[,2], 
                         Env = as.numeric(vScore_PPI[,3]), vScore = as.numeric(vScore_PPI[,4]))
library(igraph)
g = graph_from_data_frame(PPI_network, directed = FALSE)



### Fast greedy algorithm
fast_greedy = fastgreedy.community(g)
max(fast_greedy$modularity) # 0.3948459
fast_greedy_community = cbind(fast_greedy$names, fast_greedy$membership)

vSocre_community = as.numeric(vScore_protein[match(fast_greedy_community[,1], vScore_protein[,1]),2])
vScore_community_data = data.frame(protein = fast_greedy_community[,1], 
                                   community = fast_greedy_community[,2],
                                   vScore = vSocre_community)
csvWriter(vScore_community_data, "Working_data_2/PPI_network/community_fast_greedy.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Working_data_2/PPI_network/community_fast_greedy_summary.csv")

summary_order_select = summary_order[which(summary_order[,3] > 5),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]

vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                levels = c("8","6", "10", "9", "4", "7", "1", "5", "3", "2"))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.015, alpha=0.3,dotsize = 0.8,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(limits = c("8","6", "10", "9", "4", "7", "1", "5", "3", "2"))+
        xlab("Community")+
        ylab("Stability Score") +
       
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/fast_greedy_community_stability_score.pdf", width =4, height =4)

### Random walk trap
walktrap = walktrap.community(g)
walktrap_community = cbind(walk_trap$names, walk_trap$membership)
max(walk_trap$modularity) # 0.3575106
vSocre_community = as.numeric(vScore_protein[match(walktrap_community[,1], vScore_protein[,1]),2])
vScore_community_data = data.frame(protein = walktrap_community[,1], 
                                   community = walktrap_community[,2],
                                   vScore = vSocre_community)
csvWriter(vScore_community_data, "Working_data_2/PPI_network/community_walktrap.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Working_data_2/PPI_network/community_walktrap_summary.csv")

summary_order_select = summary_order[which(summary_order[,3] > 5),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]


#col_chosen = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                                levels = c("10","9", "24", "6", "17", "1", "12", "3", "21", "4"))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.015, alpha=0.3,dotsize = 0.8,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(limits = c("10","9", "24", "6", "17", "1", "12", "3", "21", "4"))+
        xlab("Community")+
        ylab("Stability Score") +
        
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/walktrap_community_stability_score.pdf", width =4, height =4)




#### Label.propagation
label = label.propagation.community(g)
max(label$modularity) #0.0004620597,  0.0004620864, 0.3242257,0.0003080685
#### edge.betweenness
edge_between = edge.betweenness.community(g, directed = FALSE)

### spinglass.community, Cannot work with unconnected graph
spinglass = spinglass.community(g)
walktrap_community = cbind(label$names, label$membership)

## leading.eigenvector.community
eigenvector = leading.eigenvector.community(g)
eigenvector_community = cbind(eigenvector$names, eigenvector$membership)
max(eigenvector$modularity) #0.3461259
vSocre_community = as.numeric(vScore_protein[match(eigenvector_community[,1], vScore_protein[,1]),2])
vScore_community_data = data.frame(protein = eigenvector_community[,1], 
                                   community = eigenvector_community[,2],
                                   vScore = vSocre_community)
csvWriter(vScore_community_data, "Working_data_2/PPI_network/community_eigenvector.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Working_data_2/PPI_network/community_eigenvector_summary.csv")

summary_order_select = summary_order[which(summary_order[,3] > 5),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]


#col_chosen = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                                levels = c("7","8", "3", "4", "1", "5"))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.015, alpha=0.3,dotsize = 0.8,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(limits = c("7","8", "3", "4", "1", "5"))+
        xlab("Community")+
        ylab("Stability Score") +
        
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/eigenvector_community_stability_score.pdf", width =4, height =4)

### informap
infomap = infomap.community(g)
infomap$modularity #0.382528
infomap_community = cbind(infomap$names, infomap$membership)

vSocre_community = as.numeric(vScore_protein[match(infomap_community[,1], vScore_protein[,1]),2])
vScore_community_data = data.frame(protein = infomap_community[,1], 
                                   community = infomap_community[,2],
                                   vScore = vSocre_community)
csvWriter(vScore_community_data, "Working_data_2/PPI_network/community_infomap.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Working_data_2/PPI_network/community_infomap_summary.csv")

summary_order_select = summary_order[which(summary_order[,3] > 10),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]


#col_chosen = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                                levels = as.character(summary_order_select$community))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.015, alpha=0.3,dotsize = 0.8,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(limits = as.character(summary_order_select$community))+
        xlab("Community")+
        ylab("Stability Score") +
        
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/Infomap_community_stability_score.pdf", width =6, height =4)

