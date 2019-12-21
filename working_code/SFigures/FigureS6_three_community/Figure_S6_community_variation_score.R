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
        vScore_protein[i,2] = mean(1/(as.numeric(vScore_PPI[index_protein,3])))
}

#### Calculate the degree for each protein
protein_degree_count = function(PPI){
        all_PPI_gene = split_string_vector(PPI)
        protein_degree = as.data.frame(table(as.character(c(all_PPI_gene[,1], all_PPI_gene[,2]))))
        protein_degree_order= protein_degree[order(protein_degree[,2], decreasing = T),]
        return(protein_degree_order)
}

degree_protein = protein_degree_count(vScore_PPI[,1])

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
degree_community = as.numeric(degree_protein[match(fast_greedy_community[,1], degree_protein[,1]),2])
vScore_community_data = data.frame(protein = fast_greedy_community[,1], 
                                   community = fast_greedy_community[,2],
                                   vScore = vSocre_community,
                                   degree = degree_community)
csvWriter(vScore_community_data, "Working_data_2/PPI_network/community_fast_greedy.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Working_data_2/PPI_network/community_fast_greedy_summary.csv")

summary_order_select = summary_order[which(summary_order[,3] > 10),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]
vScore_community_data_select$degree= log10(vScore_community_data_select$degree)
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                levels = as.character(summary_order_select$community))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.018, alpha=0.5,dotsize = 1,  
                     fill = apple_colors[5], col = apple_colors[5]) +
        scale_x_discrete(name = "Community",limits = as.character(summary_order_select$community))+
        #scale_y_continuous(name = expression('Log'[2]* '(Degree)'))+
        scale_y_continuous(name = "Stability score") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/fast_greedy_community_stability_score.pdf", width =4, height =4)

### Random walk trap
walktrap = walktrap.community(g)
walktrap_community = cbind(walktrap$names, walktrap$membership)
max(walktrap$modularity) # 0.3575106
vSocre_community = as.numeric(vScore_protein[match(walktrap_community[,1], vScore_protein[,1]),2])
degree_community = as.numeric(degree_protein[match(walktrap_community[,1], degree_protein[,1]),2])
vScore_community_data = data.frame(protein = walktrap_community[,1], 
                                   community = walktrap_community[,2],
                                   vScore = vSocre_community,
                                   degree = degree_community)
csvWriter(vScore_community_data, "Working_data_2/PPI_network/community_walktrap.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Working_data_2/PPI_network/community_walktrap_summary.csv")

summary_order_select = summary_order[which(summary_order[,3] > 10),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]

vScore_community_data_select$degree= log10(vScore_community_data_select$degree)
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                                levels = as.character(summary_order_select$community))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.025, alpha=0.5,dotsize = 1,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(name = "Community", limits = as.character(summary_order_select$community))+
        #scale_y_continuous(name = expression('Log'[2]* '(Degree)'))+
        scale_y_continuous(name = "Scalability score") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/walktrap_community_stability_score.pdf", width =4, height =4)

### informap
infomap = infomap.community(g)
infomap$modularity # 0.384403

infomap_community = cbind(infomap$names, infomap$membership)

vSocre_community = as.numeric(vScore_protein[match(infomap_community[,1], vScore_protein[,1]),2])
degree_community = as.numeric(degree_protein[match(infomap_community[,1], degree_protein[,1]),2])

vScore_community_data = data.frame(protein = infomap_community[,1], 
                                   community = infomap_community[,2],
                                   vScore = vSocre_community,
                                   degree = degree_community)
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
vScore_community_data_select$degree = log10(vScore_community_data_select$degree)
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                                levels = as.character(summary_order_select$community))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.02, alpha=0.5,dotsize = 1,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(name = "Community",limits = as.character(summary_order_select$community))+
        #scale_y_continuous(name = expression('Log'[2]* '(Degree)'))+
        scale_y_continuous(name = "Stability score") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/Infomap_community_stability_score.pdf", width =5, height =4)


#### Fastgreedy providing the maximum modularity and therefore I use this algorithm
setwd("~/Dropbox/PPiSeq_02/Working_data_2/PPI_network/")
#fastgreedy = csvReader_T("community_fast_greedy.csv")
#walktrap = csvReader_T("community_walktrap.csv")
infomap = csvReader_T("community_infomap.csv")
fastgreedy = infomap
'''
fastgreedy_stable = fastgreedy[which(as.numeric(fastgreedy[,2]) == 2),1]
walktrap_stable = walktrap[which(as.numeric(walktrap[,2]) == 4),1]
infomap_stable = infomap[which(as.numeric(infomap[,2]) == 1),1]

fastgreedy_median = fastgreedy[which(as.numeric(fastgreedy[,2]) == 3),1]
walktrap_median = walktrap[which(as.numeric(walktrap[,2]) == 3),1]
infomap_median = infomap[which(as.numeric(infomap[,2]) == 3),1]

fastgreedy_unstable = fastgreedy[which(as.numeric(fastgreedy[,2]) == 1),1]
walktrap_unstable = walktrap[which(as.numeric(walktrap[,2]) == 1),1]
infomap_unstable = infomap[which(as.numeric(infomap[,2]) == 2),1]

length(fastgreedy_stable) # 357
length(walktrap_stable) # 332
length(infomap_stable) # 305
overlap_stable = intersect(intersect(fastgreedy_stable, walktrap_stable), infomap_stable)
length(overlap_stable) # 275
length(intersect(fastgreedy_stable, walktrap_stable)) # 295

length(fastgreedy_median) # 425
length(walktrap_median) # 555
length(infomap_median) # 147
overlap_median = intersect(intersect(fastgreedy_median, walktrap_median), infomap_median)
length(overlap_median) # 110
length(intersect(fastgreedy_median, walktrap_median)) # 257

length(fastgreedy_unstable) # 984
length(walktrap_unstable) # 825
length(infomap_unstable) # 463
overlap_unstable = intersect(intersect(fastgreedy_unstable, walktrap_unstable), infomap_unstable)
length(overlap_unstable) # 425
length(intersect(fastgreedy_unstable, walktrap_unstable)) # 693

all_community_label = c(overlap_stable, overlap_median, overlap_unstable) #810

overlap_stable_matrix = fastgreedy[which(fastgreedy[,1] %in% overlap_stable),]
overlap_stable_matrix[,2] = "1"#"Stable"

overlap_median_matrix = fastgreedy[which(fastgreedy[,1] %in% overlap_median),]
overlap_median_matrix[,2] = "2"#Median"

overlap_unstable_matrix = fastgreedy[which(fastgreedy[,1] %in% overlap_unstable),]
overlap_unstable_matrix[,2] = "3"#"Unstable"

non_overlap_matrix = fastgreedy[which(!fastgreedy[,1] %in% all_community_label),]
non_overlap_matrix[,2] = "4" #"Others"
'''
## Label the protein with overlaped communities
## fastgreedy community 1: accessory, community 3: Core-1, community 2: Core
## infoMAP community 1: core, community 3: core-1, community 2: accessory
fastgreedy_stable = fastgreedy[which(as.numeric(fastgreedy[,2]) == 1),]
fastgreedy_median = fastgreedy[which(as.numeric(fastgreedy[,2]) == 3),]
fastgreedy_unstable = fastgreedy[which(as.numeric(fastgreedy[,2]) == 2),]
fastgreedy_label = rbind(fastgreedy_stable, fastgreedy_median, fastgreedy_unstable) # 948
fastgreedy_others = fastgreedy[which(!fastgreedy[,1] %in% fastgreedy_label[,1]),] # 1134

fastgreedy_stable[,2] = "1"
fastgreedy_median[,2] = "2"
fastgreedy_unstable[,2] = "3"
fastgreedy_others[,2] = "4"

community_label_matrix = rbind(fastgreedy_stable, fastgreedy_median,
                               fastgreedy_unstable, fastgreedy_others)

csvWriter(community_label_matrix, "Community_label_number.csv")

fastgreedy_stable[,2] = "Core"
fastgreedy_median[,2] = "Median"
fastgreedy_unstable[,2] = "Accessory"
community_three = rbind(fastgreedy_stable, fastgreedy_median, fastgreedy_unstable)
csvWriter(community_three, 'Community_label_three.csv')

community_label = dataFrameReader_T("Community_label_three.csv")
community_label$degree = log10(community_label$degree)

library(ggplot2)
ggplot(community_label, aes(x = vScore, y = degree, fill = community, color = community))+
        geom_point(pch = 16, alpha = 0.5) +
        scale_color_manual(values = c(apple_colors[c(5,7,3)], apple_colors[9])) +
        scale_fill_manual(values = c(apple_colors[c(5,7,3)], apple_colors[9])) +
        xlab("Stability score") +
        ylab(expression('Log'[10]* '(Degree)')) +
        theme(legend.key=element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("~/Dropbox/PPiseq_02/Working_figure/Figure3_accessory_PPIs/accessory_PPI/Three_communities.pdf", width =4, height =3)

stable_protein = fastgreedy_stable[,1]
median_protein = fastgreedy_median[,1]
unstable_protein = fastgreedy_unstable[,1]


### Do GO enrichment and check which CC, BP, MF are enriched in each community
library("org.Sc.sgd.db")
library("GOstats")

params <- new("GOHyperGParams",
              geneIds=stable_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="CC",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Stable_CC_enriched.csv")

params <- new("GOHyperGParams",
              geneIds=stable_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="BP",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Stable_BP_enriched.csv")


params <- new("GOHyperGParams",
              geneIds=stable_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="MF",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Stable_MF_enriched.csv")

###### median proteins enrichment

params <- new("GOHyperGParams",
              geneIds=median_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="CC",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Median_CC_enriched.csv")

params <- new("GOHyperGParams",
              geneIds=median_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="BP",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Median_BP_enriched.csv")


params <- new("GOHyperGParams",
              geneIds=median_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="MF",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Median_MF_enriched.csv")

###### Unstable proteins enrichment

params <- new("GOHyperGParams",
              geneIds=unstable_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="CC",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Unstable_CC_enriched.csv")

params <- new("GOHyperGParams",
              geneIds=unstable_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="BP",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Unstable_BP_enriched.csv")


params <- new("GOHyperGParams",
              geneIds=unstable_protein,
              universeGeneIds=protein_unique,
              annotation="org.Sc.sgd.db",
              ontology="MF",
              pvalueCutoff=0.05,
              conditional=FALSE,
              testDirection="over")
hypGO= hyperGTest(params)
tgo = summary(hypGO)
csvWriter(tgo, "Unstable_MF_enriched.csv")




###################################################################################

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

