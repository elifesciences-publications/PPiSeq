#### Detect communities and check their distribution of mean variation score
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
vScore_PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")

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
csvWriter(vScore_community_data, "Figure3_related_network_data/community_fast_greedy.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Figure3_related_network_data/community_fast_greedy_summary.csv")

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
csvWriter(vScore_community_data, "Figure3_related_network_data/community_walktrap.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Figure3_related_network_data/community_walktrap_summary.csv")

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
csvWriter(vScore_community_data, "Figure3_related_network_data/community_infomap.csv")

summary_community = aggregate(vScore~community, data=vScore_community_data, 
                              FUN=function(x) c(mean=mean(x), count=length(x)))
summary_community = data.frame(community = summary_community$community,
                               vScore = summary_community$vScore[,1],
                               count = summary_community$vScore[,2])
summary_order = summary_community[order(summary_community[,2], decreasing = T),]
csvWriter(summary_order, "Figure3_related_network_data/PPI_network/community_infomap_summary.csv")

#### Communities detected by Infomap algorithm was used 
#### Do GO enrichment and check which CC, BP, MF are enriched in each community
infomap = csvReader_T("Figure3_related_network_data/community_infomap.csv")
## infoMAP community 1: core, community 3: core-1, community 2: accessory
infomap_stable = infomap[which(as.numeric(infomap[,2]) == 1),]
infomap_median = infomap[which(as.numeric(infomap[,2]) == 3),]
infomap_unstable = infomap[which(as.numeric(infomap[,2]) == 2),]

stable_protein = infomap_stable[,1]
median_protein = infomap_median[,1]
unstable_protein = infomap_unstable[,1]

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
csvWriter(tgo, "Figure3_related_network_data/Stable_CC_enriched.csv")

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
csvWriter(tgo, "Figure3_related_network_data/Stable_BP_enriched.csv")


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
csvWriter(tgo, "Figure3_related_network_data/Stable_MF_enriched.csv")

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
csvWriter(tgo, "Figure3_related_network_data/Median_CC_enriched.csv")

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
csvWriter(tgo, "Figure3_related_network_data/Median_BP_enriched.csv")


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
csvWriter(tgo, "Figure3_related_network_data/Median_MF_enriched.csv")

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
csvWriter(tgo, "Figure3_related_network_data/Unstable_CC_enriched.csv")

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
csvWriter(tgo, "Figure3_related_network_data/Unstable_BP_enriched.csv")


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
csvWriter(tgo, "Figure3_related_network_data/Unstable_MF_enriched.csv")
