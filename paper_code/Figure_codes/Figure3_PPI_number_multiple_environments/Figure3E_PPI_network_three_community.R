setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
infomap = csvReader_T("Figure3_related_network_data/community_infomap.csv")
infomap_stable = infomap[which(as.numeric(infomap[,2]) == 1),]
infomap_median = infomap[which(as.numeric(infomap[,2]) == 3),]
infomap_unstable = infomap[which(as.numeric(infomap[,2]) == 2),]
infomap_label = rbind(infomap_stable, infomap_median, infomap_unstable) # 948
infomap_others = infomap[which(!infomap[,1] %in% infomap_label[,1]),] # 1134

## Calculate the mean stability score for each community
mean_stable = median(as.numeric(infomap_stable[,3])) # 0.6197287
mean_median = median(as.numeric(infomap_median[,3])) # 0.8063364
mean_unstable = median(as.numeric(infomap_unstable[,3])) # 1.371961

infomap_stable[,2] = "Core (low mutability)"
infomap_median[,2] = "Core (intermediate mutability)"
infomap_unstable[,2] = "Accessory (high mutability)"
community_three = rbind(infomap_stable, infomap_median, infomap_unstable)
colnames(community_three) = c("Protein", "Community", "vScore", "Degree")
csvWriter(community_three, 'Figure3_related_network_data/Figure3E_Community_mutability_degree.csv')

community_label = dataFrameReader_T("Figure3_related_network_data/Figure3E_Community_mutability_degree.csv")
community_label$Degree = log10(community_label$Degree)

community_label$Community<- factor(community_label$Community, 
                                   levels = c("Core (low mutability)", "Core (intermediate mutability)",
                                              "Accessory (high mutability)"))
#community_label$vScore = 1/community_label$vScore
library(ggplot2)
col_chosen = c("red", "springgreen4", "blue")
ggplot(community_label, aes(x = vScore, y = Degree, fill = Community, 
                            color = Community, shape = Community))+
        geom_point(alpha = 0.5) +
        geom_vline(xintercept = mean_stable, col = "red", linetype = "dashed") +
        geom_vline(xintercept = mean_median, col = "springgreen4", linetype = "dashed") +
        geom_vline(xintercept = mean_unstable, col = "blue", linetype = "dashed") +
        scale_color_manual(values = c(col_chosen)) +
        scale_fill_manual(values = c(col_chosen)) +
        scale_shape_manual(values = c(16,17,16))+
        xlab("Mutability score") +
        ylab(expression('Log'[10]* '(Degree)')) +
        theme(legend.position = c(0.75, 0.9),
              legend.key=element_blank()) +
        #guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) 
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Figures/Figure3/Figure3E_three_communities.pdf", width =5, height =5)
