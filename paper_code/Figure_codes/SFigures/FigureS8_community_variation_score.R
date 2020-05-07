setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
### Fast greedy algorithm
vScore_community_data = dataFrameReader_T("Figure3_related_network_data/community_fast_greedy.csv")
summary_order = dataFrameReader_T("Figure3_related_network_data/community_fast_greedy_summary.csv")
summary_order_select = summary_order[which(summary_order[,3] > 10),]
vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]
vScore_community_data_select$degree= log10(vScore_community_data_select$degree)
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                levels = as.character(summary_order_select$community))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.015, alpha=0.5,dotsize = 1,  
                     fill = apple_colors[5], col = apple_colors[5]) +
        scale_x_discrete(name = "Community",limits = as.character(summary_order_select$community))+
        #scale_y_continuous(name = expression('Log'[2]* '(Degree)'))+
        scale_y_continuous(name = "Average variability score of PPIs") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Figures/SFigures/SFigure8/FigureS8A_fast_greedy_community_variability_score.pdf", width =4, height =4)

### Random walk trap
vScore_community_data = dataFrameReader_T("Figure3_related_network_data/community_walktrap.csv")
summary_order = dataFrameReader_T("Figure3_related_network_data/community_walktrap_summary.csv")
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
        scale_y_continuous(name = "Average variability score of PPIs") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Figures/SFigures/SFigure8//FigureS8B_walktrap_community_variability_score.pdf", width =4, height =4)

### informap
vScore_community_data = dataFrameReader_T("Figure3_related_network_data/community_infomap.csv")
summary_order = dataFrameReader_T("Figure3_related_network_data/community_infomap_summary.csv")
summary_order_select = summary_order[which(summary_order[,3] > 10),]

vScore_community_data_select = vScore_community_data[which(as.character(vScore_community_data$community) %in% 
                                                                   as.character(summary_order_select[,1])),]
vScore_community_data_select$degree = log10(vScore_community_data_select$degree)
vScore_community_data_select$community = factor(vScore_community_data_select$community , 
                                                levels = as.character(summary_order_select$community))
library(ggplot2)
ggplot(vScore_community_data_select, aes(x = community, y = vScore))+
        geom_boxplot(outlier.shape=NA)+
        geom_dotplot(binaxis="y",stackdir="center",binwidth=0.016, alpha=0.5,dotsize = 1,  
                     fill = apple_colors[5], col = apple_colors[5])+
        scale_x_discrete(name = "Community",limits = as.character(summary_order_select$community))+
        #scale_y_continuous(name = expression('Log'[2]* '(Degree)'))+
        scale_y_continuous(name = "Average variability score of PPIs") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Figures/SFigures/SFigure8/SFigure8C_Infomap_community_variability_score.pdf", width =5, height =4)

