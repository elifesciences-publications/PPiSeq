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

setwd("~/Dropbox/PPiSeq_02/Working_data_2/PPI_network/")
#fastgreedy = csvReader_T("community_fast_greedy.csv")
#walktrap = csvReader_T("community_walktrap.csv")
infomap = csvReader_T("community_infomap.csv")

infomap_stable = infomap[which(as.numeric(infomap[,2]) == 1),]
infomap_median = infomap[which(as.numeric(infomap[,2]) == 3),]
infomap_unstable = infomap[which(as.numeric(infomap[,2]) == 2),]
infomap_label = rbind(infomap_stable, infomap_median, infomap_unstable) # 948
infomap_others = infomap[which(!infomap[,1] %in% infomap_label[,1]),] # 1134

## Calculate the mean stability score for each community
mean_stable = median(as.numeric(infomap_stable[,3])) # 0.6197287
mean_median = median(as.numeric(infomap_median[,3])) # 0.8063364
mean_unstable = median(as.numeric(infomap_unstable[,3])) # 1.371961

infomap_stable[,2] = "Core (low variability)"
infomap_median[,2] = "Core (intermediate variability)"
infomap_unstable[,2] = "Accessory (high variability)"
community_three = rbind(infomap_stable, infomap_median, infomap_unstable)
colnames(community_three) = c("Protein", "Community", "vScore", "Degree")
csvWriter(community_three, 'Community_label_three.csv')

community_label = dataFrameReader_T("Community_label_three.csv")
community_label$Degree = log10(community_label$Degree)

community_label$Community<- factor(community_label$Community, 
                                   levels = c("Core (low variability)", "Core (intermediate variability)",
                                              "Accessory (high variability)"))
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
        #xlab("Stability score") +
        #ylab(expression('Log'[10]* '(Degree)')) +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank()) +
        #theme(legend.key=element_blank()) +
        #guides(fill=guide_legend(nrow=1,byrow=TRUE))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) 
        #theme(axis.text.x = element_text(size = 10, color = "black"),
              #axis.text.y.left = element_text(size = 10, color = "black"),
              #axis.title.y=element_text(size=10)) + 
        #theme(text = element_text(size=10))

ggsave("~/Dropbox/PPiseq_02/Working_figure/Figure3_accessory_PPIs/Figure3E_three_communities.pdf", width =3, height =2.8)
