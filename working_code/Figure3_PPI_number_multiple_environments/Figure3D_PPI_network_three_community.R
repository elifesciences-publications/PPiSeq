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
fastgreedy = infomap
fastgreedy_stable = fastgreedy[which(as.numeric(fastgreedy[,2]) == 1),]
fastgreedy_median = fastgreedy[which(as.numeric(fastgreedy[,2]) == 3),]
fastgreedy_unstable = fastgreedy[which(as.numeric(fastgreedy[,2]) == 2),]
fastgreedy_label = rbind(fastgreedy_stable, fastgreedy_median, fastgreedy_unstable) # 948
fastgreedy_others = fastgreedy[which(!fastgreedy[,1] %in% fastgreedy_label[,1]),] # 1134

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