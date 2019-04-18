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

### Figure S2B check the network density for different groups of PPIs
setwd("~/Dropbox/PPiSeq_02/")
PPI_env_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_1 = PPI_env_count[which(PPI_env_count[,2] == "1"),1]
PPI_2 = PPI_env_count[which(PPI_env_count[,2] == "2"),1]
PPI_3 = PPI_env_count[which(PPI_env_count[,2] == "3"),1]
PPI_4 = PPI_env_count[which(PPI_env_count[,2] == "4"),1]
PPI_5 = PPI_env_count[which(PPI_env_count[,2] == "5"),1]
PPI_6 = PPI_env_count[which(PPI_env_count[,2] == "6"),1]
PPI_7 = PPI_env_count[which(PPI_env_count[,2] == "7"),1]
PPI_8 = PPI_env_count[which(PPI_env_count[,2] == "8"),1]
PPI_9 = PPI_env_count[which(PPI_env_count[,2] == "9"),1]
network_density = function(PPI_1){
        PP_pair = split_string_vector(PPI_1)
        protein_unique = length(unique(c(PP_pair[,1], PP_pair[,2]))) # 2089
        potential_PPI = protein_unique * (protein_unique-1)/2
        actual_PPI = nrow(PP_pair)
        density = actual_PPI/potential_PPI
        return(density)
}
protein_degree_count = function(PPI){
        all_PPI_gene = split_string_vector(PPI)
        protein_degree = as.data.frame(table(as.character(c(all_PPI_gene[,1], all_PPI_gene[,2]))))
        protein_degree_order= protein_degree[order(protein_degree[,2], decreasing = T),]
        return(protein_degree_order)
}
mean_degree = function(PPI_1){
        PPI_degree = protein_degree_count(PPI_1)
        degree_mean = mean(PPI_degree[,2])
        return(degree_mean)
}
density = network_density(PPI_env_count[,1]) # 0.006209406
degree = mean_degree(PPI_env_count[,1]) # 13.10185


density_1 = network_density(PPI_1) # 0.003834169
density_2 = network_density(PPI_2) # 0.004776649
density_2
density_3 = network_density(PPI_3) # 0.006150755
density_3
density_4 = network_density(PPI_4) # 0.006875281
density_4
density_5 = network_density(PPI_5) # 0.007168598
density_5
density_6 = network_density(PPI_6) # 0.006317131
density_6
density_7 = network_density(PPI_7) # 0.007514337
density_7
density_8 = network_density(PPI_8) # 0.007914589
density_8
density_9 = network_density(PPI_9) # 0.01865794
density_9

degree_1 = mean_degree(PPI_1) # 8.005744
degree_1
degree_2 = mean_degree(PPI_2) # 3.596817
degree_2
degree_3 = mean_degree(PPI_3) # 2.835498
degree_3
degree_4 = mean_degree(PPI_4) # 2.708861
degree_4
degree_5 = mean_degree(PPI_5) # 2.817259
degree_5
degree_6 = mean_degree(PPI_6) # 2.893246
degree_6
degree_7 = mean_degree(PPI_7) # 2.923077
degree_7
degree_8 = mean_degree(PPI_8) # 3.126263
degree_8
degree_9 = mean_degree(PPI_9) # 4.365957
degree_9


# Only check these reported PPIs
reported_PPI = csvReader_T("Working_data/multiple_validated_PPI.csv")
PPI_env_count_reported = match_both_direction(PPI_env_count, reported_PPI[,1]) # 1853
PPI_1_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "1"),1]
PPI_2_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "2"),1]
PPI_3_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "3"),1]
PPI_4_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "4"),1]
PPI_5_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "5"),1]
PPI_6_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "6"),1]
PPI_7_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "7"),1]
PPI_8_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "8"),1]
PPI_9_reported = PPI_env_count_reported[which(PPI_env_count_reported[,2] == "9"),1]

density_1_r = network_density(PPI_1_reported) # 0.005344439
density_1_r
density_2_r = network_density(PPI_2_reported) # 0.00802005
density_2_r
density_3_r = network_density(PPI_3_reported) # 0.009145844
density_3_r
density_4_r = network_density(PPI_4_reported) # 0.00829666
density_4_r
density_5_r = network_density(PPI_5_reported) # 0.008790687
density_5_r
density_6_r = network_density(PPI_6_reported) # 0.008790687
density_6_r
density_7_r = network_density(PPI_7_reported) # 0.007967878
density_7_r
density_8_r = network_density(PPI_8_reported) # 0.007486157
density_8_r
density_9_r = network_density(PPI_9_reported) # 0.01723602
density_9_r

degree_1_r = mean_degree(PPI_1_reported) # 1.811765
degree_1_r
degree_2_r = mean_degree(PPI_2_reported) # 1.515789
degree_2_r
degree_3_r = mean_degree(PPI_3_reported) # 1.426752
degree_3_r
degree_4_r = mean_degree(PPI_4_reported) # 1.377246
degree_4_r
degree_5_r = mean_degree(PPI_5_reported) # 1.608696
degree_5_r
degree_6_r = mean_degree(PPI_6_reported) # 1.643411
degree_6_r
degree_7_r = mean_degree(PPI_7_reported) # 2.007905
degree_7_r
degree_8_r = mean_degree(PPI_8_reported) # 2.245847
degree_8_r
degree_9_r = mean_degree(PPI_9_reported) # 2.757764
degree_9_r

density_all = c(density_1, density_2, density_3, density_4, density_5, density_6,
            density_7, density_8, density_9)
density_reported = c(density_1_r, density_2_r, density_3_r, density_4_r, 
                     density_5_r, density_6_r, density_7_r, density_8_r, density_9_r)

### Network density of all PPIs in different number of environments
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2B_network_density_mean_degree/FigureS2B_network_density_bar_plot.pdf", width= 6, height=5)
barCenter = barplot(density_all*100, horiz=F, beside=F, ylim = c(0,2), ylab="Network density (%)",
                    space= c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                    col= apple_colors[5] , axisnames=F, border=NA)
text(x= barCenter, y = -0.15, labels = as.character(1:9), xpd = TRUE)
text(x = mean(barCenter), y = -0.3, labels = "Number of environments in which a PPI is observed", xpd = TRUE)
dev.off()
### Network density of reported PPIs in different number of environments
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2B_network_density_mean_degree/FigureS2B_network_density_bar_plot_reported.pdf", width= 6, height=5)
barCenter = barplot(density_reported*100, horiz=F, beside=F, ylim = c(0,2), ylab="Network density (%)",
                    space= c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                    col= apple_colors[5] , axisnames=F, border=NA)
text(x= barCenter, y = -0.15, labels = as.character(1:9), xpd = TRUE)
text(x = mean(barCenter), y = -0.3, labels = "Number of environments in which a PPI is observed", xpd = TRUE)
dev.off()

mean_degree_all = c(degree_1, degree_2, degree_3, degree_4, degree_5, degree_6,
                    degree_7, degree_8, degree_9)
mean_degree_reported = c(degree_1_r, degree_2_r, degree_3_r, degree_4_r, degree_5_r,
                         degree_6_r, degree_7_r, degree_8_r, degree_9_r)
### Mean degree of all PPIs in different number of environments
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2B_network_density_mean_degree/FigureS2B_mean_degree_bar_plot.pdf", width= 6, height=5)
barCenter = barplot(mean_degree_all, horiz=F, beside=F, ylim = c(0,10), ylab="Mean degree",
                    space= c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                    col= apple_colors[5] , axisnames=F, border=NA)
text(x= barCenter, y = -1, labels = as.character(1:9), xpd = TRUE)
text(x = mean(barCenter), y = -2, labels = "Number of environments in which a PPI is observed", xpd = TRUE)
dev.off()

### Mean degree of reported PPIs in different number of environments
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2B_network_density_mean_degree/FigureS2B_mean_degree_bar_plot_reported.pdf", width= 6, height=5)
barCenter = barplot(mean_degree_reported, horiz=F, beside=F, ylim = c(0,4), ylab="Mean degree",
                    space= c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                    col= apple_colors[5] , axisnames=F, border=NA)
text(x= barCenter, y = -0.3, labels = as.character(1:9), xpd = TRUE)
text(x = mean(barCenter), y = -0.6, labels = "Number of environments in which a PPI is observed", xpd = TRUE)
dev.off()


