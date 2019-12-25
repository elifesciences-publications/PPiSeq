########################### Make a heatmap to show the normalized fitness values for Carbohydrate transport PPIs
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

# Input the normalized fitness values for all PPIs in each environment
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit = csvReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
PPI_count_filter = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")

GO_slim = as.matrix(read.table("Paper_data/Outside_datasets/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))

Gene_carbon = unique(GO_slim[which(GO_slim[,6] == "GO:0008643"), 1])
# A function to extract specific PPIs
check_specific_protein = function(PPI, Gene_Carbon){
        PPI_chosen = "0"
        protein_pair = split_string_vector(PPI[,1])
        if(length(Gene_Carbon) > 1){
                for(i in 1:nrow(protein_pair)){
                        if(protein_pair[i,1] %in% Gene_Carbon | protein_pair[i,2] %in% Gene_Carbon){
                                PPI_chosen = c(PPI_chosen, PPI[i,1])
                        }
                }  
        }else {
                for(i in 1:nrow(protein_pair)){
                        if(protein_pair[i,1] == Gene_Carbon | protein_pair[i,2] == Gene_Carbon){
                                PPI_chosen = c(PPI_chosen, PPI[i,1])
                        }
                }  
        }
        
        PPI_chosen = PPI_chosen[2:length(PPI_chosen)]
        return(PPI_chosen)
}
PPI_carbon = check_specific_protein(PPI_fit, Gene_carbon) 
PPI_carbon_fit = PPI_count_filter[which(PPI_fit[,1] %in% PPI_carbon),]
PPI_carbon_fit_count = as.data.frame(table(PPI_carbon_fit[,2]))

PPI_carbon_fit[which(PPI_carbon_fit[,2] == "1"),]
pdf("Working_figure/SFigures/paper/FigureS9_carbohydrate_transport/Figure_S9A_PPI_stability.pdf", width = 5, height = 5)
par(mar= c(4,5,1,1))
barCenter = barplot(PPI_carbon_fit_count$Freq, ylim = c(0,100), ylab = "Number of PPIs",
                    col = apple_colors[5])
text(x = barCenter, y = -5, labels = PPI_carbon_fit_count$Var1, xpd = TRUE)
text(x = median(barCenter), y = -10, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

## Check what are the GO terms for proteins that interact with these carbonhydrate transport proteins
PPI_carbon_protein = split_string_vector(PPI_carbon)
protein_label = matrix(0, nrow(PPI_carbon_protein), ncol(PPI_carbon_protein))
protein_non_carbon = "0"
for(i in 1:nrow(protein_label)){
        for(j in 1:ncol(protein_label)){
                if (PPI_carbon_protein[i,j] %in% Gene_carbon){
                        protein_label[i,j] = 1
                }else{
                        protein_non_carbon= c(protein_non_carbon, PPI_carbon_protein[i,j])   
                }
        }
}
protein_label = cbind(PPI_carbon, protein_label)
colnames(protein_label) = c("PPI","protein_A_GO_carbon", "Protein_B_GO_carbon")
csvWriter(protein_label, "Working_data_2/Glucose_transporter/PPI_carbo_GO_check.csv")
GO_check = rep(0, nrow(protein_label))
for(i in 1:nrow(protein_label)){
        GO_check[i] = sum(as.numeric(protein_label[i,2:3]))
}
GO_check_table = as.data.frame(table(GO_check))
pdf("Working_figure/SFigures/paper/FigureS9_carbohydrate_transport/Figure_S9B_GO_interacting_partner.pdf", width = 5, height = 5)
par(mar= c(4,5,1,1))
barCenter = barplot(GO_check_table$Freq, ylim = c(0,300), ylab = "Number of interacting partners",
                    col = apple_colors[5])
text(x = barCenter, y = -20, labels = c("No", "Yes"), xpd = TRUE)
text(x = median(barCenter), y = -40, labels = "Invovled in carbohydrate transport", xpd = TRUE)
dev.off()

GO_slim = GO_slim[which(!GO_slim[,5] %in% c("cellular_component", "biological_process","signaling",
                                           "molecular_function", "not_yet_annotated", "other")),]
protein_non_carbon = protein_non_carbon[2:length(protein_non_carbon)]
GO_slim_non_carbon = GO_slim[which(GO_slim[,1] %in% protein_non_carbon),]
GO_slim_non_carbon_BP = GO_slim_non_carbon[which(GO_slim_non_carbon[,4] == "P"),]
GO_slim_non_carbon_CC = GO_slim_non_carbon[which(GO_slim_non_carbon[,4] == "C"),]
GO_slim_non_carbon_MF = GO_slim_non_carbon[which(GO_slim_non_carbon[,4] == "F"),]

GO_slim_non_carbon_BP_table = as.data.frame(table(GO_slim_non_carbon_BP[,5]))
GO_slim_non_carbon_BP_table = GO_slim_non_carbon_BP_table[order(GO_slim_non_carbon_BP_table[,2], decreasing = T),]

pdf("Working_figure/SFigures/paper/FigureS9_carbohydrate_transport/Figure_S9C_interacting_partner_BP_10.pdf", width = 5, height = 4)
par(mar= c(5,12,1,1))
barCenter = barplot(rev(GO_slim_non_carbon_BP_table$Freq[1:10]), xlim = c(0,30), xlab = "Frequency",
                    col = apple_colors[5], horiz = T)
text(x = -15, y = barCenter, labels = rev(GO_slim_non_carbon_BP_table$Var1[1:10]), xpd = TRUE)
#text(x = median(barCenter), y = -40, labels = "Invovled in carbohydrate transport", xpd = TRUE)
dev.off()

GO_slim_non_carbon_CC_table = as.data.frame(table(GO_slim_non_carbon_CC[,5]))
GO_slim_non_carbon_CC_table = GO_slim_non_carbon_CC_table[order(GO_slim_non_carbon_CC_table[,2], decreasing = T),]

pdf("Working_figure/SFigures/paper/FigureS9_carbohydrate_transport/Figure_S9D_interacting_partner_CC_10.pdf", width = 5, height = 4)
par(mar= c(5,10,1,1))
barCenter = barplot(rev(GO_slim_non_carbon_CC_table$Freq[1:10]), xlim = c(0,120), xlab = "Frequency",
                    col = apple_colors[5], horiz = T)
text(x = -45, y = barCenter, labels = rev(GO_slim_non_carbon_CC_table$Var1[1:10]), xpd = TRUE)
#text(x = median(barCenter), y = -40, labels = "Invovled in carbohydrate transport", xpd = TRUE)
dev.off()

GO_slim_non_carbon_MF_table = as.data.frame(table(GO_slim_non_carbon_MF[,5]))
GO_slim_non_carbon_MF_table = GO_slim_non_carbon_MF_table[order(GO_slim_non_carbon_MF_table[,2], decreasing = T),]

pdf("Working_figure/SFigures/paper/FigureS9_carbohydrate_transport/Figure_S9E_interacting_partner_MF_10.pdf", width = 5, height = 4)
par(mar= c(5,14,1,1))
barCenter = barplot(rev(GO_slim_non_carbon_MF_table$Freq[1:10]), xlim = c(0,30), xlab = "Frequency",
                    col = apple_colors[5], horiz = T)
text(x = -20, y = barCenter, labels = rev(GO_slim_non_carbon_MF_table$Var1[1:10]), xpd = TRUE)
#text(x = median(barCenter), y = -40, labels = "Invovled in carbohydrate transport", xpd = TRUE)
dev.off()


