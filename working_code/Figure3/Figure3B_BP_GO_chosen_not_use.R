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

### Check the biological process of carbonhydrate transport
setwd("~/Dropbox/PPiSeq_02/working_data/Positive_PPI_environment/PPI_pair_GO/environment")
DMSO_matrix = csvReader_T("DMSO/Network_density_p_value_matrix_BP.csv")
Forskolin_matrix = csvReader_T("Forskolin/Network_density_p_value_matrix_BP.csv")
FK506_matrix = csvReader_T("FK506/Network_density_p_value_matrix_BP.csv")
NaCl_matrix = csvReader_T("NaCl/Network_density_p_value_matrix_BP.csv")
Raffinose_matrix = csvReader_T("Raffinose/Network_density_p_value_matrix_BP.csv")
HU_matrix = csvReader_T("HU/Network_density_p_value_matrix_BP.csv")
H2O2_matrix = csvReader_T("H2O2/Network_density_p_value_matrix_BP.csv")
Dox_matrix = csvReader_T("Dox/Network_density_p_value_matrix_BP.csv")
cold_matrix = csvReader_T("16C/Network_density_p_value_matrix_BP.csv")

### Write a function to extract specific BP GO term. Meanwhile, filter out some GO terms in the other hand
extract_specific_BP = function(DMSO_matrix, chosen_BP, filter_BP){
        if(length(chosen_BP) == 1){
                chosen_BP_matrix = DMSO_matrix[which(DMSO_matrix[,1] == chosen_BP & DMSO_matrix[,2] %in% filter_BP),]
        }else{
                chosen_BP_matrix = DMSO_matrix[which(DMSO_matrix[,1] %in% chosen_BP & DMSO_matrix[,2] %in% filter_BP),]
        }
        return(chosen_BP_matrix)
}
chosen_BP = c("carbohydrate transport", "DNA templated transcription  initiation", "translational initiation")
filter_BP = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/chosen_BP_GO_59.txt",
                                 header = F, sep = "\t"))
DMSO_chosen = extract_specific_BP(DMSO_matrix, chosen_BP, filter_BP)
Forskolin_chosen = extract_specific_BP(Forskolin_matrix, chosen_BP, filter_BP)
FK506_chosen = extract_specific_BP(FK506_matrix, chosen_BP, filter_BP)
NaCl_chosen = extract_specific_BP(NaCl_matrix, chosen_BP, filter_BP)
Raffinose_chosen = extract_specific_BP(Raffinose_matrix, chosen_BP, filter_BP)
HU_chosen = extract_specific_BP(HU_matrix, chosen_BP, filter_BP)
H2O2_chosen = extract_specific_BP(H2O2_matrix, chosen_BP, filter_BP)
Dox_chosen = extract_specific_BP(Dox_matrix, chosen_BP, filter_BP)
cold_chosen = extract_specific_BP(cold_matrix, chosen_BP, filter_BP)


DMSO_transform = function(DMSO_chosen, chosen_BP, environment_label){
        BP_label = rep("0", nrow(DMSO_chosen))
        BP_label[which(DMSO_chosen[,1] == chosen_BP[1])] = chosen_BP[1]
        DMSO_chosen[which(DMSO_chosen[,1] == chosen_BP[1]),1] = environment_label[1]
        BP_label[which(DMSO_chosen[,1] == chosen_BP[2])] = chosen_BP[2]
        DMSO_chosen[which(DMSO_chosen[,1] == chosen_BP[2]),1] = environment_label[2]
        BP_label[which(DMSO_chosen[,1] == chosen_BP[3])] = chosen_BP[3]
        DMSO_chosen[which(DMSO_chosen[,1] == chosen_BP[3]),1] = environment_label[3]
        DMSO_final = cbind(DMSO_chosen, BP_label)
        return(DMSO_final)
}
DMSO_final = DMSO_transform(DMSO_chosen, chosen_BP, c("SD_transport", "SD_transcription", "SD_translation"))
Forskolin_final = DMSO_transform(Forskolin_chosen, chosen_BP, 
                                 c("Forskolin_transport", "Forskolin_transcription", "Forskolin_translation"))
FK506_final = DMSO_transform(FK506_chosen, chosen_BP, 
                                 c("FK506_transport", "FK506_transcription", "FK506_translation"))
NaCl_final = DMSO_transform(NaCl_chosen, chosen_BP, 
                                 c("NaCl_transport", "NaCl_transcription", "NaCl_translation"))
Raffinose_final = DMSO_transform(Raffinose_chosen, chosen_BP, 
                                 c("Raffinose_transport", "Raffinose_transcription", "Raffinose_translation"))
HU_final = DMSO_transform(HU_chosen, chosen_BP, 
                                 c("Hydroxyurea_transport", "Hydroxyurea_transcription", "Hydroxyurea_translation"))
H2O2_final = DMSO_transform(H2O2_chosen, chosen_BP, 
                                 c("H2O2_transport", "H2O2_transcription", "H2O2_translation"))
Dox_final = DMSO_transform(Dox_chosen, chosen_BP, 
                                 c("Doxorubicin_transport", "Doxorubicin_transcription", "Doxorubicin_translation"))
cold_final = DMSO_transform(cold_chosen, chosen_BP, 
                                 c("16C_transport", "16C_transcription", "16C_translation"))


all_matrix_chosen = rbind(DMSO_final, Forskolin_final, FK506_final, NaCl_final,
                          Raffinose_final, HU_final, H2O2_final, Dox_final, cold_final)
all_frame_chosen = data.frame(all_matrix_chosen[,1], all_matrix_chosen[,2], 
                              as.numeric(all_matrix_chosen[,3]), as.numeric(all_matrix_chosen[,4]),
                              all_matrix_chosen[,5])
colnames(all_frame_chosen) = c("rowv", "columnv", "Network_density", "label", "BP_label")

GO_order = read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                      header = T, sep = "\t")
environment_order = c( "16C_transcription", "Doxorubicin_transcription", "H2O2_transcription", "Hydroxyurea_transcription",
                      "Raffinose_transcription", "NaCl_transcription", "FK506_transcription","Forskolin_transcription",
                      "SD_transcription",
                      "16C_translation", "Doxorubicin_translation","H2O2_translation", "Hydroxyurea_translation",
                      "Raffinose_translation", "NaCl_translation", "FK506_translation",   "Forskolin_translation", 
                      "SD_translation",
                      "16C_transport", "Doxorubicin_transport", "H2O2_transport", "Hydroxyurea_transport",
                      "Raffinose_transport", "NaCl_transport", "FK506_transport",  "Forskolin_transport", "SD_transport")


shapes = c("carbohydrate transport" = 16, "DNA templated transcription  initiation" = 15, 
           "translational initiation" = 17)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = columnv, y = rowv, size =Network_density, color = label, shape = BP_label), all_frame_chosen)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(7,10)], limits = c(0, 0.1),
                              oob = scales::squish)+ 
        scale_y_discrete(limits = environment_order) + 
        scale_x_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = c(0, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11), 
                   label = c("0", "1%", "3%", "5%", "7%", "9%", "11%"), range = c(0,4))+ ## to tune the size of circles
        scale_shape_manual(name = "GO term", values = shapes) +
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1), plot.margin = unit(c(0.2,0.2,0.5,2.5), "cm"),
              axis.text.y.left = element_text(size = 9, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/Figure3B_BP_GO_chosen.pdf", width = 14, height = 9)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure3/Figure3B_BP_GO_chosen.pdf", width = 14, height =9)
