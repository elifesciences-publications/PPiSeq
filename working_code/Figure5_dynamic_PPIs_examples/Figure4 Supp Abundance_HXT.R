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
protein_abundance = csvReader_T("~/Dropbox/PPiseq_02/Working_data/protein_abundance/table_S4.csv")
protein_median_abundance = log2(as.numeric(protein_abundance[,5]))
HXT = c("HXT1", "HXT2", "HXT3", "HXT4", "HXT5", "HXT6", "HXT7")
protein_abundance_chosen = protein_abundance[which(protein_abundance[,2] %in% HXT),]
method = ncol(protein_abundance_chosen) -6
index_good = 0
for(i in 7:ncol(protein_abundance_chosen)){
       if(sum(is.na(protein_abundance_chosen[,i])) <= 3){
              index_good = c(index_good, i) 
       }
}
protein_chosen = protein_abundance_chosen[,index_good]
mean = rep(0,7)
median = rep(0,7)
for(i in 1:nrow(protein_chosen)){
        mean[i] = mean(na.omit(as.numeric(protein_chosen[i,])))
        median[i] = median(na.omit(as.numeric(protein_chosen[i,])))
}
#mean = apply(protein_chosen,2, )
mean # 1594.800  2755.000  1820.000  8880.667  4548.000  5360.333 10823.750
median # 741.0 2294.0 1698.0 5925.0 3786.0 4572.5 6508.5
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure4/HXT_protein_abundance.pdf", width =5, height =5)
barcenter = barplot(as.numeric(protein_abundance_chosen[,5]), col = apple_colors[5],
                    ylab = "Molecules number per cell")
text(barcenter, -2000, protein_abundance_chose[,2], srt = 45,xpd = TRUE)
dev.off()

