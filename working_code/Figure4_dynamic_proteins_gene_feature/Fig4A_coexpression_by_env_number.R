
#### Caculate an variation score for each PPI based on its fitness values across
#### different environments. Use the normalized value and Keep the primary normalized fitess even a PPI is called negative
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

setwd("~/Dropbox/PPiseq_02/")

#A plot of co-expression mutual rank by environment number
require(ggplot2)
#load("~/Dropbox/PPiSeq_02/Working_data_2/variation_score.Rfile") #variation_score
load("~/Dropbox/PPiSeq_02/Working_data_2/Coexpression/CoExpressDB.Rfile") #coex

#load data
ppi_coex = coex
v = csvReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
a = sapply(as.character(v[,1]), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
pairs = x
dynamicity = as.numeric(v[,3])
env.number = as.numeric(v[,2])


#Get CoExpression for all PPIs
x = 1:length(dynamicity)
x[] = NA
for(i in 1:length(x)){
  if(pairs[i,1] %in% colnames(coex) & pairs[i,2] %in% colnames(coex) )x[i] = coex[pairs[i,1], pairs[i,2]]
}
ppi_coex = x

#make plot of co-expression by environment number
df = data.frame(ppi_coex, env.number)
df$env.number = as.factor(df$env.number)

pdf(file = "~/Dropbox/PPiSeq_02/Working_figure/Figure4_date_party_hubs/Figure4A_coexpression_by_env_number.pdf", height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=45), 
         axis.text.x = element_text(size=10), 
         panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

