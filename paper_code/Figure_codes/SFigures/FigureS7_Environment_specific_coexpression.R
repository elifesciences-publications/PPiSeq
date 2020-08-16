#A plot of co-expression mutual rank by environment number
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
load("~/Desktop/PPiSeq_additional_data/Outsourced_datasets/Coexpression/CoExpressDB.Rfile") # Load coexpression database
ppi_coex = coex


# SD
coexpression_env = function(column_num){
  PPI_count = csvReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv")
  v = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
  v = v[which(PPI_count[,column_num] == "1"),]
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
  return(df)
}
library(ggplot2)
## SD
column_num = 3
pdf_name = 'Figures/SFigures/SFigure7/SD_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=45), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

##H2O2
column_num = 4
pdf_name = 'Figures/SFigures/SFigure7/H2O2_coexpression_by_env_number.pdf'
df = coexpression_env(column_num, pdf_name)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

##HU
column_num = 5
pdf_name = 'Figures/SFigures/SFigure7/HU_coexpression_by_env_number.pdf'
df = coexpression_env(column_num, pdf_name)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

##Dox
column_num = 6
pdf_name = 'Figures/SFigures/SFigure7/Dox_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()


## Forskolin
column_num = 7
pdf_name = 'Figures/SFigures/SFigure7/Forskolin_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

## Raffinose
column_num = 8
pdf_name = 'Figures/SFigures/SFigure7/Raffinose_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()


## NaCl
column_num = 9
pdf_name = 'Figures/SFigures/SFigure7/NaCl_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

## 16C
column_num = 10
pdf_name = 'Figures/SFigures/SFigure7/16C_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()

## FK506
column_num = 11
pdf_name = 'Figures/SFigures/SFigure7/FK506_coexpression_by_env_number.pdf'
df = coexpression_env(column_num)
pdf(file = pdf_name, height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=90), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()
