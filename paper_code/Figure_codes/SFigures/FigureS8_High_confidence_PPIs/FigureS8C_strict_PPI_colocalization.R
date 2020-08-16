# This script compares PPIs discoverd by PPiSeq to Localization data from Yolanda Chong et al. and Gene Ontology
# It finds that PPIs that are found in  more environments are more likely to colocalize 

setwd("~/Desktop/PPiSeq_additional_data/")
load("Outsourced_datasets/Colocalization/yolonda_benign_loc_matrix.Rfile") # Load GFP colocalizaiton data
load("Outsourced_datasets/Colocalization/GO_loc_matrix.Rfile") # Load GO colocalization data
GO_loc_matrix = GO_loc_matrix[,c(1, 5:24)] #remove general terms  "other", "membrane",  "cellular_component" 

#Import PPIs by environemnt
PPI_count_filter = read.csv("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv")
PPI_variation = read.csv("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter_strict_threshold.csv")
ppi_by_env = PPI_count_filter[,3:ncol(PPI_count_filter)]
rownames(ppi_by_env) <- PPI_count_filter[,1]

a = sapply(as.character(rownames(ppi_by_env)), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
ppi = x #detected PPIs across all environments


#Find Colacalization b/w PPI functions
find_colacalized = function(ppi, yolonda_benign_loc_matrix){
  wt = yolonda_benign_loc_matrix
  x = 1:nrow(ppi)
  x[] = NA
  for(i in 1:length(x)){
    a = ppi[i,1]
    b = ppi[i,2]
    if(a %in% rownames(wt) & b %in% rownames(wt)){
      t = wt[c(a,b),]
      t = t > 0
      u = apply(t, 2, sum)
      if(sum(u == 2) > 0){
        x[i] = 1
      }else{
        x[i] = 0
      }
    }
  }
  return(x)
}

percent_colacalized = function(x){
  #print(sum(x == 1, na.rm = T))
  #print(sum(x == 0, na.rm = T))
  sum(x == 1, na.rm = T)/ (sum(x == 1, na.rm = T) + sum(x == 0, na.rm = T))
}

env.number = apply(ppi_by_env, 1, sum, na.rm = T)
number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(number.of.env.colcalized.in.DMSO) = 1:9
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == i ,] #PPIs seen in i env
  y = find_colacalized(x, GO_loc_matrix)
  number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
}
number.of.env.colcalized.bootstrap = matrix(NA, 1000, ncol(ppi_by_env))
colnames(number.of.env.colcalized.bootstrap) = 1:9
for(j in 1:1000){
  b = sample(1:nrow(ppi_by_env), nrow(ppi_by_env), replace = T) #the new bootstrapped indexes
  ppi_by_env_b = ppi_by_env[b,] #the new bootstrapped ppis
  ppi_b = ppi[b,]
  env.number_b = apply(ppi_by_env_b, 1, sum, na.rm = T)
  for(i in 1:ncol(ppi_by_env_b)){
    x = ppi_b[env.number_b == i ,] #PPIs seen in i env
    y = find_colacalized(x, yolonda_benign_loc_matrix)
    number.of.env.colcalized.bootstrap[j,i] = percent_colacalized(y)
  }
}
m = apply(number.of.env.colcalized.bootstrap, 2, mean, na.rm = T)*100
sd = apply(number.of.env.colcalized.bootstrap, 2, sd, na.rm = T)*100
sem = sd/sqrt(1000)
u = m+sd
l = m-sd

#Do a bootstrap with GO to get error bars -- sample PPIs with replacement
number.of.env.colcalized.bootstrap = matrix(NA, 1000, ncol(ppi_by_env))
colnames(number.of.env.colcalized.bootstrap) = 1:9
for(j in 1:1000){
  b = sample(1:nrow(ppi_by_env), nrow(ppi_by_env), replace = T) #the new bootstrapped indexes
  ppi_by_env_b = ppi_by_env[b,] #the new bootstrapped ppis
  ppi_b = ppi[b,]
  env.number_b = apply(ppi_by_env_b, 1, sum, na.rm = T)
  for(i in 1:ncol(ppi_by_env_b)){
    x = ppi_b[env.number_b == i ,] #PPIs seen in i env
    y = find_colacalized(x, GO_loc_matrix)
    number.of.env.colcalized.bootstrap[j,i] = percent_colacalized(y)
  }
}
mgo = apply(number.of.env.colcalized.bootstrap, 2, mean, na.rm = T)*100
sdgo = apply(number.of.env.colcalized.bootstrap, 2, sd, na.rm = T)*100
semgo = sd/sqrt(1000)
ugo = mgo+sdgo
lgo = mgo-sdgo

pdf(file = "Figures/SFigures/SFigure8/FigureS8C_colocalization_rate.pdf", height = 7, width = 4)
plot(1:9, m, xlab = "Environments in which a PPI is observed", 
     ylab = "Percent colocalized", type = 'l', ylim = c(30,100), pch = 16, cex = 1.7, lwd =2, col = "black" ) 
points(1:9, m,  ylim = c(30,90), pch = 16, cex = 1.7, col = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027"))
points(1:9, mgo,type = 'l', lwd = 2, lty = "dashed")
points(1:9, mgo,   pch = 16, cex = 1.7, 
       col = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027"))
arrows( 1:9, u, 1:9,l,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
arrows( 1:9, ugo, 1:9,lgo,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
legend(3, 40, legend = c("Gene Ontology", "Fluoresence"),lty=c("dashed", "solid"), bty= "n" )
dev.off()



