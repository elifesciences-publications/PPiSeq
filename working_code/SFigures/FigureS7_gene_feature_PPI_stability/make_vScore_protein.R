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

setwd("~/Dropbox/PPiSeq_02/")
variation_score = csvReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
#variation_score = PPI_norm_matrix_final
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) # 2083
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
  index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                           which(PP_pair[,2] == protein_unique[i])))
  vScore_protein[i,2] = mean(1/as.numeric(variation_score[index_protein,3]))
}
save(vScore_protein, file = "Working_data_2/vScore_protein.Rfile") #for quick retrieval