
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

'''
### Figure2D A network to show the enriched biological process for PPI network
setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/")
PPI_env_summary = csvReader_T("PPI_environment_count_summary.csv")
PPI_env_count = PPI_env_summary[,1:2]
## First split these PPIs into individual proteins and check the functional enrichment of these proteins
## on a genetic interaction network. Split PPI into different groups: all_env, one_env, two_env, three_env, four_more_env
PPI_unique_protein = function(PPI_env_count){
all_PPI= split_string_vector(PPI_env_count[,1])
all_PPI_unique = unique(c(all_PPI[,1], all_PPI[,2]))
return(all_PPI_unique)
}
all_env = PPI_unique_protein(PPI_env_count) # 2111
csvWriter(all_env, "unique_protein_from_PPI/all_environment_unique_PPI.csv")
#One_env
PPI_env_one = PPI_env_count[which(as.numeric(PPI_env_count[,2]) == 1),]
one_env = PPI_unique_protein(PPI_env_one) # 2089
csvWriter(one_env, "unique_protein_from_PPI/One_environment_unique_PPI.csv")
#Two_env
PPI_env_two = PPI_env_count[which(as.numeric(PPI_env_count[,2]) == 2),]
two_env = PPI_unique_protein(PPI_env_two) # 754
csvWriter(two_env, "unique_protein_from_PPI/Two_environment_unique_PPI.csv")
#Three_env
PPI_env_three = PPI_env_count[which(as.numeric(PPI_env_count[,2]) == 3),]
three_env = PPI_unique_protein(PPI_env_three) # 462
csvWriter(three_env, "unique_protein_from_PPI/Three_environment_unique_PPI.csv")
#Four_more_env
PPI_env_four = PPI_env_count[which(as.numeric(PPI_env_count[,2]) >= 4),]
four_env = PPI_unique_protein(PPI_env_four) # 800
csvWriter(four_env, "unique_protein_from_PPI/Four_more_environment_unique_PPI.csv")

# First, I will build a network by using igraph (just show the reported and unreported PPIs no annatation)
library("igraph")
reported_PPI = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/multiple_validated_PPI.csv")
name_exchange = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Systematic_standard_protein.csv")
PPI_env_reported = match_both_direction(PPI_env_count, reported_PPI[,1]) # 623
PPI_env_split = split_string_vector(PPI_env_count[,1]) 
a = name_exchange[match(PPI_env_split[,1], name_exchange[,1]),2]
b = name_exchange[match(PPI_env_split[,2], name_exchange[,1]),2]
type = rep(1, nrow(PPI_env_count))
for(i in 1:nrow(PPI_env_count)){
if(PPI_env_count[i,1] %in% PPI_env_reported[,1]){
type[i] = 2
}
}
PPI_net = data.frame(a, b, type)
#### Plot network
net <- graph_from_data_frame(d= PPI_net, directed = F)
net <- simplify(net, remove.multiple = T, remove.loops = T)
deg <- igraph::degree(net, normalized = T)
#deg_count <- log2(igraph::degree(net))
#V(net)$size <- deg*50 # Node size stands for degree of protein
#E(net)$width <- E(net)$weight * 4 # Edge thickness represents the dynamics
edge_color = apple_colors[c(6, 7)]
E(net)$color = edge_color[E(net)$type]
l <- layout_with_dh(net)
#plot(net, layout = l)
pdf("~/Desktop/PPI_network.pdf", height = 5, width = 5)
par(mar = c(2,1,0,0))
plot(net,  vertex.frame.color = apple_colors[3], vertex.label.cex = 0.25,vertex.color = apple_colors[3], layout = l)
legend(x =0.25, y = -1, c("Previously unreported", "Previously reported"), 
lty = c(1,1), col = apple_colors[c(6,7)], bty= "n", lwd = 2, cex = 0.8)
dev.off()
'''
#######################################
## Figure 2D make a barplot to show the ratio of coannations for PPIs that are 
## detected in different number of environments
## Use python code to make a dictionary of go term for each gene 
## and then check the overlap co-annoations between two interacting partners
## I will generate matrix data in python, and then make plot in R
