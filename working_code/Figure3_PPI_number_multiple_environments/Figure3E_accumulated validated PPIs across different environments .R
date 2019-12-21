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

#(2) Get all the possible orders of 9 environments, and check the PPI number changes
setwd("~/Dropbox/PPiSeq_02/")
count_summary = csvReader_T("Paper_data/Useful_datasets/Validated_PPI_environment_count_summary_SD_merge_filter.csv") # 12981
### Followsing this order DMSO, H2O2, HU, Dox, Forskolin, Raffinose, NaCl, 16C, FK506
row_num = c(factorial(9), factorial(8), factorial(7), factorial(6), factorial(5), factorial(4), factorial(3), factorial(2),factorial(1))

matrix_count_1 = matrix(0, factorial(9), 9)


matrix_constructer = rep(NA, nrow(count_summary))
all_env = 3:11
# First column
bulk_1 = row_num[2]
index_1 = bulk_1 -1 
for(i in 1:9){
        column_chosen = as.numeric(count_summary[,i +2])
        matrix_chosen = data.frame(column_chosen, matrix_constructer)
        mean_row = rowMeans(matrix_chosen, na.rm = T)
        environ_1_1 = sum(mean_row, na.rm = T)
        index= (bulk_1 * i - index_1):(bulk_1*i)
        matrix_count_1[index,1] = environ_1_1
}

# Second column
bulk_2 = row_num[3]
index_2 = bulk_2-1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h + 2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                column_chosen_2 = as.numeric(count_summary[,env_remaining_1[i]])
                matrix_chosen = data.frame(column_chosen_1, column_chosen_2)
                mean_row = rowMeans(matrix_chosen, na.rm = T)
                environ_2_1 = sum(mean_row, na.rm = T)
                index = (bulk_2 * i - index_2):(bulk_2*i) + env1_index
                matrix_count_1[index,2] = environ_2_1
        }
}


# Third column
bulk_3 = row_num[4] #720
index_3 = bulk_3 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                        column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                        column_chosen_3 = as.numeric(count_summary[,env_remaining_2[j]])
                        matrix_chosen = data.frame(column_chosen_1, column_chosen_2, column_chosen_3)
                        mean_row = rowMeans(matrix_chosen, na.rm = T)
                        environ_3_1 = sum(mean_row, na.rm = T)
                        index = (bulk_3 * j - index_3):(bulk_3*j) + env2_index
                        matrix_count_1[index, 3] = environ_3_1 
                }
                
        } 
}


# Fourth column
bulk_4 = row_num[5] # 120
index_4 = bulk_4 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                column_chosen_4 = as.numeric(count_summary[,env_remaining_3[k]])
                                matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                           column_chosen_3, column_chosen_4)
                                mean_row = rowMeans(matrix_chosen, na.rm = T)
                                environ_4_1 = sum(mean_row, na.rm = T)
                                index = (bulk_4 * k - index_4):(bulk_4*k) + env3_index
                                matrix_count_1[index, 4] = environ_4_1 
                        }
                          
                }
                
        } 
}


matrix_count_1[362880 -120,]
matrix_count_1[362880-120 + 1,]
length(unique(matrix_count_1[,4])) # 121
sum(matrix_count_1[,4] > matrix_count_1[,3])

# Fifth column
bulk_5 = row_num[6] # 24
index_5 = bulk_5 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                        column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                        column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                        column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                        column_chosen_5 = as.numeric(count_summary[,env_remaining_4[l]])
                                        matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                   column_chosen_3, column_chosen_4,
                                                                   column_chosen_5)
                                        mean_row = rowMeans(matrix_chosen, na.rm = T)
                                        environ_5_1 = sum(mean_row, na.rm = T)
                                        index = (bulk_5 * l - index_5):(bulk_5*l) + env4_index
                                        matrix_count_1[index, 5] = environ_5_1
                                       
                                }
                                    
                        }
                        
                }
                
        } 
}
matrix_count_1[362880-24,]
matrix_count_1[362880 - 24 +1,]
matrix_count_1[362880,]
sum(matrix_count_1[,5] > matrix_count_1[,4])


# Sixth column
bulk_6 = row_num[7] # 6
index_6 = bulk_6 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        env5_index = (l -1) * bulk_5 + env4_index
                                        chosen_env_5 = env_remaining_4[l]
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        for(m in 1:4){
                                                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                                column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                                column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                                column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                                column_chosen_5 = as.numeric(count_summary[,chosen_env_5])
                                                column_chosen_6 = as.numeric(count_summary[,env_remaining_5[m]])
                                                matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                           column_chosen_3, column_chosen_4,
                                                                           column_chosen_5, column_chosen_6)
                                                mean_row = rowMeans(matrix_chosen, na.rm = T)
                                                environ_6_1 = sum(mean_row, na.rm = T)
                                                
                                                index = (bulk_6 * m - index_6):(bulk_6*m) + env5_index
                                                matrix_count_1[index, 6] = environ_6_1 
                                                
                                        }
                                            
                                }
                                
                        }
                        
                }
                
        } 
}
matrix_count_1[1:7,]
sum(matrix_count_1[,6] > matrix_count_1[,5]) # 362880


# Seventh column
bulk_7 = row_num[8] # 2
index_7 = bulk_7 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        env5_index = (l -1) * bulk_5 + env4_index
                                        chosen_env_5 = env_remaining_4[l]
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        for(m in 1:4){
                                                env6_index = (m -1) * bulk_6 + env5_index
                                                chosen_env_6 = env_remaining_5[m]
                                                env_remaining_6 = env_remaining_5[which(env_remaining_5 != chosen_env_6)]
                                                for(n in 1:3){
                                                        column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                                        column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                                        column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                                        column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                                        column_chosen_5 = as.numeric(count_summary[,chosen_env_5])
                                                        column_chosen_6 = as.numeric(count_summary[,chosen_env_6])
                                                        column_chosen_7 = as.numeric(count_summary[,env_remaining_6[n]])
                                                        matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                                   column_chosen_3, column_chosen_4,
                                                                                   column_chosen_5, column_chosen_6,
                                                                                   column_chosen_7)
                                                        mean_row = rowMeans(matrix_chosen, na.rm = T)
                                                        environ_7_1 = sum(mean_row, na.rm = T)
                                                        index = (bulk_7 * n - index_7):(bulk_7*n) + env6_index
                                                        matrix_count_1[index, 7] = environ_7_1 
                                                       
                                                }
                                                 
                                        }
                                        
                                }
                                
                        }
                        
                }
                
        } 
}
matrix_count_1[1:7,]
matrix_count_1[40300:40320,]
sum(matrix_count_1[,7] > matrix_count_1[,6])

# Eighth column
bulk_8 = row_num[9] # 1
index_8 = bulk_8 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        env5_index = (l -1) * bulk_5 + env4_index
                                        chosen_env_5 = env_remaining_4[l]
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        for(m in 1:4){
                                                env6_index = (m -1) * bulk_6 + env5_index
                                                chosen_env_6 = env_remaining_5[m]
                                                env_remaining_6 = env_remaining_5[which(env_remaining_5 != chosen_env_6)]
                                                for(n in 1:3){
                                                        env7_index = (n -1) * bulk_7 + env6_index
                                                        chosen_env_7 = env_remaining_6[n]
                                                        env_remaining_7 = env_remaining_6[which(env_remaining_6 != chosen_env_7)]
                                                        for(o in 1:2){
                                                                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                                                column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                                                column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                                                column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                                                column_chosen_5 = as.numeric(count_summary[,chosen_env_5])
                                                                column_chosen_6 = as.numeric(count_summary[,chosen_env_6])
                                                                column_chosen_7 = as.numeric(count_summary[,chosen_env_7])
                                                                column_chosen_8 = as.numeric(count_summary[,env_remaining_7[o]])
                                                                matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                                           column_chosen_3, column_chosen_4,
                                                                                           column_chosen_5, column_chosen_6,
                                                                                           column_chosen_7, column_chosen_8)
                                                                mean_row = rowMeans(matrix_chosen, na.rm = T)
                                                                environ_8_1 = sum(mean_row, na.rm = T)
                                                                
                                                                index = (bulk_8 * o - index_8):(bulk_8*o) + env7_index
                                                                matrix_count_1[index, 8] = environ_8_1
                                                                
                                                        }
                                                        
                                                }
                                                
                                        }
                                        
                                }
                                
                        }
                        
                }
                
        } 
}

matrix_count_1[362870:362880,]
sum(matrix_count_1[,8] > matrix_count_1[,7])
length(unique(matrix_count_1[,8])) #9


# Ninth column, all have the same value by taking the mean value for 
matrix_chosen = data.frame(as.numeric(count_summary[,3]), as.numeric(count_summary[,4]), as.numeric(count_summary[,5]),
                           as.numeric(count_summary[,6]), as.numeric(count_summary[,7]), as.numeric(count_summary[,8]),
                           as.numeric(count_summary[,9]), as.numeric(count_summary[,10]), as.numeric(count_summary[,11]))

mean_all = rowMeans(matrix_chosen, na.rm = T)
matrix_count_1[,9] = sum(mean_all, na.rm = T)

colnames(matrix_count_1) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")

csvWriter(matrix_count_1, "Working_data_2/Accumulated_PPI_all/PPI_number_environment_accumulation_all_order_1.csv")

##### Make a plot to show these accumulated counts
setwd("~/Dropbox/PPiSeq_02/Working_data_2/Accumulated_PPI_all/")
matrix_count_1 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_1.csv")
col_mean_1 = colMeans(matrix_count_1)

library(scales)
col_group = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
col_chosen = alpha(col_group, 0.05)
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure3_accessory_PPIs/Figure3D_PPI_number_accumulation_environments_validated.pdf", width= 5.5, height =5.5)
#pdf("~/Desktop/Figure2B_PPI_number_accumulation_environments_combine.pdf", width =5, height =5)
par(mar = c(4,4,2,4))
plot(1:9, as.numeric(matrix_count_1[1,]), xlim = c(1,9), ylim = c(0,12000), type = "l",
     col = col_chosen[1],lwd = 0.5, axes=F,
     xlab = "Number of environments assayed",
     ylab = "Number of PPIs", main = NA, bty = "n")
axis(2, at= seq(0, 12000, by = 2000), labels = seq(0, 12000, by = 2000))
axis(1, at= 1:9, labels = 1:9)
for(i in 2:nrow(matrix_count_1)){
        lines(1:9, as.numeric(matrix_count_1[i,]), col = col_chosen[1], lwd = 0.5)
}

lines(1:9, col_mean_1, col = apple_colors[11], lwd = 2, lty = 2)

#legend(1,18000, as.character(1:8),lty = 1, col = col_group[1:8], ncol = 3, bty= "n")
#text(4, 18500, "Minimum number of environments\nin which the PPI is observed", xpd = TRUE)
dev.off()



