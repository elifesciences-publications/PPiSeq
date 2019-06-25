load(file = "~/Dropbox/PPiSeq_02/Working_data/GOSlim_CC.Rfile") #cc
load(file = "~/Dropbox/PPiSeq_02/Working_data/ppi_by_env_list.Rfile") 
load(file = "~/Dropbox/PPiSeq_02/Working_data/fitness_by_env_list.Rfile")
load(file = "~/Dropbox/PPiSeq_02/Working_data/bait_fitness_all_env.Rfile") 
load(file = "~/Dropbox/PPiSeq_02/Working_data/prey_fitness_all_env.Rfile")
load(file = "~/Dropbox/PPiSeq_02/Working_data/cc_bins.Rfile")
go = cc
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4",   "#E6C7C2",  "#3AA6C5", "#003366", "#AAE4E5", "#4B9C61","#DFCCB7")
env_names = names(fitness_by_env_list)

make.fit.plot1 = function(bait, env){
  m = fitness_by_env_list[[which(names(fitness_by_env_list) == env)]]
  n = ppi_by_env_list[[which(names(ppi_by_env_list) == env)]]
  p =  m[bait, 1:max(cc_bins$prey_bins)] #just first env
  plot(1:length(p), p , ylab = "Fitness", xlab = "", ylim = c(-0.4, 1), main = paste(bait, env))
  for(j in 1:(length(cc_bins$prey_bins) -1)){
    rect(cc_bins$prey_bins[j], -0.4, cc_bins$prey_bins[j+1], -0.35, col = apple_colors[j], lwd = 0, border = F)
    segments(cc_bins$prey_bins[j], median(p[(cc_bins$prey_bins[j] +1):(cc_bins$prey_bins[j+1])], na.rm = T), 
             cc_bins$prey_bins[j+1], median(p[(cc_bins$prey_bins[j] +1):(cc_bins$prey_bins[j+1])], na.rm = T), col = "red", lwd = 3)
  }
  abline(v = cc_bins$prey_bins, col = "grey")
  k = which(cc_bins$prey_cellular_component == go[rownames(m)[i]])
  rect(cc_bins$prey_bins[k], -1, cc_bins$prey_bins[k+1],2, col = rgb(0, 0, 0, alpha = 0.2))
  if(bait %in% rownames(n)){
    pos = which(colnames(m) %in% colnames(n)[which(n[bait,] == 1)])
    points(pos, p[pos], col = "blue", pch = 19 )
  }
}

par(mfrow = c(3,3))
for(i in 211:218){
  env = env_names[3]
  bait = rownames(ppi_by_env_list[[1]])[i]
  make.fit.plot1(bait, env)
}

dev.off()