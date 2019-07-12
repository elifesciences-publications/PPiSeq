load(file = "~/Dropbox/PPiSeq_02/Working_data/GOSlim_CC.Rfile") #cc
load(file = "~/Dropbox/PPiSeq_02/Working_data/ppi_by_env_list.Rfile") 
load(file = "~/Dropbox/PPiSeq_02/Working_data/fitness_by_env_list.Rfile")
load(file = "~/Dropbox/PPiSeq_02/Working_data/bait_fitness_all_env.Rfile") 
load(file = "~/Dropbox/PPiSeq_02/Working_data/prey_fitness_all_env.Rfile")
load(file = "~/Dropbox/PPiSeq_02/Working_data/cc_bins.Rfile")
go = cc
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4",   "#E6C7C2",  "#3AA6C5", "#003366", "#AAE4E5", "#4B9C61","#DFCCB7")
compartment.colors = c(apple_colors, apple_colors, apple_colors, apple_colors, apple_colors, apple_colors, apple_colors,
                       apple_colors, apple_colors)
env.colors = rainbow(9, s = 0.5)
names(env.colors) = names(fitness_by_env_list)
env_names = names(fitness_by_env_list)
bcomp.names = cc_bins[[1]]




#Compare background fitness of bait with or without a PPI in that compartment
neg_fitness_by_comp = ppi_by_env_list #this will be the fitness values in each compartment for all bait that do not have a PPI in that env
for(j in 1:length(ppi_by_env_list)){
  neg_fitness_by_comp_env = list()
  env =  env_names[j]
  m = fitness_by_env_list[[which(names(fitness_by_env_list) == env)]]
  n = ppi_by_env_list[[which(names(ppi_by_env_list) == env)]]
  pdf(file = paste("~/Dropbox/PPiSeq_02/Working_figure/Sasha/density_plot_compartment_", env, ".pdf", sep = ""), height = 14, width = 14)
  par(mfrow = c(4,4))
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE,ann=FALSE)
  legend(0, 10, c("No PPIs within this compartment", "At least one PPI within this compartment"), fill = c( rgb(1,0,0,0.2), rgb(0,0,1,0.2) ))
  for(i in 1:length(bcomp.names)){
    comp = bcomp.names[i]
    if(comp %in% cc_bins$prey_cellular_component){
      prey.bins = cc_bins$prey_bins[which(cc_bins$prey_cellular_component == comp):(which(cc_bins$prey_cellular_component == comp)+1) ]
      mc = m[,(prey.bins[1]+1):prey.bins[2]]
      nc = n[,(prey.bins[1]+1):prey.bins[2]]
      pos.int = apply(nc, 1, sum, na.rm = T) # Number of interactions in that environment
      mc[which(nc == 1)] = NA #remove all discovered interactions from the calculation
      p = mc[which(pos.int > 0), ]
      ne = mc[which(pos.int == 0), ]
      pos = density(p[!is.na(p)])
      neg = density(ne[!is.na(ne)])
      neg_fitness_by_comp_env[[i]] = ne[!is.na(ne)]
      plot(neg, main = comp, xlab = "Fitness", xlim = c(-0.4, 0.4), ylim = c(0,10))
      polygon(neg, col = rgb(1,0,0,0.2))
      points(pos, type = 'l')
      polygon(pos, col = rgb(0,0,1,0.2))
      legend(-0.4, 6, c(round(mean(ne, na.rm = T), digits = 3),round( mean(p, na.rm = T), digits = 3)), 
             fill = c( rgb(1,0,0,0.2), rgb(0,0,1,0.2)), title = "mean")
    }
  }
  names(neg_fitness_by_comp_env) = cc_bins$bait_cellular_component[1:13]
  neg_fitness_by_comp[[j]] = neg_fitness_by_comp_env
  dev.off()
}
save(neg_fitness_by_comp, file = "~/Dropbox/PPiSeq_02/Working_data/neg_fitness_by_comp.Rfile")


#Compare background fitness of bait with or without a PPI in any compartment
pdf(file = "~/Dropbox/PPiSeq_02/Working_figure/Sasha/density_plot_compartment_control.pdf", height = 14, width = 14)
par(mfrow = c(4,4))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE,ann=FALSE)
legend(0, 10, c("No PPIs within any compartment", "At least one PPI within any compartment"), fill = c( rgb(1,0,0,0.2), rgb(0,0,1,0.2) ))
for(i in 1:length(bcomp.names)){
  comp = bcomp.names[i]
  env = env_names[1]
  prey.bins = cc_bins$prey_bins[which(cc_bins$prey_cellular_component == comp):(which(cc_bins$prey_cellular_component == comp)+1) ]
  mc = m[,(prey.bins[1]+1):prey.bins[2]]
  nc = n[,(prey.bins[1]+1):prey.bins[2]]
  pos.int = apply(n, 1, sum, na.rm = T) # Number of interactions in any environment
  mc[which(nc == 1)] = NA #remove all discovered interactions from the calculation
  p = mc[which(pos.int > 0), ]
  ne = mc[which(pos.int == 0), ]
  pos = density(p[!is.na(p)])
  neg = density(ne[!is.na(ne)])
  plot(neg, main = comp, xlab = "Fitness", xlim = c(-0.4, 0.4))
  polygon(neg, col = rgb(1,0,0,0.2))
  points(pos, type = 'l')
  polygon(pos, col = rgb(0,0,1,0.2))
  legend(-0.4, 6, c(round(mean(ne, na.rm = T), digits = 3),round( mean(p, na.rm = T), digits = 3)), fill = c( rgb(1,0,0,0.2), rgb(0,0,1,0.2)), title = "mean")
}
dev.off()


#Get a p-value for each bait in each compartment vs. negative
comp_pvalue_list = list()
for(k in 1:length(env_names)){
  env =  env_names[k]
  m = fitness_by_env_list[[which(names(fitness_by_env_list) == env)]]
  x = matrix(NA, nrow(m),  sum(cc_bins$bait_cellular_component %in% cc_bins$prey_cellular_component))
  colnames(x) = cc_bins$bait_cellular_component[cc_bins$bait_cellular_component %in% cc_bins$prey_cellular_component]
  rownames(x) = rownames(ppi_by_env_list[[1]])
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      comp = colnames(x)[j]
      bait = rownames(x)[i]
      prey.bins = cc_bins$prey_bins[which(cc_bins$prey_cellular_component == comp):(which(cc_bins$prey_cellular_component == comp)+1) ]
      test = m[bait,(prey.bins[1]+1):prey.bins[2]]
      neg = neg_fitness_by_comp[[env]][[comp]]
      if(sum(!is.na(test))>1){
        x[i,j] = t.test(test, neg, alternative = "greater", na.action = "omit")[[3]]
      }
    }
  }
  comp_pvalue_list[[k]] = x
}
names(comp_pvalue_list) = env_names
save(comp_pvalue_list, file = "~/Dropbox/PPiSeq_02/Working_data/comp_pvalue_list.Rfile")

#Compare our predictions to GoSlim in DMSO
x = comp_pvalue_list$DMSO
y = x[,1:3]
colnames(y) = c("GO", "PPiSeq Primary", "PPiSeq others")
y[] = NA
for(i in rownames(y)){
  y[i,1] = go[i]
  if(length(colnames(x)[x[i,] < 0.05] > 0)){
    y[i,2] = paste(colnames(x)[x[i,] < 0.05], collapse = ", ")
  }
}

bait = rownames(x)[1]
go[bait]
x[bait,]


make.fit.plot1 = function(bait, env){
  m = fitness_by_env_list[[which(names(fitness_by_env_list) == env)]]
  n = ppi_by_env_list[[which(names(ppi_by_env_list) == env)]]
  p =  m[bait, 1:max(cc_bins$prey_bins)] #just first env
  plot(1:length(p), p , ylab = "Fitness", xlab = "", ylim = c(-0.4, 0.6), main = paste(bait, env))
  for(j in 1:(length(cc_bins$prey_bins) -1)){
    rect(cc_bins$prey_bins[j], -0.4, cc_bins$prey_bins[j+1], -0.35, col = apple_colors[j], lwd = 0, border = F)
    segments(cc_bins$prey_bins[j], median(p[(cc_bins$prey_bins[j] +1):(cc_bins$prey_bins[j+1])], na.rm = T), 
             cc_bins$prey_bins[j+1], median(p[(cc_bins$prey_bins[j] +1):(cc_bins$prey_bins[j+1])], na.rm = T), col = "red", lwd = 3)
  }
  abline(v = cc_bins$prey_bins, col = "grey")
  k = which(cc_bins$prey_cellular_component == go[bait])
  rect(cc_bins$prey_bins[k], -1, cc_bins$prey_bins[k+1],2, col = rgb(0, 0, 0, alpha = 0.2))
  if(bait %in% rownames(n)){
    pos = which(colnames(m) %in% colnames(n)[which(n[bait,] == 1)])
    points(pos, p[pos], col = "blue", pch = 19 )
  }
}

make.fit.plot2 = function(bait){
  m = fitness_by_env_list
  n = ppi_by_env_list[names(m)]
  p =  m[[1]][bait,] #just first env
  q = m[[1]][bait,]
  q[] = NA
  q[colnames(n[[1]])] = n[[1]][bait,]
  r = q
  compartment.bins = cc_bins$prey_bins
  for(i in 2:length(m)){
    p = c(p,m[[i]][bait,])
    q = m[[i]][bait,]
    q[] = NA
    q[colnames(n[[i]])] = n[[i]][bait,]
    r = c(r,q)
    c = cc_bins$prey_bins + (i-1)*1732
    compartment.bins = c(compartment.bins, c)
  }
  plot(1:length(p), p , ylab = "Fitness", xlab = "", ylim = c(-0.4, 1), main = paste(bait, env), type = 'l')
  env.bins = c(0, (1:9)*ncol(m[[1]]))
  for(j in 1:(length(compartment.bins) -1)){
    rect(compartment.bins[j], -0.4, compartment.bins[j+1], -0.35, col = compartment.colors[j], lwd = 0, border = F)
  }
  for(j in 1:(length(env.bins) -1)){
    rect(env.bins[j], -0.34, env.bins[j+1], -0.3, col = env.colors[j], lwd = 0, border = F)
  }
  abline(v = env.bins, col = "grey")
  if(bait %in% rownames(n[[1]])){
    pos = which(r == 1)
    points(pos, p[pos], col = "blue", pch = 19 )
  }
  p.no.na = p[!is.na(p)]
  l = lowess(which(!is.na(p)), p.no.na, f = .000005)
  points(l$x, l$y, col = "red", type = "l", lwd = 3)
}

make.fit.plot3 = function(bait){
  m = fitness_by_env_list[[1]]
  plot(1:ncol(m), m[1,] , ylab = "Fitness", xlab = "", ylim = c(-0.4, 1), main = paste(bait), col = "white")
  for(env in env_names){
    m = fitness_by_env_list[[which(names(fitness_by_env_list) == env)]]
    n = ppi_by_env_list[[which(names(ppi_by_env_list) == env)]]
    p =  m[bait, 1:max(cc_bins$prey_bins)] #the data
    #plot(1:length(p), p , ylab = "Fitness", xlab = "", ylim = c(-0.4, 0.6), main = paste(bait, env))
    for(j in 1:(length(cc_bins$prey_bins) -1)){
      rect(cc_bins$prey_bins[j], -0.4, cc_bins$prey_bins[j+1], -0.35, col = apple_colors[j], lwd = 0, border = F)
      segments(cc_bins$prey_bins[j], median(p[(cc_bins$prey_bins[j] +1):(cc_bins$prey_bins[j+1])], na.rm = T), 
               cc_bins$prey_bins[j+1], median(p[(cc_bins$prey_bins[j] +1):(cc_bins$prey_bins[j+1])], na.rm = T), col = env.colors[env], lwd = 2)
    }
    abline(v = cc_bins$prey_bins, col = "grey")
    k = which(cc_bins$prey_cellular_component == go[bait])
    #rect(cc_bins$prey_bins[k], -1, cc_bins$prey_bins[k+1],2, col = rgb(0, 0, 0, alpha = 0.2))
    #p.no.na = p[!is.na(p)]
    #l = lowess(which(!is.na(p)), p.no.na, f = .05)
    #points(l$x, l$y, col = "red", type = "l", lwd = 3)
    if(bait %in% rownames(n)){
      pos = which(colnames(m) %in% colnames(n)[which(n[bait,] == 1)])
      points(pos, p[pos], col = env.colors[env], pch = 19 )
    }
  }
}

plot.prey.legend = function(){
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE,ann=FALSE)
  legend(0,10, cc_bins$prey_cellular_component[1:7], fill =  apple_colors[1:7], bty = "n")
  legend(5,10, cc_bins$prey_cellular_component[8:14], fill =  apple_colors[8:14], bty = "n")
}

plot.env.legend = function(){
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE,ann=FALSE)
  legend(0,10, env_names[1:5], fill =  env.colors[1:5], bty = "n")
  legend(5,10, env_names[6:9], fill =  env.colors[6:9], bty = "n")
}

par(mfrow = c(3,3))
for(i in sample(nrow(ppi_by_env_list[[1]]), 9)){
  env = env_names[1]
  bait = rownames(ppi_by_env_list[[1]])[i]
  make.fit.plot1(bait, env)
}

par(mfrow = c(3,3))
for(i in sample(nrow(ppi_by_env_list[[1]]), 9)){
  bait = rownames(ppi_by_env_list[[1]])[i]
  make.fit.plot2(bait)
}

par(mfrow = c(3,3))
plot.env.legend()
for(i in sample(nrow(ppi_by_env_list[[1]]), 8)){
  bait = rownames(ppi_by_env_list[[1]])[i]
  make.fit.plot3(bait)
}

pdf(file = "~/Desktop/fitplot.pdf", width = 16, height = 10)
par(mfrow = c(3,3))
plot.prey.legend()
b = c("YAL053W", "YKR048C", "YJR042W", "YMR109W", "YIL027C", "YHR094C", "YMR010W", "YBL007C")
for(i in b){
  env = env_names[1]
  bait = i
  make.fit.plot1(bait, env)
}
dev.off()

b = c("YAL053W", "YKR048C", "YJR042W", "YMR109W", "YIL027C", "YHR094C", "YMR010W", "YBL007C")
for(i in b){
  bait = i
  make.fit.plot3(bait)
}

f = fitness_by_env_list
p = ppi_by_env_list
c = cc_cins


