## This script is to check the correlation between any pair of two environments
setwd("~/Desktop/PPiSeq_additional_data/")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
call_PPI = read.csv("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv")
vScore_PPI = read.csv("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
vScore_PPI = data.frame(vScore_PPI$PPI, vScore_PPI$Environment_number, vScore_PPI$Variation_score,
                        vScore_PPI$SD, vScore_PPI$H2O2, vScore_PPI$Forskolin, vScore_PPI$HU,
                        vScore_PPI$FK506, vScore_PPI$Raffinose, vScore_PPI$NaCl,
                        vScore_PPI$Dox, vScore_PPI$X16C)
call_PPI = data.frame(call_PPI$X, call_PPI$environment_number, 
                        call_PPI$SD, call_PPI$H2O2, call_PPI$Forskolin, call_PPI$HU,
                        call_PPI$FK506, call_PPI$Raffinose, call_PPI$NaCl,
                        call_PPI$Dox, call_PPI$X16C)

colnames(vScore_PPI) = c("PPI", "Environment_number", "vScore", "SD", "H2O2", "Hydroxyurea", "Doxorubicin",
                       "Forskolin", "Raffinose", "NaCl", "16 \u00B0C", "FK506")
colnames(call_PPI) = c("PPI", "Environment_number",  "SD",
                         "H2O2", "Forskolin", "Hydroxyurea", "FK506", "Raffinose",
                         "NaCl", "Doxorubicin", "16 \u00B0C")

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 2)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  library(scales)
  points(x,y, pch = 19, 
         col = alpha(col_group[vScore_PPI$Environment_number], 0.2))
  lines(x,x, col = "black", lty = 1)
}

png("Figures/SFigures/SFigure5/FigureS5C_pairwisefitness_scatter_plot.png", height = 15, width = 15, units = 'in', res = 600)
col_group = c("#4575b4","#74add1","#abd9e9","#e0f3f8",
              "#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
pairs(vScore_PPI[,4:12], 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      oma=c(3,3,15,3))
par(xpd = TRUE)
legend("top", title = "Number of environments in which a PPI is detected",
       legend = as.character(1:9), fill = col_group, ncol = 9, bty = "n", cex = 1.5)
dev.off()

##### Plot each pair of environments (Same figure with the above one)
PPI_count = call_PPI[,3:11]
PPI_fit = vScore_PPI[,4:12]
environment = colnames(PPI_count)
col_group = alpha(c("red", "purple", "pink", "grey"), alpha = 0.2)
for(i in 1:8){
  fit_a = PPI_fit[,i]
  count_a = PPI_count[,i]
  env_a = environment[i]
  for (j in (i+1):9){
    env_b = environment[j]
    fit_b = PPI_fit[,j]
    corr = round(cor(fit_a, fit_b), digits = 2)
    count_b = PPI_count[,j]
    group_1 = which(count_a == 1 & count_b == 1)
    group_2 = which(count_a == 1 & count_b == 0)
    group_3 = which(count_a == 0 & count_b == 1)
    group_4 = which(count_a == 0 & count_b == 0)
    plot_name = paste(env_a, env_b, ".pdf")
    pdf(plot_name, width=5, height =5)
    par(mar=c(2.5,2.5,0.5,0.5))
    plot(fit_b[group_4], fit_a[group_4], type = "p", pch = 19, col= col_group[4],
         xlim= c(0, 1.5), ylim = c(0, 1.5), xlab = NA, ylab = NA, cex.axis = 1.5)
    points(fit_a[group_3], fit_b[group_3], pch = 19, col= col_group[3])
    points(fit_a[group_2], fit_b[group_2], pch = 19, col= col_group[2])
    points(fit_a[group_1], fit_b[group_1], pch = 19, col = col_group[1])
    lines(seq(0, 1.5, by = 0.3), seq(0, 1.5, by = 0.3), col = "black")
    text(x = 0.3, y = 1.4, labels = paste("R"," = ", as.character(corr)), cex = 2.5)
    dev.off()
  }
}
