
go = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/go_slim_mapping_tab_20190405.txt",sep = "\t", header = F)) 
c = which(go[,4] == "C")
cc = go[c,5]
names(cc) = go[c,1]
save(cc, file="~/Dropbox/PPiSeq_02/Working_data/GOSlim_CC.Rfile" )
