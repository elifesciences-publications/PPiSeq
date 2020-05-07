## title: preprocessing localization data
## author: darach

library(tidyverse)
library(DBI)

db_ld <- dbConnect(RSQLite::SQLite(),"localization_data.sqlite")

# reading in the chong2015 files exported from the excel file, manually
localization <- list.files("./","wt.\\.csv") %>% tibble(filenames=.) %>% 
    mutate(rawfile=map(filenames,read_csv)) %>% unnest(rawfile) %>% 
    mutate(`Gene name`=ifelse(`Gene name`=="1-Oct","OCT1",`Gene name`)) %>% 
    arrange(`Gene name`) %>% group_by(ORF,`Gene name`) %>% 
    summarize_at(vars(Bud:`Vacuole/Vacuolar Membrane`),sum) %>% 
    pivot_longer(cols=-c(ORF,`Gene name`)) %>%
    group_by(ORF,`Gene name`) %>% 
    mutate(value=value/sum(value,na.rm=T))

localization_matrix <- localization %>% 
    ungroup() %>% select(ORF,name,value) %>% 
    mutate(value=ifelse(value<0.1,0,value)) %>%
    pivot_wider(names_from="name",values_from="value") %>% 
    {tmp <- as.matrix(.[,-1]);rownames(tmp)<-.$ORF;tmp} 

jacardz <- sapply(1:nrow(localization_matrix),
        function(i){sapply(1:nrow(localization_matrix),
            function(j){
                if (j >= i){
                    return(NA)
                } else {
                    index <- localization_matrix[i,]|localization_matrix[j,]
                    return(
                        sum(
                            localization_matrix[i,index]&localization_matrix[j,index]
                            ) / 
                            sum(index)
                        )
                }
        })  })
rownames(jacardz) <- rownames(localization_matrix)
colnames(jacardz) <- rownames(localization_matrix)

long_jacardz <- jacardz %>% 
    {bind_cols(ORF=rownames(.),as.tibble(.))} %>% 
    pivot_longer(-ORF) %>% filter(!is.na(value)) %>%
    rename(ORF1=ORF,ORF2=name,Colocalization=value) 

dbWriteTable(conn=db_ld,value=as.data.frame(long_jacardz),name="localization_jacard")

localization %>% 
    rename(YORF=ORF) %>%
    ungroup() %>%
    select(-`Gene name`) %>%
    pivot_wider(names_from="name",values_from="value") %>%
    {dbWriteTable(conn=db_ld,value=as.data.frame(.),name="localizations")}

dbDisconnect(db_ld)
