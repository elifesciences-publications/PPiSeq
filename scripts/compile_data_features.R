## title: compiling names and go terms into tables
## author: darach

library(tidyverse)
library(DBI)

db_fe <- dbConnect(RSQLite::SQLite(),"features.sqlite")

read.delim("geneFeatures_022415_EK.txt",
        sep="\t") %>%
    rename(YORF=X) %>%
    {dbWriteTable(conn=db_fe,value=as.data.frame(.),name="costanza2016")}

read_csv("sgd_homologues_191027.csv",col_names=F) %>%
    rename(ORF1=X1,ORF2=X2) %>%
    mutate(homolog=1) %>%
    {dbWriteTable(conn=db_fe,value=as.data.frame(.),name="homologs")}


read_csv("marchant2019_tableS3.csv",col_names=T) %>%
    rowwise() %>%
    mutate(
        namez=list(
            list(
                data.frame(ORF1=P1,ORF2=P2),
                data.frame(ORF1=P2,ORF2=P1)
        )   )   ) %>% 
    select(namez,Duplication,Origin.of.WGDs) %>%
    unnest(namez) %>% unnest(namez) %>% 
    {dbWriteTable(conn=db_fe,value=as.data.frame(.),name="marchant_paralogs")}

read_csv("ygob_ohnologs_191028.csv",col_names=T) %>%
    rowwise() %>%
    mutate(
        namez=list(
            list(
                data.frame(ORF1=`Gene 1`,ORF2=`Gene 2`),
                data.frame(ORF1=`Gene 2`,ORF2=`Gene 1`)
        )   )   ) %>% 
    mutate(ohnolog=1) %>%
    select(namez,ohnolog) %>%
    unnest(namez) %>% unnest(namez) %>% 
    {dbWriteTable(conn=db_fe,value=as.data.frame(.),name="ohnologs")}

read_csv("ho2018unified_ProteinAbundances.csv") %>%
    rename(YORF=`Systematic Name`,
        MeanMol=`Mean molecules per cell`,
        MedianMol=`Median molecules per cell`
        ) %>%
    {dbWriteTable(conn=db_fe,value=as.data.frame(.),name="ho2018")}

read_delim("BIOGRID-PTM-3.5.177_ptm_yeast.txt",
        col_types="---c----ccc-------",delim="\t") %>% 
    rename(YORF=`Systematic Name`) %>% 
    select(YORF,`Post Translational Modification`,Residue) %>%
    pivot_wider(names_from="Post Translational Modification",
        values_from="Residue",values_fn=list(Residue=length),
        values_fill=list(Residue=0)
        ) %>%
    {dbWriteTable(conn=db_fe,value=as.data.frame(.),name="biogrid_ptms")}

dbDisconnect(db_fe)
