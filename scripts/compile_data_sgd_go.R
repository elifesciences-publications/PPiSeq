## title: compiling names and go terms into sqlite tables
## author: darach

library(tidyverse)
library(DBI)

db_sg <- dbConnect(RSQLite::SQLite(),"sgd_go.sqlite")

read_tsv("go_terms_191111.tab",
        col_names=F,comment="!") %>%
    mutate(GO=paste0("GO:",sprintf("%07d",X1))) %>%
    rename(Domain=X3,GO_desc=X2,GO_long_desc=X4) %>%
    select(-X1) %>%
    {dbWriteTable(conn=db_sg,value=as.data.frame(.),name="go_terms")}

read_tsv("gene_association_191111.sgd",
        col_names=F,comment="!") %>%
    rename(
        Systematic=X2,
        Common=X3,
        GO=X5,
        Domain=X9
        ) %>%
    mutate(Systematic=unlist(map(X11,function(x){unlist(strsplit(x,"\\|"))[1]}))) %>% 
    {dbWriteTable(conn=db_sg,value=as.data.frame(.),name="sgd_gene_association")}

dbExecute(db_sg,"CREATE TABLE sgd_go_map AS 
    SELECT * FROM sgd_gene_association 
    LEFT OUTER JOIN go_terms
    ON 
        ( go_terms.GO == sgd_gene_association.GO ) 
    AND
        ( go_terms.Domain == sgd_gene_association.Domain ) ")

read_tsv("go_slim_mapping_191111.tab",
        col_names=F,comment="!") %>%
    rename(
        Systematic=X1,
        Common=X2,
        GO=X6,
        GO_desc=X5,
        Domain=X4
        ) %>%
    {dbWriteTable(conn=db_sg,value=as.data.frame(.),name="go_slim")}

read_tsv("go_protein_complex_slim_191111.tab",col_names=F) %>% 
    mutate(complex=map(X1,function(x){strsplit(x,"Component: ")[[1]][2]})) %>% 
    unnest(complex) %>% 
    mutate(genes=map(X2,function(x){strsplit(x,"\\|")[[1]][1]})) %>% 
    mutate(Systematic=map(genes,function(x){strsplit(x,"\\/")[[1]][3]})) %>% 
    unnest(Systematic) %>% select(Systematic,complex) %>%
    mutate(Complex=sub("^(.*)\\/.*?$","\\1",complex),
        ComplexID=sub("^(.*)\\/(.*?)$","\\2",complex)
        ) %>%
    select(-complex) %>% 
    {dbWriteTable(conn=db_sg,value=as.data.frame(.),name="go_complexes")}

dbDisconnect(db_sg)
