library(tidyverse)
library(ggrepel)
library(egg)

# setwd("../")

read_tsv("data/sgd_go_slim_171013.txt",col_names=F,comment="!"
        ) -> SGDGOSlim
read_tsv("data/sgd_go_terms_171013.txt",col_names=F
        ,comment="!") %>% 
    mutate(GOID=str_c("GO:" ,str_pad(string=X1,width=7,side="left",pad="0"))
        ) -> SGDGOTerms
read_tsv("data/sgd_go_full_171013.txt",col_names=F,comment="!"
        ) -> SGDGOFull

allBiologicalPPIMeanFitness <- 
    readRDS("tmp/all_environments_all_biological_ppi_mean_fitness.RData") %>%
    filter(
        !grepl("HO:TEF1",PPI),
        !grepl("^Pos",PPI),!grepl("^Neg",PPI),
        !grepl("^positive",PPI),!grepl("^negative",PPI)
        )
allBioPPIMeansForConditions <- allBiologicalPPIMeanFitness %>% 
    spread(PPI,MeanFitness) %>% as.matrix() %>%
    { rownames(.) <- .[,1]; .[,-1] }
rm(allBiologicalPPIMeanFitness)
gc()
these_conditions <- rownames(allBioPPIMeansForConditions)
pca_allBioPPIMeansForConditions <- 
    allBioPPIMeansForConditions %>% 
    apply(.,2,as.numeric) %>%
    { prcomp(
        x=.[,
            apply(.,2,sum,na.rm=T)>0&
            apply(.,2,function(x){all(x>0.2)})&
            apply(.,2,function(x){!any(is.na(x))})]
        ,scale.=T,center=T
        ) }
#rm(allBioPPIMeansForConditions)
gc()

g <- 
    bind_cols(
        Conditions=these_conditions,
        as_tibble(pca_allBioPPIMeansForConditions$x)
        ) %>%
    ggplot()+theme_bw()+
    aes(label=Conditions)+
    geom_label_repel()
g + aes(x=PC1,y=PC2)
g + aes(x=PC3,y=PC4)
g + aes(x=PC5,y=PC6)

varz <- pca_allBioPPIMeansForConditions$sdev /
    sum(pca_allBioPPIMeansForConditions$sdev)

aplot <- ggarrange(
    g + aes(x=PC1,y=PC2)+
        xlab(str_c("PC1"," ",signif(varz[1])))+ylab(str_c("PC2"," ",signif(varz[2]))),
    g + aes(x=PC3,y=PC4)+
        xlab(str_c("PC3"," ",signif(varz[3])))+ylab(str_c("PC4"," ",signif(varz[4]))),
    g + aes(x=PC5,y=PC6)+
        xlab(str_c("PC5"," ",signif(varz[5])))+ylab(str_c("PC6"," ",signif(varz[6]))),
    g + aes(x=PC7,y=PC8)+
        xlab(str_c("PC7"," ",signif(varz[7])))+ylab(str_c("PC8"," ",signif(varz[8]))),
    nrow=2
    )

ggsave("output/prcomp_of_samples_pc12345678.svg",aplot,width=8,height=8)

load("tmp/NameIDList.RData")

term_lookup <- SGDGOSlim %>% dplyr::select(X1,X6) %>% group_by(X1) %>% 
    summarize(all_terms=list(X6)) %>%
    { setNames(.$all_terms,.$X1) }


term_to_name <- SGDGOSlim %>% dplyr::select(X6,X5) %>% 
    rename(TERM=X6,NAME=X5) %>% data.frame()

ppis_by_loading <- bind_cols(
        PPI=rownames(pca_allBioPPIMeansForConditions$rotation),
        as_tibble(pca_allBioPPIMeansForConditions$rotation)
        ) %>%
    gather(PC,Value,starts_with("PC")) 

ppis_terms <- ppis_by_loading %>% 
    dplyr::select(-PC,-Value) %>% unique() %>%
    separate(PPI,into=c("A","B"),sep="_",remove=F) %>%
    mutate(terms_for_ppi=pmap(list(A,B),function(a,b){
                (unique(c(unlist(term_lookup[a]),unlist(term_lookup[b]))))
            })) %>% 
    unnest() %>%
    dplyr::select(-A,-B) %>% 
    dplyr::rename(TERM=terms_for_ppi,GENE=PPI) %>% 
    unique() %>%
    data.frame()

library(clusterProfiler)
#library(org.Sc.sgd.db)

rm(allBioPPIMeansForConditions)
gc()

nPermutations <- 1e8
avg_rank_gsea <- ppis_by_loading %>% 
#    dplyr::filter(PC=="PC1") %>%
    group_by(PC) %>% 
    mutate(ranked=rank(-Value)) %>% 
    arrange(ranked) %>%
    summarize(raw_gsea=list(
            GSEA(gene=setNames(Value,nm=PPI),
                TERM2GENE=ppis_terms[,c("TERM","GENE")],
                TERM2NAME=term_to_name,
                minGSSize=2,maxGSSize=500,
                pAdjustMethod="fdr",
                pvalueCutoff=0.2,seed=190523,
                nPerm=nPermutations
                )
            )
        )
avg_rank_gsea$raw_gsea[[1]]@result[,1:9]

avg_rank_gsea %>% { map(.$raw_gsea,~.x@result) }

