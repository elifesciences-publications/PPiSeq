library(tidyverse)
library(quantreg)
library(vegan)
library(ggplot2)
library(DBI)
library(permute)
library(parallel)                                                                         

db_bd <- dbConnect(RSQLite::SQLite(),"scratch.sqlite")  
dbExecute(db_bd,"ATTACH DATABASE 'fitness.sqlite' AS fitness") 
all_pos <- dbGetQuery(db_bd,"
    SELECT PPI,ORF1,ORF2,Experiment 
    FROM fitness WHERE Positive==1"
    ) %>%
    rowwise() %>%
    mutate(ppi_sort_name=paste0(sort(c(ORF1,ORF2)),collapse="_")) 
ppi_env <- all_pos %>% 
    select(ppi_sort_name,Experiment) %>%
    filter(!Experiment%in%c("SD","SD2")) %>%
    filter(!grepl("os",ppi_sort_name)) %>%
    filter(!grepl("eg",ppi_sort_name)) %>%
    filter(grepl("^Y",ppi_sort_name)) %>%
    distinct()
whole_ppi_list <- ppi_env %>% group_by(Experiment) %>% 
    summarize(ppis=list(ppi_sort_name)) %>%
    {setNames(.$ppis,.$Experiment)}
str(whole_ppi_list)

ppi_accumulations <- function(
        a_ppi_list,
        a_ppi_env=ppi_env,
        limit=NULL,
        n_cores=1,
        calc_the_vegan_specaccum=T,
        calc_the_vegan_specpool=T,
        calc_the_resampled_accumulation=T
    ) {
    if (calc_the_resampled_accumulation) {
        print(paste("resampling environments to get the sums, using",
                n_cores,"cores..."))
        list_of_permutations <- 
            apply(
                rbind(1:(length(a_ppi_list)),
                    permute::allPerms(1:(length(a_ppi_list)),how(maxperm = 999999))
                    ),
                1,function(x){list(x)}
            )
        resampled_accumulation <- 
            mclapply(
                X=if(is.null(limit)) {
                    list_of_permutations
                } else {
                    sample(list_of_permutations,limit)
                },
                FUN=function(x){
                    x <- unlist(x)
                    set <- vector("list",length(x))
                    set[[1]] <- a_ppi_list[x[1]][[1]]
                    for ( i in 2:length(x) ) {
                        set[[i]] <- union(set[[i-1]],a_ppi_list[x[i]][[1]])
                    }
                    return(sapply(set,length))
                },
                mc.cores=n_cores
            ) %>% 
            tibble(resultz=.) %>%
            mutate(idz=1:nrow(.)) %>%
            {suppressMessages(unnest_wider(.,resultz))} %>%
            pivot_longer(names_to="name",values_to="ppis",-idz) %>%
            mutate(samplez=str_replace(name,"\\.\\.\\.","")) %>%
            select(-name) 
        print("deduplicating the accumulation results...")
        dedup_resampled_accumulation <- resampled_accumulation %>%
            group_by(ppis,samplez) %>%
            summarize(weight=length(idz))
        print("summarizing the accumulation results...")
        summary_resampled_accumulation <- resampled_accumulation %>%
            group_by(samplez) %>%
            summarize(quantilez=list(quantile(ppis,na.rm=T))) %>%
            mutate(quantilez=map(quantilez,function(x){
                        return(tibble(value=x,quantz=names(x)))
                    })
            ) %>% unnest(quantilez) %>% 
            pivot_wider(names_from="quantz",values_from="value")
    } else {
        resampled_accumulation <- NULL
        dedup_resampled_accumulation <- NULL
        summary_resampled_accumulation <- NULL
    }
    if (calc_the_vegan_specpool | calc_the_vegan_specaccum) {
        print("creating a matrix for vegan to use...")
        all_the_ppis <- unique(unlist(a_ppi_list))
        ppi_matrix <- ppi_env %>% 
            rename(PPI=ppi_sort_name) %>% 
            rowwise() %>%
            filter(PPI %in% all_the_ppis) %>%
            mutate(code=1) %>%
            ungroup() %>% distinct() %>%
            select(PPI,Experiment,code) %>%
            arrange(PPI) %>%
            pivot_wider(names_from="Experiment",values_from="code",
                values_fill=list(code=0)
            ) %>%
            {matrix(t(.[,-1]),nrow=length(names(.))-1,
                dimnames=list(names(.)[-1],pull(.,PPI)))}
        rownames(ppi_matrix) <- NULL
    } else {
        ppi_matrix <- NULL
    }
    if (calc_the_vegan_specaccum) {
        print("calculating species accumulation curve...")
        species_accumulation <- specaccum(ppi_matrix,ci=1.96)
        gc()
    } else {
        species_accumulation <- NULL
    }
    if (calc_the_vegan_specpool) {
        print("calculating total species estimate...")
        species_pool <- specpool(ppi_matrix)
        gc()
    } else {
        species_pool <- NULL
    }
    return(
        list(
            species_pool=species_pool,
            species_accumulation=species_accumulation,
            ppi_matrix=ppi_matrix,
            resampled_accumulation=resampled_accumulation,
            dedup_resampled_accumulation=dedup_resampled_accumulation,
            summary_resampled_accumulation=summary_resampled_accumulation,
            idz=runif(1)
        )
    )
}




#####
#####
##### calc them
#####
#####

whole_list_accumulation <- ppi_accumulations(whole_ppi_list,
        a_ppi_env=ppi_env,
        n_cores=ifelse(detectCores()==16,14,6),
        limit=10000,
        calc_the_vegan_specaccum=T,
        calc_the_vegan_specpool=T,
        calc_the_resampled_accumulation=F 
    )

saveRDS(whole_list_accumulation,"raw_accumulations.RData")

rm(whole_list_accumulation)
gc()

#####
#####
##### re-write below to take the ppi by ppi resampling, do that 10 times
#####
#####

validation_rates_per_measure <- 
    read_csv("Validated_PPI_environment_count_summary_SD_merge_filter.csv") %>%
    select(-environment_number) %>%
    pivot_longer(names_to="Experiment",values_to="validation_pred",-PPI) %>%
    filter(!is.na(validation_pred)) %>%
    separate(PPI,into=c("ORF1","ORF2"),sep="_") %>%
    rowwise() %>%
    mutate(ppi_sort_name=paste0(sort(c(ORF1,ORF2)),collapse="_")) 

set.seed(1234)
resampled_list_accumulation <- 
    mclapply(
        c(1:32),
        function(x){
            validation_rates_per_measure %>% 
                ungroup() %>%
                filter(runif(nrow(.)) < validation_pred) %>%
                select(-validation_pred) %>%
                group_by(Experiment) %>% 
                summarize(ppis=list(ppi_sort_name)) %>%
                {setNames(.$ppis,.$Experiment)} %>% 
                ppi_accumulations(
                    limit=10000,
                    calc_the_resampled_accumulation=T,
                    calc_the_vegan_specaccum=T,
                    calc_the_vegan_specpool=T,
                ) %>%
                return()
        }
        ,mc.cores=ifelse(detectCores()==16,8,2)
        ,mc.cleanup=T
    )

saveRDS(resampled_list_accumulation,"thirtytwo_resampled_accumulations.RData")

rm(resampled_list_accumulation)
gc()

#####
#####
##### Doing just one resample with all of this
#####
#####


set.seed(1234)
one_resampled_list_accumulation <- 
    validation_rates_per_measure %>% 
        ungroup() %>%
        filter(runif(nrow(.)) < validation_pred) %>%
        select(-validation_pred) %>%
        group_by(Experiment) %>% 
        summarize(ppis=list(ppi_sort_name)) %>%
        {setNames(.$ppis,.$Experiment)} %>% 
        ppi_accumulations(
            limit=10000,
            n_cores=ifelse(detectCores()==16,14,6),
            calc_the_resampled_accumulation=T,
            calc_the_vegan_specaccum=T,
            calc_the_vegan_specpool=T
        ) 

saveRDS(one_resampled_list_accumulation,"one_resampled_accumulations.RData")

rm(one_resampled_list_accumulation)
gc()
