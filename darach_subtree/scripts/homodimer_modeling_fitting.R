## title: Modeling how the homodimers explain the heterodimers
## author: darach

# Here we aim to model how well the PPISeq signal of each protein's homodimer 
# explains the heterodimer of the two proteins. This happens for quite a few 
# proteins in yeast. So, this is the closest with this dataset we get to 
# controlling or accounting for abundance changes.

# This file gets that data in, munges it around, fits a couple of models,
# then saves it. We interpret the fits in the next file...

# Libraries

library(tidyverse)
library(DBI)
library(magrittr)

db_sc <- dbConnect(RSQLite::SQLite(),"fitness.sqlite")

# TEMP
# STRICT THRESHOLDS
#   So the ultra-strict thresholds were added late in the game to double check,
#   and I lazily just inserted the below snippet and the one marked in the 
#   chunk below. We saw that it didn't change much except that there were a
#   lot less, but the patterns of correlation of abundance with PPI signal was
#   still increasing with environment number so ...
#   Yes it's sloppy engineering, but this is where we're at.
#
#strict_calls <- read_csv("PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv") %>%
#    rename(PPI=X) %>% rename(`16C`=X16C) %>% 
#    pivot_longer(names_to="Experiment",values_to="Positive",3:11)
# TEMP

normed_datar <- dbGetQuery(db_sc,'SELECT * FROM fitness') %>% as_tibble() %>%
    filter(!grepl("SD2?$",Experiment)) %>% ungroup() %>% 
# STRICT THRESHOLDS
#   uncomment the below to do stricter thresholds ("higher confidence" set)
    #select(-Positive) %>% left_join(strict_calls,by=c("PPI","Experiment")) %>% 
    ungroup() %>% 
    mutate(Positive=ifelse(is.na(Positive),0,Positive))
normed_datar

# Selecting out the homo and hetero dimers from in that normalized data

# Note that we're picking homodimers seen at least four times, with at least
# one positive interaction

homodimers <- normed_datar %>% 
    filter(ORF1==ORF2) %>%
    mutate(sort_id=unlist(
            pmap(
                list(o1=ORF1,o2=ORF2),
                function(o1,o2){
                    paste0(sort(c(o1,o2)),collapse="_")
                    }
                )
            )
        ) %>% 
    group_by(sort_id) %>% 
    filter(length(Positive)>=4,sum(Positive)>=1) %>% 
    mutate(YORF=ORF1) %>% 
    select(Experiment,YORF,Normalized_Fitness) %>%
    rename(homodimer_score=Normalized_Fitness)
homodimers

# Then we use joins to pull out all heterodimers with the same criteria

homodimer_datar <- normed_datar %>% 
    select(Experiment,ORF1,ORF2,Normalized_Fitness,Positive) %>%
    left_join( homodimers %>% ungroup() %>% select(-sort_id) %>%
            rename(ORF1=YORF,homodimer_score_1=homodimer_score),
        by=c("ORF1","Experiment")
        ) %>%
    left_join( homodimers %>% ungroup() %>% select(-sort_id) %>%
            rename(ORF2=YORF,homodimer_score_2=homodimer_score),
        by=c("ORF2","Experiment")
        ) %>%
    filter(!( is.na(homodimer_score_1)|is.na(homodimer_score_2) )) %>%
    filter(ORF1!=ORF2) %>% 
    group_by(ORF1,ORF2) %>% 
    nest(data=c(Experiment,Normalized_Fitness, Positive, 
        homodimer_score_1,homodimer_score_2) 
        ) %>%
    mutate(sort_id=unlist(
            pmap(
                list(o1=ORF1,o2=ORF2),
                function(o1,o2){
                    paste0(sort(c(o1,o2)),collapse="_")
                    }
                )
            )
        ) 


amongst_homodimers <- homodimer_datar %>% 
    unnest(data) %>%
    group_by(ORF1,ORF2) %>% 
    filter(length(Positive)>=4,sum(Positive)>=1) 

amongst_homodimers_pooled_tags <- homodimer_datar %>% 
    group_by(sort_id) %>%
    summarize(data=list(bind_rows((data)))) %>%
    separate(sort_id,into=c("ORF1","ORF2"),sep="_",remove=F) %>%
    unnest(data) %>%
    group_by(ORF1,ORF2) %>% 
    filter(length(Positive)>=4,sum(Positive)>=1) 

print(head(amongst_homodimers_pooled_tags,20))

# Doing the modeling

#*#Here's six models:
#*#- Normalized_Fitness~I(homodimer_score_1+homodimer_score_2)
#*#- Normalized_Fitness~I(homodimer_score_1\*homodimer_score_2) 
#*#- Normalized_Fitness~I(sqrt(homodimer_score_1\*homodimer_score_2)) 
#*#- Normalized_Fitness~0+I(sqrt(homodimer_score_1\*homodimer_score_2)) 
#*#- Normalized_Fitness~I(log(sqrt(exp(homodimer_score_1+homodimer_score_2)))) 
#*#- Normalized_Fitness~max_homodimer_score 
#*#- Normalized_Fitness~min_homodimer_score 

#`0+` means no intercept, all else have implicit intercept fit. All are linear
#models.

within_ppi_pre <- amongst_homodimers_pooled_tags %>% # just taking amongst homodimers
    arrange(ORF1,ORF2) %>% # sorting so we can spot duplicates easily, just cuz
    group_by(ORF1,ORF2) %>% # for each PPI
    mutate(
        min_homodimer_score=ifelse( # we test if
            rep(
                mean(homodimer_score_1)<=mean(homodimer_score_2), # is 1 <= 2?
                length(homodimer_score_1) # this so it's rep'd the length needed
                ),
            homodimer_score_1,  # 1 <= 2, so is the minimum
            homodimer_score_2), # otherwise
        max_homodimer_score=ifelse(
            rep(
                mean(homodimer_score_1)>=mean(homodimer_score_2), # the converse
                length(homodimer_score_1)
                ),
            homodimer_score_1,  # 1 >= 2, so is max
            homodimer_score_2)  # otherwise
        ) %>%
    nest(data=c(Experiment,Normalized_Fitness, Positive, 
        homodimer_score_1,homodimer_score_2, # nesting to pack data up into rows
        min_homodimer_score,max_homodimer_score)
        ) %>%
    rowwise() %>% # and we proceed along the rows
    mutate(
        model_list=list(list(      # fitting a couple of models, so that's 
            additive=lm(data=data, #     packaged in lists of lists
                Normalized_Fitness~I(homodimer_score_1+homodimer_score_2) # additive
                ),
            multiplicative=lm(data=data,
                Normalized_Fitness~I(homodimer_score_1*homodimer_score_2) # multiplicative
                ),
            sqrtmult=lm(data=data,
                Normalized_Fitness~I(sqrt(homodimer_score_1*homodimer_score_2)) # 
                ),
            sqrtmult_nointer=lm(data=data,
                Normalized_Fitness~0+I(sqrt(homodimer_score_1*homodimer_score_2)) # 
                ),
            sqrtmultexp=lm(data=data,
                Normalized_Fitness~I(log(sqrt(exp(homodimer_score_1+homodimer_score_2)))) # 
                ),
            max_determines=lm(data=data,
                Normalized_Fitness~max_homodimer_score # max score above
                ),
            min_determines=lm(data=data,
                Normalized_Fitness~min_homodimer_score # min score above
                )
            ))
        ) %>%
    rowwise() %>% # then, we go along each row (each PPI)
    mutate(
        summaries=list(pmap(list(model=model_list,name=names(model_list)),
            function(model,name){ # we take each model, and each name of it
                tmpcoef <- summary(model)$coef # take the fit coefficients
                z <- as_tibble(tmpcoef) # munge
                names(z) <- c("estimate","stderr","tvalue","pvalue") # munge
                z <- z %>% mutate(term=rownames(tmpcoef)) # munge
                return(tibble(z=list(z),r_squared=summary(model)$r.squared,name=name))
                # and we're returning the coefficients with pvalues, and the R^2
            })) 
        ) %>% 
    select(-model_list) %>%  # drop that
    unnest(summaries) %>% unnest(summaries) %>% unnest(z) %>% 
    group_by(name,term) %>%
    mutate(qvalz=p.adjust(pvalue,method="fdr")) %>% # benjamini-hochberg
    unnest(data) %>% ungroup()


# Commented out, but here just to remember the nesting scheme I'll use later
# when I read the file back in
#within_ppi <- within_ppi_pre %>%
#    nest(data=c(Positive,Experiment,Normalized_Fitness,
#        homodimer_score_1,homodimer_score_2,
#        min_homodimer_score,max_homodimer_score)
#        ) 

write_csv(within_ppi_pre,path="homodimer_all_models.csv")
