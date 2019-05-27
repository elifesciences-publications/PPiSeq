library(tidyverse)

# Reading in the processed csv calls

# Means and statistics
meanz <- tibble(
        filenames=list.files("data",pattern=".*mean_fitness.*\\.csv",
            full.names=T)
        ) %>% 
    mutate(rawfiles=map(filenames,function(x){print(x);read_csv(x)})) %>% 
    mutate(
        filenames=sub("data\\/(.*)_mean_fitness.*","\\1",filenames)
        ) %>%
    unnest() %>%
    mutate(PPI=ifelse(
            grepl("Neg",PPI)|grepl("Pos",PPI)|grepl("negative",PPI),
            paste0(gsub("_","-",PPI),"_neg"),
            PPI)
        ) %>% 
    separate(PPI,into=c("ORF1","ORF2"),sep="_")
saveRDS(meanz,"tmp/all_environments_scores_and_stats.RData")
rm(meanz)
gc()

