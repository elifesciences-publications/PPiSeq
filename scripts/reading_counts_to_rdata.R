library(tidyverse)

# Reading in the processed csv calls

# Counts
countz <- tibble(
        filenames=list.files("data",pattern=".*fitness_counts\\.csv",
            full.names=T)
        ) %>% 
    mutate(rawfiles=map(filenames,function(x){print(x);read_csv(x)})) %>% 
    mutate(
        filenames=sub("data\\/(.*)_PPI_.*","\\1",filenames)
        ) %>%
    unnest() %>%
    rename(Fitness_estimation_error=Fitness_estimaion_error) 

saveRDS(countz,"tmp/all_environments_counts.RData")

