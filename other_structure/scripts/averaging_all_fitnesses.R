library(tidyverse)

allMeanFitness <- readRDS("tmp/all_environments_counts.RData") %>%
    filter(Number_of_Barcodes > 2) %>% 
    group_by(filenames,PPI) %>% 
    summarize(MeanFitness=mean(Fitness,na.rm=T)) %>% 
    ungroup() %>% mutate(filenames=sub("data\\/(.*)_PPI_.*","\\1",filenames))

saveRDS(allMeanFitness,
    file="tmp/all_environments_all_biological_ppi_mean_fitness.RData")

