library(tidyverse)
library(egg)

wdtweak <- ""

# wdtweak <- "../"

stats <- readRDS(str_c(wdtweak,"tmp/all_environments_scores_and_stats.RData"))

tecan_datar <- read_csv(str_c(wdtweak,"tmp/tecan_validation_statistics.csv")) 

#
#
#
#
#

validation_datar <- 
    left_join(tecan_datar,datar%>%filter(filenames=="DMSO"),by="PPI") %>% 
    group_by(PPI,`p-value`,one_tail_Q,statistic) %>%
    summarize(
        MeanFitness=mean(Fitness),MinFitness=min(Fitness),
        MaxFitness=max(Fitness),
        MeanFitness_error=mean(Fitness_estimaion_error),
        MinFitness_error=min(Fitness_estimaion_error),
        MaxFitness_error=max(Fitness_estimaion_error),
        MinLog10Counts=max(
            log10(c(Counts_G0,Counts_G3,Counts_G6,Counts_G9,Counts_G12,Counts_G18)),
            na.rm=T),
        MeanLog10Counts=mean(
            log10(c(Counts_G0,Counts_G3,Counts_G6,Counts_G9,Counts_G12,Counts_G18)),
            na.rm=T),
        MedianLog10Counts=median(
            c(Counts_G0,Counts_G3,Counts_G6,Counts_G9,Counts_G12,Counts_G18),
            na.rm=T)
        ) %>% unique() %>%
    rename(t_statistic=statistic) %>%
    mutate(
        MinLog10Counts=ifelse(is.finite(MinLog10Counts),MinLog10Counts,NA),
        MedianLog10Counts=ifelse(is.finite(MedianLog10Counts),MedianLog10Counts,NA),
        MeanLog10Counts=ifelse(is.finite(MeanLog10Counts),MeanLog10Counts,NA)
        )

lm(data=validation_datar[,-c(1:3)],t_statistic~.) %>% summary()

lm(data=validation_datar[,-1],
    one_tail_Q~MeanFitness+MeanFitness_error+
        MedianLog10Counts+MinLog10Counts+MeanLog10Counts
    ) %>% summary()

glm(data=validation_datar[,-1],
    I(one_tail_Q<0.05)~MeanFitness+MeanFitness_error+
        MedianLog10Counts+MinLog10Counts+MeanLog10Counts,
    family=binomial(link="logit")
    ) %>% summary()

glm(data=validation_datar[,-1],
    I(one_tail_Q<0.05)~MeanFitness+MeanFitness_error+MeanLog10Counts,
    family=binomial(link="logit")
    ) %>% summary()

lm(data=validation_datar[,-1],
    t_statistic~MeanFitness+MeanFitness_error+
        MedianLog10Counts+MinLog10Counts+MeanLog10Counts,
    ) %>% summary()

lm(data=validation_datar[,-1],
    t_statistic~MeanFitness+MeanFitness_error+MeanLog10Counts
    ) %>% summary()

lm(data=validation_datar[,-1],
    t_statistic~log10(MeanFitness)+log10(MeanFitness_error)+MeanLog10Counts
    ) %>% summary()

g <- ggplot(validation_datar)+theme_bw()+geom_point(size=0.5)

aplot <- ggarrange(
    g+aes(x=MeanFitness,y=t_statistic)+stat_smooth(method="lm")+
        annotate(geom="label",
            label=signif(cor.test(validation_datar$t_statistic,validation_datar$MeanFitness,use="complete.obs",method="spearman")$estimate,4),
            x=0.8,y=70)
    ,
    g+aes(x=MeanFitness_error,y=t_statistic)+stat_smooth(method="lm")+
        annotate(geom="label",
            label=signif(cor.test(validation_datar$t_statistic,validation_datar$MeanFitness_error,use="complete.obs",method="spearman")$estimate,4),
            x=7.5,y=70)
    ,
    g+aes(x=MeanLog10Counts,y=t_statistic)+stat_smooth(method="lm")+
        annotate(geom="label",
            label=signif(cor.test(validation_datar$t_statistic,validation_datar$MeanLog10Counts,use="complete.obs",method="spearman")$estimate,4),
            x=4.0,y=70)
    ,
    g+aes(x=MeanLog10Counts,y=MeanFitness_error)+stat_smooth(method="lm")+
        annotate(geom="label",
            label=signif(cor.test(validation_datar$MeanFitness_error,validation_datar$MeanLog10Counts,use="complete.obs",method="spearman")$estimate,4),
            x=3.0,y=10)
    ,
    nrow=2
    )

ggsave("output/validations_statistic_against_properties.png",
    aplot,width=6,height=6)

#
#
#
#
#

validation_rates <- 
    tecan_datar %>% 
    filter(grepl("^Y",PPI)) %>%
    left_join(
        stats %>% mutate(PPI=str_c(ORF1,ORF2,sep="_")) %>% 
            filter(Positive==1) %>%
            select(PPI,filenames,Positive)
        ,
        by="PPI"
        ) %>%
    mutate(passes=one_tail_Q < 0.05) %>%
    group_by(PPI,passes) %>%
    summarize(environments=length(unique(filenames))) %>%
    group_by(environments) %>%
    summarize(tab_passes=list(tibble(
                validate=names(table(passes)),
                count_of=table(passes)))
        ) %>% 
    unnest() %>% 
    spread(validate,count_of) %>%
    mutate(ratio=`TRUE`/(`TRUE`+`FALSE`))
validation_rates

ppiz <- stats %>% 
    filter(Positive==1) %>%
    mutate(PPI=str_c(ORF1,ORF2,sep="_")) %>% 
    filter(grepl("^Y",PPI)) %>%
    group_by(PPI) %>%
    summarize(environments=length(unique(filenames))) %>%
    group_by(environments) %>%
    summarize(count_of=length(unique(PPI)))
ppiz 

g <- full_join(validation_rates,ppiz,by="environments") %>%
    mutate(adjusted_count=count_of*ratio) %>%
    gather(type,value,count_of,adjusted_count) %>%
    mutate(type=factor(type,levels=c("count_of","adjusted_count"))) %>%
    ggplot()+theme_bw()+
    geom_bar(stat="identity",position="dodge")+
    aes(x=factor(environments),fill=type,y=value)+
    scale_fill_discrete("",
        labels=c(`adjusted_count`="Count scaled\nby validation",`count_of`="Raw count")
        )+
    ylab("Count")+xlab("Environments observed in")
g

ggsave(str_c(wdtweak,"output/ppi_counts_adjusted_by_validation.png"),
    g,
    width=6,height=4)

full_join(validation_rates,ppiz,by="environments") %>%
    mutate(adjusted_count=count_of*ratio) %>%
    pull(adjusted_count) %>% sum
