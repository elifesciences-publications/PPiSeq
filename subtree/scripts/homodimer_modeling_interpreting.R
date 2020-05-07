## title: checking homodimer modeling result, determining the ones we think are explained
## author: darach

# Libraries

library(tidyverse)
library(magrittr)

# A palette for the plots, from
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

within_ppi <- 
    read_csv("homodimer_all_models.csv") %>%
    nest(data=c(Positive,Experiment,Normalized_Fitness,
        homodimer_score_1,homodimer_score_2,
        min_homodimer_score,max_homodimer_score)
        ) 

# Checking out the models

# How'd we do? First, p-value histograms.

within_ppi %>% 
    filter(grepl("homodimer",term)) %>% 
    group_by(name)%>%
    summarize(z=list(hist(pvalue,100,plot=F))) %>%
    do(
        pmap_dfr(list(x=.$z,name=.$name),function(x,name){
            plot(x,main=paste0(name," model"))
            abline(v=0.05,col="red")})
        )

#Okay. Let's see them all. Here's all the observations as points, then line
#is fit just by `stat_smooth` (which should be identical), and it's faint so you
#can just see the "grain" of what's happening.
#
#First additive model

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="additive")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=homodimer_score_1+homodimer_score_2,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1)

# Multiplicative model

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="multiplicative")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=homodimer_score_1*homodimer_score_2,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1)

# sqrt multiplicative

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="sqrtmult")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1,
        formula="y~x+1")

# sqrt multiplicative with no intercept

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="sqrtmult_nointer")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1,
        formula="y~x+0")

# log of square root of exponentiated sum (pretty similar to sqrt multiplicative

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="sqrtmultexp")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=log(sqrt(exp(homodimer_score_1+homodimer_score_2))),
        y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1,
        formula="y~x+1")

#Minimum average score of the pair

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="min_determines")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=min_homodimer_score,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1)

#Maximum average score of the pair

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="max_determines")%>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    aes(x=max_homodimer_score,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.1)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.03,alpha=0.1)

#Okay, and now the same but just for the significant ones. Printed for
#each one is the total number of combinations tested (same for all, should be),
#and then the number with pval<0.05 significance.

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="additive") %T>% 
    {print("additive model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=homodimer_score_1+homodimer_score_2,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1)

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="multiplicative") %T>% 
    {print("multiplicative model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=homodimer_score_1*homodimer_score_2,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1)

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="sqrtmult") %T>% 
    {print("sqrtmult model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+1")

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="sqrtmult_nointer") %T>% 
    {print("sqrtmult_nointer model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+0")

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="sqrtmultexp") %T>% 
    {print("sqrtmultexp model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=log(sqrt(exp(homodimer_score_1+homodimer_score_2))),
        y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+1")

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="min_determines") %T>% 
    {print("min_determines model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=min_homodimer_score,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1)

within_ppi %>% filter(grepl("homodimer",term)) %>%
    filter(name=="max_determines") %T>% 
    {print("max_determines model")} %T>% {print(nrow(.))} %>%
    filter(qvalz<0.05) %T>% {print(nrow(.))} %>%
    unite("PPI",ORF1,ORF2) %>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=max_homodimer_score,y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1)

#Okay, so which model explains a lot of the data for a lot of PPIs? How does the
#R^2 values look?

within_ppi %>% filter(grepl("homodimer",term)) %>%
    unite("PPI",ORF1,ORF2) %>%
    ggplot()+
    facet_grid(name~.)+
    theme_bw()+
    aes(x=r_squared)+
    geom_histogram(bins=50)

within_ppi %>% filter(grepl("homodimer",term)) %>%
    unite("PPI",ORF1,ORF2) %>%
    ggplot()+
    theme_bw()+
    aes(x=name,y=r_squared)+
    geom_boxplot()

#Well, it looks like the max/min models don't do well, which makes sense. I 
#think the sqrt multiplicative model does the best, with intercept.

#Okay, well, let's plot out 12 examples of actual data, first against the
#additive data/model and then against the multiplicative. 
#Only for significant ones.
#Note the scales are free

set.seed(123)
within_ppi %>% filter(grepl("homodimer",term)) %>% 
    filter(name=="additive") %>%
    filter(qvalz<0.05) %>%
    unite("PPI",ORF1,ORF2) %>%
    sample_n(12)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    ggtitle("additive model")+
    aes(x=homodimer_score_1+homodimer_score_2,y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1)

set.seed(123)
within_ppi %>% filter(grepl("homodimer",term)) %>% 
    filter(name=="multiplicative") %>%
    filter(qvalz<0.05) %>%
    unite("PPI",ORF1,ORF2) %>%
    sample_n(12)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    ggtitle("multiplicative model")+
    aes(x=homodimer_score_1*homodimer_score_2,y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1)

set.seed(123)
within_ppi %>% filter(grepl("homodimer",term)) %>% 
    filter(name=="sqrtmult") %>%
    filter(qvalz<0.05) %>%
    unite("PPI",ORF1,ORF2) %>%
    sample_n(12)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    ggtitle("sqrtmult model")+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+1")

set.seed(123)
within_ppi %>% filter(grepl("homodimer",term)) %>% 
    filter(name=="sqrtmult_nointer") %>%
    filter(pvalue<0.05) %>%
    unite("PPI",ORF1,ORF2) %>%
    sample_n(12)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    ggtitle("sqrtmult_nointer model")+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    ggtitle("THIS ONE IS PVAL<0.05 because qvalues are bonkers")+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+0")

set.seed(123)
within_ppi %>% filter(grepl("homodimer",term)) %>% 
    filter(name=="sqrtmultexp") %>%
    filter(qvalz<0.05) %>%
    unite("PPI",ORF1,ORF2) %>%
    sample_n(12)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    ggtitle("sqrtmultexp model")+
    aes(x=log(sqrt(exp(homodimer_score_1+homodimer_score_2))),y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+1")

# Selecting set we really believe

# How does low observation affect our confidence here?
# I've fiddled with the thresolds. If we require four or more observations,
# it looks pretty good! 

within_ppi %>%
    filter(name=="sqrtmult",term!="(Intercept)") %>% 
    mutate(
        nobs=unlist(map(data,function(x){length(x$Positive)})),
        npos=unlist(map(data,function(x){sum(x$Positive)}))
        ) %T>%
    {print(hist(.$nobs))} %T>%
    {print(hist(.$npos))} %>%
    filter(nobs>=4,npos>=1) %>%
    arrange(-r_squared) %>%
    head(16) %>%
    unnest(data) %>%
    ggplot()+theme_bw()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness)+
    stat_smooth(method="lm",col="grey")+
    geom_point(aes(col=Positive))+
    facet_wrap(~ORF1+ORF2)

within_ppi %>%
    filter(name=="sqrtmult",term!="(Intercept)") %>% 
    mutate(
        nobs=unlist(map(data,function(x){length(x$Positive)})),
        npos=unlist(map(data,function(x){sum(x$Positive)}))
        ) %>%
    filter(nobs>=4,npos>=1) %>%
    arrange(r_squared) %>%
    head(16) %>%
    unnest(data) %>%
    ggplot()+theme_bw()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),y=Normalized_Fitness)+
    stat_smooth(method="lm",col="grey")+
    geom_point(aes(col=Positive))+
    facet_wrap(~ORF1+ORF2)

# Applying those, also
# The no-intercept model is the proper one to use , but everything fits that 
# because that's then the only term and there's a trend there. Instead, I'll 
# take the fits with the linear model and use that to filter out PPIs where 
# the intercept is significant (p-value single test).

# My metric of a binary good fit is this, then for a quantative metric of
# fit will be R^2 to no intercept sqrt multiplicative model.
# Actually, this loses the sign of the slope, so instead I'll use the
# pearson correlation

full_set <- within_ppi %>%
    filter(name=="sqrtmult") %>%
    mutate(value=ifelse(term=="(Intercept)",pvalue,qvalz)) %>% 
    select(ORF1,ORF2,name,term,value) %>% 
    pivot_wider(names_from="term") %>%
    mutate(Explained=
            (`I(sqrt(homodimer_score_1 * homodimer_score_2))` < 0.05) &
            (`(Intercept)` > 0.05)
        ) %>%
    select(ORF1,ORF2,name,Explained) %>% 
    distinct() %>%
    left_join(within_ppi,
        by=c("ORF1","ORF2","name")
        ) %>%
    filter(term!="(Intercept)") %>%
    unite("PPI",ORF1,ORF2,remove=F) %>%
    mutate(
        nobs=unlist(map(data,function(x){length(x$Positive)})),
        npos=unlist(map(data,function(x){sum(x$Positive)})),
        pearson=unlist(map(data,function(x){
                    cor(method="pearson",
                        x$Normalized_Fitness,
                        (x$homodimer_score_1*x$homodimer_score_2)^0.5
                        )
                    })),
        ) %>%
    filter(nobs>=4,npos>=1)

# Write this out:
full_set %>% 
    unnest(data) %>%
    write_csv(path="explained_ppis_modeling.csv")

# How do Explained or not PPIs look?

full_set %>%
    filter(qvalz<0.05) %>%
    unnest(data) %>%
    ggplot()+
    theme_classic()+
    facet_wrap(~Explained) +
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),
        y=Normalized_Fitness,group=PPI,col=Experiment)+
    geom_point(size=0.5,alpha=0.5)+
    guides(col=guide_legend(override.aes=list(alpha=1)))+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=1.0,
        formula="y~x+1")

# Okay, let's pick 12 and see

full_set %>%
    filter(Explained) %>% 
    sample_n(12)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),
        y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+1")

# And some that aren't explained?

full_set %>%
    filter(Explained) %>% 
    sample_n(6)%>%
    unnest(data) %>%
    ggplot()+
    theme_bw()+
    aes(x=sqrt(homodimer_score_1*homodimer_score_2),
        y=Normalized_Fitness,col=Experiment)+
    facet_wrap(~PPI,scales="free",ncol=3)+
    geom_point()+
    stat_smooth(method="lm",se=F,color="black",size=0.1,alpha=0.1,
        formula="y~x+1")

# How many are in each category?

full_set %>%
    group_by(qvalz<0.05,Explained) %>%
    summarize(count=length(PPI),prop=length(PPI)/nrow(.))


#Okay. So, of the explained ones, how are different ORFs represented?
#How many do we see of each ORF?

full_set %>% select(ORF1,ORF2,r_squared,Explained) %>% 
    gather(which,YORF,ORF1,ORF2) %>% group_by(YORF) %>% 
    summarize(count=length(YORF)) %>% 
    ggplot()+theme_bw()+aes(x=count)+geom_histogram(stat="count") 

#Long-tail distribution...

#How might the fits be distributed by the observability of the PPI, and/or the
#number of fits?

full_set %>% 
    select(ORF1,ORF2,r_squared,pearson,Explained,nobs) %>% 
    gather(which,YORF,ORF1,ORF2) %>% 
    group_by(YORF) %>% 
    mutate(count=length(YORF)) %>% 
    arrange(count) %>% 
    ungroup() %>%
    mutate(YORF=factor(YORF,levels=unique(YORF))) %>% 
    group_by(YORF) %T>% 
    {
    print(
    ggplot(select(.,YORF,nobs) %>% distinct())+theme_bw()+
    aes(x=as.numeric(YORF),y=nobs)+
    geom_point()+
    xlab("Ranked from few PPIs per ORF to many on right")+
    ylab("Observed in this many environments")+
    theme(axis.text.x=element_blank())
    )
    } %>%
    {
    print(
    ggplot(.)+theme_bw()+
    aes(x=as.numeric(YORF),y=pearson)+
    geom_dotplot(binwidth=0.01,binaxis="y",stackdir="center",aes(group=YORF))+
    geom_line(
        data=summarize(.,mean_pearson=mean(pearson))%>%ungroup(),
        aes(y=mean_pearson)
        )+
    xlab("Ranked from few PPIs per ORF to many on right")+
    ylab("Pearson of a PPI or mean pearson per ORF")+
    theme(axis.text.x=element_blank())
    )
    }

#First, number of PPIs participating in compared to number of environments its
#observed in. It would appear that the ORFs where there's only one or two PPIs
#aren't just there because of noisy rare observation.

#Then, it appears that the more PPIs the wider spread and lower mean of R^2 that
#seem to appear.

#What about just within the explained ones?

full_set %>% 
    filter(Explained) %>%
    select(ORF1,ORF2,r_squared,pearson,Explained,nobs) %>% 
    gather(which,YORF,ORF1,ORF2) %>% 
    group_by(YORF) %>% 
    mutate(count=length(YORF)) %>% 
    arrange(count) %>% 
    ungroup() %>%
    mutate(YORF=factor(YORF,levels=unique(YORF))) %>% 
    group_by(YORF) %T>% 
    {
    print(
    ggplot(select(.,YORF,nobs) %>% distinct())+theme_bw()+
    aes(x=as.numeric(YORF),y=nobs)+
    geom_point()+
    xlab("Ranked from few PPIs per ORF to many on right")+
    ylab("Observed in this many environments")+
    theme(axis.text.x=element_blank())
    )
    } %>%
    {
    print(
    ggplot(.)+theme_bw()+
    aes(x=as.numeric(YORF),y=pearson)+
    geom_dotplot(binwidth=0.005,binaxis="y",stackdir="center",aes(group=YORF))+
    geom_line(
        data=summarize(.,mean_pearson=mean(pearson))%>%ungroup(),
        aes(y=mean_pearson)
        )+
    xlab("Ranked from few PPIs per ORF to many on right")+
    ylab("R^2 of a PPI or mean R^2 per ORF")+
    theme(axis.text.x=element_blank())
    )
    }


