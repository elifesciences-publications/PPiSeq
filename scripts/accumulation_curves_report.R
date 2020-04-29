library(tidyverse)
library(ggplot2)

ggsave_both <- function(filestub,plot_obj,width=5,height=5) {
    message(paste0("Saving ",filestub,".png"))
    ggsave(paste0(filestub,".png"),plot_obj,width=width,height=height)
    message(paste0("Saving ",filestub,".pdf"))
    ggsave(paste0(filestub,".pdf"),plot_obj,device=cairo_pdf,width=width,height=height)
}

whole_list_accumulation <- readRDS("raw_accumulations.RData")

pdatar <- 
    tibble(
        estimates=list(whole_list_accumulation$species_pool)
        ,
        accum=list(with(whole_list_accumulation$species_accumulation,
                tibble(samples=sites,ppis=richness)
            ))
        ,
        idz=whole_list_accumulation$idz
    ) %>% unnest(c(accum,estimates))
pdatar

g <- pdatar %>%
    ggplot()+theme_classic()+
    aes(x=samples,y=ppis,group=idz)+
    geom_line()+geom_point()+
    geom_hline(aes(yintercept=boot),col="red",linetype="dashed")+
    scale_x_continuous(breaks=1:9)+
    scale_y_continuous(limits=c(0,NA),breaks=c(0,seq(1000,21e3,2e3)))+
    xlab("Number of environments sampled")+
    ylab("Number of unique PPIs,\nuncorrected for validation rate")
g

#ggsave(paste0(Sys.Date(),"_panel_ppis_accumulated_uncorrected.png"),g,
#    width=6,height=4)

#####
#####
#####
#####
#####

resampled_list_accumulation <- readRDS("thirtytwo_resampled_accumulations.RData")

pdatar <- 
    tibble(
        estimates=lapply(resampled_list_accumulation,
                function(list_of_accumulation) {
                    return(list_of_accumulation$species_pool)
                }
            )
        ,
        accum=lapply(resampled_list_accumulation,
                function(list_of_accumulation) {
                    with(list_of_accumulation$species_accumulation,
                        return(tibble(samples=sites,ppis=richness))
                    )
                }
            )
        ,
        idz=lapply(resampled_list_accumulation,
                function(list_of_accumulation) {
                    return(list_of_accumulation$idz)
                }
            ),
        summaryz=lapply(resampled_list_accumulation,
                function(list_of_accumulation) {
                    return(list_of_accumulation$summary_resampled_accumulation)
                }
            ),
        deduped=lapply(resampled_list_accumulation,
                function(list_of_accumulation) {
                    return(list_of_accumulation$dedup_resampled_accumulation)
                }
            )
    ) %>%
    unnest(cols=c(estimates,idz))
pdatar

g <- pdatar %>%
    {
    ggplot()+theme_classic()+
    geom_boxplot(data=select(.,deduped) %>% unnest(deduped),
        alpha=0.2,width=0.5,outlier.size=0.1,
        aes(x=as.numeric(samplez),
            group=as.numeric(samplez),
            y=ppis,weight=weight
        )
    ) +
    geom_line(data=select(.,accum,idz)%>%unnest(accum),
        alpha=0.2,col="red",aes(x=samples,y=ppis,group=idz))+
#    geom_point(data=select(.,accum,idz)%>%unnest(accum),
#        alpha=0.2,col="red",aes(x=samples,y=ppis,group=idz))+
    geom_hline(data=select(.,boot),#%>%summarize(boot=median(boot)),
        aes(yintercept=boot),col="red",linetype="dashed",alpha=0.1,size=0.5)+
    scale_x_continuous(breaks=1:9)+
    scale_y_continuous(limits=c(0,NA),breaks=c(0,seq(1000,21e3,2e3)))+
    xlab("Number of environments sampled")+
    ylab("Number of unique PPIs,\nresampled by predicted validation rate")
    }
g

ggsave_both(paste0(Sys.Date(),"_panel_ppi_accumulation_resampled_32_times"),
    g,width=6,height=4)


#####
#####
##### Doing just one resample with all of this
#####
#####

one_resampled_list_accumulation <- readRDS("one_resampled_accumulations.RData")

g <- one_resampled_list_accumulation %>%
    {
    ggplot(.$summary_resampled_accumulation)+theme_classic()+
    geom_boxplot(stat="identity",alpha=0.2,width=0.5,
        aes(
            x=as.numeric(samplez),
            group=as.numeric(samplez),
            ymin=`0%`,lower=`25%`,middle=`50%`,upper=`75%`,ymax=`100%`
        )
    )+
    geom_hline(aes(yintercept=.$species_pool$boot),col="red",linetype="dashed")+
    geom_line(data=with(
            one_resampled_list_accumulation$species_accumulation,
            tibble(samples=sites,ppis=richness)
        ),aes(x=samples,y=ppis),alpha=0.5,size=1.5,col="red"
    )+
    scale_x_continuous(breaks=1:9)+
    scale_y_continuous(limits=c(0,NA),breaks=c(0,seq(1000,21e3,2e3)))+
    xlab("Number of environments")+
    ylab("Number of unique PPIs,\nresampled by validation rates")
    }
g

#ggsave_both(paste0(Sys.Date(),"_accumulation_of_resampled_ppis_by_permuted"),
#    g,width=6,height=4)

