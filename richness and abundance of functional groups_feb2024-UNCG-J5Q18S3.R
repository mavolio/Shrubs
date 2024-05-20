### Calculate and plot -- woody, forb, and grass richness and cover for upland lowland and slope plots
###
### Kevin Wilcox (k_wilcox@uncg.edu)
### created Feb 28, 2024

### SEt up workspace
library(tidyverse)
library(ggthemes)
library(vegan)

rm(list=ls())
setwd("C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\Konza projects_other\\ShrubRecolonization\\")

SE_function<-function(x,na.rm=na.rm){
  SD = sd(x,na.rm=TRUE)
  denom = sqrt(length(x[!is.na(x)]))
  SE=SD/denom
  return(SE)
}

### Read in and prepare data 
### NOTE: there are zeros in some of the years of data which inflates richness and diversity calcs -- need to change these to NAs for analysis
splist <- read.csv("KNZ_species_list.csv") %>%
  rename(sp_code=code)

topo_key <- read.csv("topo_key.csv")

spcomp_df <- read.csv("species_compositon_allyears_Feb2024.csv") %>%
  left_join(splist, by="sp_code") %>%
  mutate(fxn_group = ifelse(lifeform %in% c("f", "m", "o"), "forb",
                            ifelse(lifeform %in% c("g", "s"), "graminoid",
                                   ifelse(lifeform=="w", "woody", NA)))
  ) %>%
  left_join(topo_key, by="transect") %>%
  mutate(abs_cover = replace(abs_cover, abs_cover==0, NA)) %>%
  filter_at(vars(abs_cover), all_vars(!is.na(.)))


unique(spcomp_df$growthform)
unique(spcomp_df$lifeform)
unique(spcomp_df$fxn_group)

abun_rich_plot_df <- spcomp_df %>%
  group_by(year, transect, topo, plot, cover, shrubSP, fxn_group) %>%
  summarize(
    abs_cover = sum(abs_cover),
    richness = specnumber(genus_species)
  )

abun_rich_plot_df %>%
  filter(year==2019&transect==1 & plot==1)

### quick look at the raw plot data
ggplot(abun_rich_plot_df, aes(x=year, y=abs_cover, col=topo)) +
  geom_jitter(width=0.2) +
  facet_wrap(~fxn_group*cover) +
  theme_few()

ggplot(abun_rich_plot_df, aes(x=year, y=richness, col=topo)) +
  geom_jitter(width=0.2, alpha=0.5) +
  facet_wrap(~fxn_group*cover) +
  theme_few()

abun_rich_df %>% filter(richness>40)


