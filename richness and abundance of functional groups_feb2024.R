### Calculate and plot -- woody, forb, and grass richness and cover for upland lowland and slope plots, 
### Also organize and plot abundances across transect location (from shrub center) by functional group, also maybe split by topo
###
### Kevin Wilcox (k_wilcox@uncg.edu)
### created Feb 28, 2024

### Set up workspace
library(tidyverse)
library(ggthemes)
library(vegan)
library(ggpattern)
rm(list=ls())
setwd("C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\Konza projects_other\\ShrubRecolonization\\")
setwd("C:\\Users\\wilco\\OneDrive - UNCG\\Current projects\\Konza projects_other\\ShrubRecolonization\\") # personal laptop

SE_function<-function(x,na.rm=na.rm){
  SD = sd(x,na.rm=TRUE)
  denom = sqrt(length(x[!is.na(x)]))
  SE=SD/denom
  return(SE)
}

### Read in and prepare data 
### NOTE: there are zeros in some of the years of data which inflates richness and diversity calcs -- need to change these to NAs for analysis
### NOTE: 2018 data doesn't line up with later years. Also this is when they were just stems, so removing these data from all analyses
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
  filter_at(vars(abs_cover), all_vars(!is.na(.))) %>%
  filter(year!=2018) %>%
  mutate(growthform2 = ifelse(growthform=="p","p",
                              ifelse(growthform %in% c("a","b"),"a/b",NA)))


unique(spcomp_df$growthform)
unique(spcomp_df$lifeform)
unique(spcomp_df$fxn_group)

### calculate abundance and richness, then pivot to get 0s for times when there is no funcitonal abundance - need to do abundance and richness separately then merge them back together
abun_plot_df <- spcomp_df %>%
  group_by(year, transect, topo, plot, cover, shrubSP, fxn_group, growthform2) %>%
  summarize(
    abs_cover = sum(abs_cover)
  ) %>%
  pivot_wider(names_from=fxn_group, values_from = abs_cover) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols=forb:woody, names_to="fxn_group", values_to = "abs_cover") %>%
  ungroup()

rich_plot_df <- spcomp_df %>%
  group_by(year, transect, topo, plot, cover, shrubSP, fxn_group, growthform2) %>%
  summarize(
    richness = specnumber(genus_species)
  ) %>%
  pivot_wider(names_from=fxn_group, values_from = richness) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols=forb:woody, names_to="fxn_group", values_to = "richness")

# merge back together
abun_rich_plot_df <- abun_plot_df %>%
  full_join(rich_plot_df, by=c("year", "transect", "topo", "plot", "cover", "shrubSP", "fxn_group", "growthform2"))

### quick look at the raw plot data
ggplot(abun_rich_plot_df, aes(x=year, y=abs_cover, col=topo)) +
  geom_jitter(width=0.2) +
  facet_wrap(~fxn_group*cover) +
  theme_few()

ggplot(abun_rich_plot_df, aes(x=year, y=richness, col=topo)) +
  geom_jitter(width=0.2, alpha=0.5) +
  facet_wrap(~fxn_group*cover) +
  theme_few()


### calculate abundance and richness means and se's, for plotting
abun_rich_means <- abun_rich_plot_df %>%
  group_by(year, topo, cover, fxn_group, growthform2, transect) %>%
  summarize(abs_cover = mean(abs_cover, na.rm=T),
            richness = mean(richness, na.rm=T)) %>%
  ungroup() %>%
  group_by(year, topo, cover, fxn_group, growthform2) %>%
  summarize(
    abs_cover_u = mean(abs_cover, na.rm=T),
    abs_cover_se = SE_function(abs_cover),
    richness_u = mean(richness, na.rm=T),
    richness_se = SE_function(richness)
    ) %>%
  ungroup()


abun_rich_means$cover <- factor(abun_rich_means$cover, levels=c("G", "T", "Sh"))
### full abundance plot
ggplot(filter(abun_rich_means, cover!="T" & topo=="lowland" & cover=="G"), 
       aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group, alpha=growthform2)) +
#  geom_errorbar(width=.1, position = 'stack') +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~year) +
  theme_few() +
  scale_alpha_manual(values = c(0.4,1)) +
  scale_fill_manual(values=c("purple4", "springgreen4","burlywood4")) +
  ylim(0,130)

ggplot(filter(abun_rich_means, cover!="T" & topo=="lowland" & cover=="Sh"), 
       aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group, alpha=growthform2)) +
  #  geom_errorbar(width=.1, position = 'stack') +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~year) +
  theme_few() +
  scale_alpha_manual(values = c(0.4,1)) +
  ylim(0,130)

ggplot(filter(abun_rich_means, cover!="T" & topo=="upland" & cover=="G"), 
       aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group, alpha=growthform2)) +
  #  geom_errorbar(width=.1, position = 'stack') +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~year) +
  theme_few() +
  scale_alpha_manual(values = c(0.4,1)) +
  ylim(0,130)

ggplot(filter(abun_rich_means, cover!="T" & topo=="upland" & cover=="Sh"), 
       aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group, alpha=growthform2)) +
  #  geom_errorbar(width=.1, position = 'stack') +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~year) +
  theme_few() +
  scale_alpha_manual(values = c(0.4,1)) +
  ylim(0,130)

ggplot(filter(abun_rich_means, cover!="T" & topo=="slope" & cover=="G"), 
       aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group, alpha=growthform2)) +
  #  geom_errorbar(width=.1, position = 'stack') +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~year) +
  theme_few() +
  scale_alpha_manual(values = c(0.4,1)) +
  ylim(0,130)

ggplot(filter(abun_rich_means, cover!="T" & topo=="slope" & cover=="Sh"), 
       aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group, alpha=growthform2)) +
  #  geom_errorbar(width=.1, position = 'stack') +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~year) +
  theme_few() +
  scale_alpha_manual(values = c(0.4,1)) +
  ylim(0,130)






ggsave(paste0("figures//fxn group cover through time_full_bars", Sys.Date(),".png"), width=10, height=9, units="in")

# ggplot(filter(abun_rich_means, growthform2=="a/b"), aes(x=year, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, fill=fxn_group)) +
#   #  geom_errorbar(width=.1, position=position_dodge(width=0.9)) +
#   #  geom_path() +
#   #  geom_point(size=2, position=position_dodge(width=.2)) +
#   geom_col(position=position_dodge(), col="black") +
#   facet_wrap(~cover*topo) +
#   theme_few()


ggplot(abun_rich_means, aes(x=year, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, ymax=abs_cover_u+abs_cover_se, col=topo)) +
  geom_errorbar(width=.1) +
  geom_path() +
  geom_point(size=2) +
  facet_wrap(~cover*fxn_group) +
  theme_few()
ggsave(paste0("figures//fxn group cover through time_full_", Sys.Date(),".png"), width=10, height=9, units="in")

### full richness plot
ggplot(abun_rich_means, aes(x=year, y=richness_u, ymin=richness_u-richness_se, ymax=richness_u+richness_se, col=topo)) +
  geom_errorbar(width=.1) +
  geom_path() +
  geom_point(size=2) +
  facet_wrap(~cover*fxn_group) +
  theme_few()
ggsave(paste0("figures//fxn group richness through time_full_", Sys.Date(),".png"), width=10, height=9, units="in")



###
### Organize and plot abundances across transect location (from shrub center) by functional group, also split by topo
###



topo_df <- abun_rich_plot_df %>% ungroup() %>% select(topo, transect) %>% unique(.) 

topo_vec <- topo_key %>% pull(topo)
names(topo_vec) <- topo_key$transect # have to created a named vector for the labeller to work
                                                                                                
### Look at transect curves for each transect

ggplot(filter(abun_rich_plot_df, fxn_group=="woody" & year != 2018), aes(x=plot, y=abs_cover, label=cover, col=as.factor(year))) + 
  geom_path() + 
  geom_text() + 
  facet_wrap(~transect, labeller=labeller(transect = topo_vec)) +
  ylab("Shrub absolute cover (%)")
ggsave(paste0("figures//shrub cover along transect_year as colors", Sys.Date(),".png"), width=14, height=9, units="in")


### prep data frame
abun_dist_from_shrub <- abun_rich_plot_df %>%
  
?facet_wrap






