### Calculate Analyze and plot -- woody, forb, and grass richness and cover for upland lowland and slope plots, 
### Also organize and plot abundances across transect location (from shrub center) by functional group, also maybe split by topo
###
### Kevin Wilcox (k_wilcox@uncg.edu)
### created Feb 28, 2024, last updated May 5 2025

### Set up workspace
library(tidyverse)
library(ggthemes)
library(vegan)
library(nlme)
library(lmerTest)
library(emmeans)
library(performance)
library(car)
#library(ggpattern)
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
### NOTE: there are zeros in some of the years of data which inflates richness and diversity calcs 
###   -- SPCOMP_DF NO LONGER HAS ANY ZEROS, SO THIS ISSUE SHOULD BE FIXED
### NOTE: 2018 data doesn't line up with later years. Also this is when they were just stems, so removing these data from all analyses
{
splist <- read.csv("KNZ_species_list.csv") %>%
  rename(sp_code=code)

topo_key <- read.csv("topo_key.csv")
plot_key <- read.csv("shrub transect locations.csv")

spcomp_df <- read.csv("species_compositon_allyears_Feb2024.csv") %>%
  dplyr::select(-X, -cover, -X2022.data) %>%
  full_join(plot_key, by=c("transect", "plot")) %>%
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

### AVerage across plots within transects for running models -- WE SHOULD DOUBLE CHECK THAT WE HAVE VALUES FOR EVERY CELL AND THAT THERE AREN'T ANY MISSING ROWS
abun_rich_transect_df <- abun_rich_plot_df %>%
  group_by(year, transect, topo, cover, fxn_group) %>%
  summarize(abs_cover=mean(abs_cover, na.rm=T),
            richness= mean(richness, na.rm=T)) %>%
  ungroup()
with(abun_rich_transect_df, table(year, transect))

}


###
### Calculate abundance and richness means and se's, for plotting -- 
###
{
abun_rich_means <- abun_rich_plot_df %>%
  group_by(year, topo, cover, fxn_group, transect) %>%
  summarize(abs_cover = mean(abs_cover, na.rm=T),
            richness = mean(richness, na.rm=T)) %>%
  ungroup() %>%
  group_by(year, topo, cover, fxn_group) %>%
  summarize(
    abs_cover_u = mean(abs_cover, na.rm=T),
    abs_cover_se = SE_function(abs_cover),
    richness_u = mean(richness, na.rm=T),
    richness_se = SE_function(richness)
    ) %>%
  ungroup()
  
abun_rich_means_topo_fxn <- abun_rich_plot_df %>%
  group_by(year, topo, cover, fxn_group, transect) %>%
  summarize_at(vars(abs_cover,richness),.funs = c(mean),na.rm=TRUE) %>%
  ungroup() %>%
  group_by(year, topo, cover, fxn_group) %>%
  summarize_at(vars(abs_cover,richness),.funs = c(mean=mean,se=SE_function),na.rm=TRUE) %>%
  ungroup() %>%
  group_by(topo, cover, fxn_group) %>%
  summarize(
    abs_cover_u = mean(abs_cover_mean, na.rm=T),
    abs_cover_se = mean(abs_cover_se, na.rm=T), # First calculate SE by year then average SE across years to get SE (this way it doesn't inflate your N)
    richness_u = mean(richness_mean, na.rm=T),
    richness_se = mean(richness_se, na.rm=T)
  ) %>%
  ungroup()

annual_relcov_transect <- abun_rich_plot_df %>%
  group_by(year, topo, cover, growthform2, transect) %>%
  summarize(abs_cover = mean(abs_cover, na.rm=T),
            richness = mean(richness, na.rm=T)) %>%
  ungroup() %>%
  group_by(year, topo, cover, transect) %>%
  mutate(rel_cover = abs_cover/sum(abs_cover,na.rm=T)) %>% # I now calculate relcover here for the model, and this might muck up things downstream -- just screwing future Kevin over yet again!
  ungroup() 
  

annual_means <- annual_relcov_transect %>%
  group_by(year, topo, cover, growthform2) %>%
  summarize(
    abs_cover_u = mean(abs_cover, na.rm=T),
    abs_cover_se = SE_function(abs_cover),
    richness_u = mean(richness, na.rm=T),
    richness_se = SE_function(richness)
  ) %>%
  ungroup()

annual_relcov <- abun_rich_plot_df %>%
  group_by(year, topo, cover, growthform2, transect) %>%
  summarize(abs_cover = mean(abs_cover, na.rm=T),
            richness = mean(richness, na.rm=T)) %>%
  ungroup() %>%
  dplyr::select(-richness) %>%
  pivot_wider(names_from=growthform2, values_from = abs_cover) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols=c("a/b":"p"), names_to = "growthform2", values_to = "abs_cover")

annual_relcov_means <- annual_relcov %>%
  group_by(year, topo, cover, transect) %>%
  mutate(rel_cover = abs_cover/sum(abs_cover,na.rm=T)) %>%
  ungroup() %>%
  group_by(year, topo, cover, growthform2) %>%
  summarize(
    abs_cover_u = mean(abs_cover, na.rm=T),
    abs_cover_se = SE_function(abs_cover),
    rel_cover_u = mean(rel_cover, na.rm=T),
    rel_cover_se = SE_function(rel_cover)
  ) %>%
  ungroup()

annual_relcov_noYr_means <- annual_relcov %>%
  group_by(year, topo, cover, transect) %>%
  mutate(rel_cover = abs_cover/sum(abs_cover,na.rm=T)) %>%
  ungroup() %>%
  group_by(topo, cover, growthform2, year) %>%
  summarize_at(vars(abs_cover,rel_cover),.funs = c(mean=mean,se=SE_function),na.rm=TRUE) %>%
  ungroup() %>%
  group_by(topo, cover, growthform2) %>%
  summarize(
    abs_cover_u = mean(abs_cover_mean, na.rm=T),
    abs_cover_se = mean(abs_cover_se, na.rm=T), # First calculate SE by year then average SE across years to get SE (this way it doesn't inflate your N)
    rel_cover_u = mean(rel_cover_mean, na.rm=T),
    rel_cover_se = mean(rel_cover_se, na.rm=T)
  ) %>%
  ungroup()

annual_relcov_noYr_means$topo <- factor(annual_relcov_noYr_means$topo, levels=c("upland","slope","lowland"))
}


###
### Plotting
###
{
  abun_rich_means_simple <- abun_rich_means %>% filter(cover!="T")
  abun_rich_means_topo_fxn_simple <- abun_rich_means_topo_fxn %>% filter(cover!="T")
  
  
abun_rich_means_simple$topo <- factor(abun_rich_means_simple$topo, levels=c("upland", "slope", "lowland"))
abun_rich_means_topo_fxn_simple$topo <- factor(abun_rich_means_topo_fxn_simple$topo, levels=c("upland", "slope", "lowland"))
annual_relcov_means$topo <- factor(annual_relcov_means$topo, levels=c("upland", "slope", "lowland"))

### Abundance
abun_plot <- ggplot(filter(abun_rich_means_simple, cover=="G"), aes(x=year, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, 
                                               ymax=abs_cover_u+abs_cover_se, col=fxn_group)) +
  geom_path() + 
  geom_point(size=2) +
  geom_errorbar(width=0.1) +
  scale_color_manual(values=c("purple4", "springgreen4","burlywood4")) +
  facet_grid(topo ~ .) +
  theme_few() +
  ylim(0,90)+
  ylab("Absolute cover (%)") + xlab("Year")
  
ggsave("figures//cover through time facet by cover type and fxn group_grassy.png", width=4, height=7)
pdf("figures//cover through time facet by cover type and fxn group_grassy.pdf", width=4, height=7, useDingbats = F)
print(abun_plot)
dev.off()

# Abundance grassy patches bar plots not through time
abun_rich_means_topo_fxn$fxn_group <- factor(abun_rich_means_topo_fxn$fxn_group, levels=c("graminoid","woody","forb"))
abun_rich_means_topo_fxn$topo <- factor(abun_rich_means_topo_fxn$topo, levels=c("upland","slope","lowland"))
abun_rich_means_topo_fxn_simple$fxn_group <- factor(abun_rich_means_topo_fxn_simple$fxn_group, levels=c("graminoid","woody","forb"))
abun_grassy_barplot <- ggplot(filter(abun_rich_means_topo_fxn, cover=="G"),
                             aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, 
                                 ymax=abs_cover_u+abs_cover_se, fill=fxn_group)) +
  geom_col(col="black") +
  geom_errorbar(width=0.1, col="black") +
  scale_fill_manual(values=c("springgreen4","burlywood4", "purple4")) +
  facet_grid(topo~.) +
  theme_few() +
  ylim(0,90) +
  ylab("Absolute cover (%)") + xlab("Functional Group")

ggsave("figures//grassy - cover facet by cover type and fxn group.png", width=4, height=7)
pdf("figures//grassy - cover facet by cover type and fxn group.pdf", width=4, height=7, useDingbats = F)
print(abun_grassy_barplot)
dev.off()

# Abundance shrub islands bar plots not through time
abun_rich_means_topo_fxn$fxn_group <- factor(abun_rich_means_topo_fxn$fxn_group, levels=c("woody","graminoid","forb"))
abun_rich_means_topo_fxn$topo <- factor(abun_rich_means_topo_fxn$topo, levels=c("upland","slope","lowland"))
abun_shrub_barplot <- ggplot(filter(abun_rich_means_topo_fxn, cover=="Sh"),
                                    aes(x=fxn_group, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, 
                                       ymax=abs_cover_u+abs_cover_se, fill=fxn_group)) +
  geom_col(col="black") +
  geom_errorbar(width=0.1, col="black") +
  scale_fill_manual(values=c("burlywood4", "springgreen4","purple4")) +
  facet_grid(topo~.) +
  theme_few() +
#  ylim(0,90) +
  ylab("Absolute cover (%)") + xlab("Functional Group")

ggsave("figures//cover facet by cover type and fxn group.png", width=4, height=7)
pdf("figures//cover facet by cover type and fxn group.pdf", width=4, height=7, useDingbats = F)
print(abun_shrub_barplot)
dev.off()

### Abundance not through time all together
abun_noTime_plot <- ggplot(abun_rich_means_topo_fxn_simple, aes(x=topo, y=abs_cover_u, ymin=abs_cover_u-abs_cover_se, 
                                            ymax=abs_cover_u+abs_cover_se, fill=fxn_group)) +
  geom_col(position=position_dodge(width=.8), width=.7) +
  geom_errorbar(width=0.1, position=position_dodge(width=.8), col="black") +
  scale_fill_manual(values=c("springgreen4","burlywood4","purple4")) +
  facet_wrap(~cover) +
  theme_few() +
  ylab("Absolute cover (%)") + xlab("Topographic position")

ggsave("figures//abun by cover type fxn group and topo.png", width=6.5, height=3)
pdf("figures//abun by cover type fxn group and topo_v2.pdf", width=6.5, height=3)
print(abun_noTime_plot)
dev.off()


### Richness
## Full faceting
ggplot(abun_rich_means_simple, aes(x=year, y=richness_u, ymin=richness_u-richness_se, 
                                   ymax=richness_u+richness_se, col=fxn_group)) +
  geom_path() + 
  geom_point(size=2) +
  geom_errorbar(width=0.1) +
  scale_color_manual(values=c("purple4", "springgreen4","burlywood4")) +
  facet_grid(topo ~ cover) +
  theme_few() +
  ylab("Species Richness") + xlab("Year")

ggsave("figures//richness through time facet by cover type and fxn group.png", width=6, height=7)

## Richness -- just topo and functional group
rich_plot <- ggplot(abun_rich_means_topo_fxn_simple, aes(x=topo, y=richness_u, ymin=richness_u-richness_se, 
                                   ymax=richness_u+richness_se, fill=fxn_group)) +
  geom_col(position=position_dodge(width=.8), width=.7) +
  geom_errorbar(width=0.1, position=position_dodge(width=.8), col="black") +
  scale_fill_manual(values=c("springgreen4","burlywood4","purple4")) +
  facet_wrap(~cover) +
  theme_few() +
  ylab("Species Richness") + xlab("Topographic position")

ggsave("figures//richness by cover type fxn group and topo.png", width=6.5, height=3)
pdf("figures//richness by cover type fxn group and topo.pdf", width=6.5, height=3)
print(rich_plot)
dev.off()

### Annual cover
## Split by topo cover and year
annual_cov_plot <- ggplot(filter(annual_relcov_means,cover!="T" & growthform2=="a/b"), aes(x=year, y=rel_cover_u, ymin=rel_cover_u-rel_cover_se, 
                                   ymax=rel_cover_u+rel_cover_se, pch=cover, lty=cover)) +
  geom_path() + 
  geom_point(size=2) +
  geom_errorbar(width=0.1, lty=1) +
#  scale_color_manual(values=c("purple4", "springgreen4","burlywood4")) +
  facet_grid(topo ~ .) +
  scale_linetype_manual(values=c(2,1)) +
  theme_few() +
  ylab("Relative Cover of Annuals (%)") + xlab("Year")

ggsave("figures//annual relcov through time facet by cover type and fxn group.png", width=4, height=7)
pdf("figures//annual relcov through time facet by cover type and fxn group.pdf", width=4, height=7, useDingbats = F)
print(annual_cov_plot)
dev.off()

## Split by topo and cover only
annual_cov_barplot <- ggplot(filter(annual_relcov_noYr_means,cover!="T" & growthform2=="a/b"), aes(x=cover, y=rel_cover_u*100, ymin=rel_cover_u*100-rel_cover_se*100, 
                                                                                           ymax=rel_cover_u*100+rel_cover_se*100, fill=cover)) +
  geom_col(col="black") + 
#  geom_point(size=2) +
  geom_errorbar(width=0.1, lty=1) +
  #  scale_color_manual(values=c("purple4", "springgreen4","burlywood4")) +
  facet_grid(topo ~ .) +
  scale_fill_manual(values=c("white","darkgrey")) +
  theme_few() +
  ylab("Relative Cover of Annuals (%)") + xlab("")

ggsave("figures//annual relcov facet by cover type and fxn group.png", width=3.5, height=7)
pdf("figures//annual relcov facet by cover type and fxn group.pdf", width=3.5, height=7, useDingbats = F)
print(annual_cov_barplot)
dev.off()
dev.new()
?dev.set()
### Annual cover split by topo and cover into one panel
annual_noTime_plot <- ggplot(filter(annual_relcov_noYr_means,cover!="T"&growthform2=="a/b"), aes(x=topo, y=rel_cover_u, ymin=rel_cover_u-rel_cover_se, 
                                                                ymax=rel_cover_u+rel_cover_se, fill=cover)) +
  geom_col(position=position_dodge(width=.8), width=.7, col="black") +
  geom_errorbar(width=0.1, position=position_dodge(width=.8), col="black") +
  scale_fill_manual(values=c("white","darkgrey")) +
  theme_few() +
  ylab("Relative cover of annuals (%)") + xlab("Topographic position")

ggsave("figures//annual relcov by topo_one panel.png", width=4, height=3)
pdf("figures//annual relcov by topo_one panel.pdf", width=4, height=3)
print(annual_noTime_plot)
dev.off()

}


###
### Run models
###
{
  ###
  ### Full models -- repeated measures
  ###
  
  ### Abundance in grass dominated areas
  abun_grassy_lme <- lme(abs_cover ~ factor(year)*topo*fxn_group
                      , data=filter(abun_rich_transect_df, cover=="G")
                      , random = ~1 |transect
                      , correlation=corAR1(form = ~1 |transect)
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
  
  anova.lme(abun_grassy_lme, type="marginal")
  AIC(abun_grassy_lme)
  qqnorm(abun_grassy_lme, abline = c(0,1)) ## qqplot
  check_model(abun_grassy_lme) ## residuals and normality of resids
  emmeans(abun_grassy_lme, pairwise ~ fxn_group*year, by="topo", adjust="sidak")
  emmeans(abun_grassy_lme, pairwise ~ fxn_group, by=c("topo","year"), adjust="sidak")
  
  ## Split by year to run emmeans 
  
  
  
  ### Abundance in shrub islands
  abun_shrub_lme <- lme(abs_cover ~ factor(year)*topo*fxn_group,
                         , data=filter(abun_rich_transect_df, cover=="Sh")
                         , random = ~1 |transect
                        , correlation=corAR1(form = ~1 |transect)
                        , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)

  anova(abun_shrub_lme, type="marginal")
  qqnorm(abun_shrub_lme, abline = c(0,1)) ## qqplot
  check_model(abun_shrub_lme) ## residuals and normality of resids
  
  emmeans(abun_shrub_lme, pairwise ~ fxn_group, adjust="sidak")
  emmeans(abun_shrub_lme, pairwise ~ fxn_group, by="topo", adjust="sidak")
  
  ### Richness model in grass dominated areas
  richness_grassy_lme <- lme(richness ~ factor(year)*topo*fxn_group,
                         , data=filter(abun_rich_transect_df, cover=="G")
                         , random = ~1 |transect
                         , correlation=corCompSymm(form = ~1 |transect)
                         , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)
  AIC(richness_grassy_lme)
  anova(richness_grassy_lme, type="marginal")
  qqnorm(richness_grassy_lme, abline = c(0,1)) ## qqplot
  check_model(richness_grassy_lme) ## residuals and normality of resids
  
  emmeans(richness_grassy_lme_1, pairwise ~ fxn_group, by="topo", adjust="sidak")
  emmeans(richness_grassy_lme_1, pairwise ~ fxn_group*topo, adjust="sidak")
  
  ### Richness models in shrub islands
   richness_shrub_lme <- lme(richness ~ factor(year)*topo*fxn_group,
                         , data=filter(abun_rich_transect_df, cover=="Sh")
                         , random = ~1 |transect
                         , correlation=corCompSymm(form = ~1 |transect)
                         , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)
  
   AIC(richness_shrub_lme) 
   anova(richness_shrub_lme, type="marginal") 
   qqnorm(richness_shrub_lme, abline = c(0,1)) ## qqplot
   check_model(richness_shrub_lme) ## residuals and normality of resids
   
   emmeans(richness_shrub_lme, pairwise ~ fxn_group, by="topo", adjust="sidak")
   emmeans(richness_shrub_lme, pairwise ~ fxn_group*topo, adjust="sidak")
   

  ###
  ### Annual relative cover models
  
  # overall model -  SQRT(Y) FOR MODEL ASSUMPTIONS
  annual_relcov_lme <- lme(sqrt(rel_cover) ~ factor(year)*topo*cover,
                           , data=filter(annual_relcov_transect, growthform2=="a/b" & cover %in% c("G","Sh"))
                           , random = ~1 |transect
                           , correlation=corAR1(form = ~1 |transect)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit)
  
  anova.lme(annual_relcov_lme, type="marginal")
  qqnorm(annual_relcov_lme, abline = c(0,1)) ## qqplot
  check_model(annual_relcov_lme) ## residuals and normality of resids
  
  emmeans(annual_relcov_lme, pairwise ~ cover, by="topo",  adjust="sidak")
  emmeans(annual_relcov_lme, pairwise ~ cover,  adjust="sidak")
  
  
}

