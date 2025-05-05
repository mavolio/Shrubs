### Calculate and plot -- woody, forb, and grass richness and cover for upland lowland and slope plots, 
### Also organize and plot abundances across transect location (from shrub center) by functional group, also maybe split by topo
###
### Kevin Wilcox (k_wilcox@uncg.edu)
### created Feb 28, 2024, last updated Nov 2024

### Set up workspace
library(tidyverse)
library(ggthemes)
library(vegan)
library(nlme)
library(lmerTest)
library(emmeans)
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
### NOTE: there are zeros in some of the years of data which inflates richness and diversity calcs -- need to change these to NAs for analysis
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
  summarize(abs_cover = mean(abs_cover, na.rm=T),
            richness = mean(richness, na.rm=T)) %>%
  ungroup() %>%
  group_by(topo, cover, fxn_group) %>%
  summarize(
    abs_cover_u = mean(abs_cover, na.rm=T),
    abs_cover_se = SE_function(abs_cover),
    richness_u = mean(richness, na.rm=T),
    richness_se = SE_function(richness)
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
  group_by(topo, cover, growthform2) %>%
  summarize(
    abs_cover_u = mean(abs_cover, na.rm=T),
    abs_cover_se = SE_function(abs_cover),
    rel_cover_u = mean(rel_cover, na.rm=T),
    rel_cover_se = SE_function(rel_cover)
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
  ylim(0,90) +
  ylab("Absolute cover (%)") + xlab("Functional Group")

ggsave("figures//cover facet by cover type and fxn group.png", width=4, height=7)
pdf("figures//cover facet by cover type and fxn group.pdf", width=4, height=7, useDingbats = F)
print(abun_shrub_barplot)
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
  scale_fill_manual(values=c("purple4", "springgreen4","burlywood4")) +
  facet_wrap(~cover) +
  theme_few() +
  ylab("Species Richness") + xlab("Topographic position")

ggsave("figures//richness by cover type fxn group and topo.png", width=6, height=3)
pdf("figures//richness by cover type fxn group and topo.pdf", width=6, height=3)
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


}


###
### Run models
###
{
  ###
  ### Full models
  abun_grassy_lme <- lme(abs_cover ~ factor(year)*topo*fxn_group
                      , data=filter(abun_rich_transect_df, cover=="G")
                      , random = ~1 |transect
#                      , correlation=corCompSymm(form = ~1 |transect)
                      , correlation=corAR1(form = ~1 |transect)
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
  
  anova.lme(abun_grassy_lme, type="marginal")
#  anova(abun_grassy_lme, type="2")
  summary(abun_grassy_lme)
  AIC(abun_grassy_lme)
  emmeans(abun_grassy_lme, pairwise ~ fxn_group*year, by="topo", adjust="sidak")
  
  # abun_grassy_lm <- lm(abs_cover ~ year*topo*fxn_group,
  #                     , data=filter(abun_rich_transect_df, cover=="G")
  #                     , na.action = na.omit)
  # 
  # anova(abun_grassy_lme, type="marginal")
  # summary(abun_grassy_lme)
  # Anova(abun_grassy_lm, type=3)

  abun_shrub_lme <- lme(abs_cover ~ factor(year)*topo*fxn_group,
                         , data=filter(abun_rich_transect_df, cover=="Sh")
                         , random = ~1 |transect
  #                       , correlation=corCompSymm(form = ~1 |transect)
                        , correlation=corAR1(form = ~1 |transect)
                        , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)

  anova(abun_shrub_lme, type="marginal")
  emmeans(abun_shrub_lme, pairwise ~ fxn_group, adjust="sidak")
  
  
  richness_grassy_lme_1 <- lme(richness ~ factor(year)*topo*fxn_group,
                         , data=filter(abun_rich_transect_df, cover=="G")
                         , random = ~1 |transect
                         , correlation=corCompSymm(form = ~1 |transect)
                         , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)
  AIC(richness_grassy_lme_1)
  anova(richness_grassy_lme_1, type="marginal") # year*fxn_group is highest P value so remove
  richness_grassy_lme_2 <- lme(richness ~ factor(year)*topo+factor(year)*fxn_group+topo*fxn_group,
                               , data=filter(abun_rich_transect_df, cover=="G")
                               , random = ~1 |transect
                               , correlation=corCompSymm(form = ~1 |transect)
                               , control=lmeControl(returnObject=TRUE)
                               , na.action = na.omit)
  AIC(richness_grassy_lme_2) # 187 v 188 -- Not better, use full model
  anova(richness_grassy_lme_2, type="marginal") # year*fxn_group is highest P value so remove

  emmeans(richness_grassy_lme_1, pairwise ~ fxn_group, by="topo", adjust="sidak")
  emmeans(richness_grassy_lme_1, pairwise ~ fxn_group*topo, adjust="sidak")
  
  
  
  
  
   richness_shrub_lme_1 <- lme(richness ~ factor(year)*topo*fxn_group,
                         , data=filter(abun_rich_transect_df, cover=="Sh")
                         , random = ~1 |transect
                         , correlation=corCompSymm(form = ~1 |transect)
                         , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)
  
   AIC(richness_shrub_lme_1) # Remove 3 way interaction
   anova(richness_shrub_lme_1, type="marginal") 
   richness_shrub_lme_2 <- lme(richness ~ factor(year)*topo+factor(year)*fxn_group+topo*fxn_group,
                               , data=filter(abun_rich_transect_df, cover=="Sh")
                               , random = ~1 |transect
                               , correlation=corCompSymm(form = ~1 |transect)
                               , control=lmeControl(returnObject=TRUE)
                               , na.action = na.omit)
   AIC(richness_shrub_lme_2) # 237 vs 246 - full model has lower AIC, stop and use full model
  
   emmeans(richness_shrub_lme_1, pairwise ~ fxn_group, by="topo", adjust="sidak")
   emmeans(richness_shrub_lme_1, pairwise ~ fxn_group*topo, adjust="sidak")
   

  ###
  ### Split models by topo

  # Abundance in grassy areas
  # upland
  abun_grassy_upland_lme <- lme(abs_cover ~ factor(year)*fxn_group,
                                , data=filter(abun_rich_transect_df, cover=="G" & topo=="upland")
                                , random = ~1 |transect
                                #                      , correlation=corCompSymm(form = ~1 |transect)
                                , correlation=corAR1(form = ~1 |transect)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
  
  anova.lme(abun_grassy_upland_lme, type="marginal")

  emmeans(abun_grassy_upland_lme, pairwise ~ fxn_group, by="year", adjust="sidak")

  # slope
  abun_grassy_slope_lme <- lme(abs_cover ~ factor(year)*fxn_group,
                                , data=filter(abun_rich_transect_df, cover=="G" & topo=="slope")
                                , random = ~1 |transect
                                #                      , correlation=corCompSymm(form = ~1 |transect)
                                , correlation=corAR1(form = ~1 |transect)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
  
  anova.lme(abun_grassy_slope_lme, type="marginal")
  
  emmeans(abun_grassy_slope_lme, pairwise ~ fxn_group, by="year", adjust="sidak")
  
  # lowland
  abun_grassy_lowland_lme <- lme(abs_cover ~ factor(year)*fxn_group,
                               , data=filter(abun_rich_transect_df, cover=="G" & topo=="lowland")
                               , random = ~1 |transect
                               #                      , correlation=corCompSymm(form = ~1 |transect)
                               , correlation=corAR1(form = ~1 |transect)
                               , control=lmeControl(returnObject=TRUE)
                               , na.action = na.omit)
  
  anova.lme(abun_grassy_lowland_lme, type="marginal")
  
  emmeans(abun_grassy_lowland_lme, pairwise ~ fxn_group, by="year", adjust="sidak")
  
  
  # Abundance in shrubby areas
  # upland
  abun_shrubby_upland_lme <- lme(abs_cover ~ factor(year)*fxn_group,
                                , data=filter(abun_rich_transect_df, cover=="Sh" & topo=="upland")
                                , random = ~1 |transect
                                #                      , correlation=corCompSymm(form = ~1 |transect)
                                , correlation=corAR1(form = ~1 |transect)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
  
  anova.lme(abun_shrubby_upland_lme, type="marginal")
  
  emmeans(abun_shrubby_upland_lme, pairwise ~ fxn_group,  adjust="sidak")
  
  # slope
  abun_shrubby_slope_lme <- lme(abs_cover ~ factor(year)*fxn_group,
                               , data=filter(abun_rich_transect_df, cover=="Sh" & topo=="slope")
                               , random = ~1 |transect
                               #                      , correlation=corCompSymm(form = ~1 |transect)
                               , correlation=corAR1(form = ~1 |transect)
                               , control=lmeControl(returnObject=TRUE)
                               , na.action = na.omit)
  
  anova.lme(abun_shrubby_slope_lme, type="marginal")
  
  emmeans(abun_shrubby_slope_lme, pairwise ~ fxn_group, adjust="sidak")
  
  # lowland
  abun_shrubby_lowland_lme <- lme(abs_cover ~ factor(year)*fxn_group,
                                 , data=filter(abun_rich_transect_df, cover=="Sh" & topo=="lowland")
                                 , random = ~1 |transect
                                 #                      , correlation=corCompSymm(form = ~1 |transect)
                                 , correlation=corAR1(form = ~1 |transect)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit)
  
  anova.lme(abun_shrubby_lowland_lme, type="marginal")
  
  emmeans(abun_shrubby_lowland_lme, pairwise ~ fxn_group, adjust="sidak")
  
  ###
  ### Annual relative cover models
  
  # overall model - MIGHT WANT TO TRANSFORM THIS, SQRT(Y) LOOKED OK
  annual_relcov_lme <- lme(sqrt(rel_cover) ~ factor(year)*topo*cover,
                           , data=filter(annual_relcov_transect, growthform2=="a/b" & cover %in% c("G","Sh"))
                           , random = ~1 |transect
                           #                      , correlation=corCompSymm(form = ~1 |transect)
                           , correlation=corAR1(form = ~1 |transect)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit)
  
  anova.lme(annual_relcov_lme, type="marginal")
  check_model(annual_relcov_lme)
  plot(annual_relcov_lme)
  
  emmeans(annual_relcov_lme, pairwise ~ cover, by="topo",  adjust="sidak")
  emmeans(annual_relcov_lme, pairwise ~ cover,  adjust="sidak")
  
  
  # upland
  annual_relcov_upland_lme <- lme(sqrt(rel_cover) ~ factor(year)*cover,
                                 , data=filter(annual_relcov_transect, growthform2=="a/b" & topo=="upland")
                                 , random = ~1 |transect
                                 #                      , correlation=corCompSymm(form = ~1 |transect)
                                 , correlation=corAR1(form = ~1 |transect)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit)
  
  check_model(annual_relcov_upland_lme)
  anova.lme(annual_relcov_upland_lme, type="marginal")
  
  emmeans(annual_relcov_upland_lme, pairwise ~ fxn_group,  adjust="sidak")
  
  # slope
  annual_relcov_slope_lme <- lme(sqrt(rel_cover) ~ factor(year)*cover,
                                  , data=filter(annual_relcov_transect, growthform2=="a/b" & topo=="slope")
                                  , random = ~1 |transect
                                  #                      , correlation=corCompSymm(form = ~1 |transect)
                                  , correlation=corAR1(form = ~1 |transect)
                                  , control=lmeControl(returnObject=TRUE)
                                  , na.action = na.omit)
  
  check_model(annual_relcov_slope_lme)
  anova.lme(annual_relcov_slope_lme, type="marginal")
  
  emmeans(annual_relcov_slope_lme, pairwise ~ fxn_group,  adjust="sidak")
  
  # lowland
  annual_relcov_lowland_lme <- lme(sqrt(rel_cover) ~ factor(year)*cover,
                                  , data=filter(annual_relcov_transect, growthform2=="a/b" & topo=="lowland" &cover!="T")
                                  , random = ~1 |transect
                                  #                      , correlation=corCompSymm(form = ~1 |transect)
                                  , correlation=corAR1(form = ~1 |transect)
                                  , control=lmeControl(returnObject=TRUE)
                                  , na.action = na.omit)
  
  anova.lme(annual_relcov_lowland_lme, type="marginal")
  
  emmeans(annual_relcov_lowland_lme, pairwise ~ cover, by="year",  adjust="sidak")
  check_model(annual_relcov_lowland_lme)
}














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






