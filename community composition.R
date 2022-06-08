### Calculate and visualize shrub trajectories post-burn
###
### Author: Meghan Avolio, Kevin Wilcox (wilcoxkr@gmail.com)
### Last updated June 7, 2022


### Set up workstation
setwd("C://Users/mavolio2/Dropbox/Konza Research/Shrub Islands/")
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\Konza projects_other\\ShrubRecolonization\\")

library(tidyverse)
library(codyn)
library(vegan)

### Standard Error function
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

###
### 1. for each transect, designate a plot as grassland, transitional, or burned, based on 2018 data
###
{
post <- read.csv("Woody removal plots_K20A_after.csv") %>% 
  rename(transect=Transect,
         plot=Plot..Frame.)

bareground<-post%>%
  filter(Species=="Other- Bare Ground") # there is so much bareground I am not sure this is going to work

hist(bareground$Abundance)

pre<-read.csv("Species_Comp_k20a_before.csv")

# splist<-pre %>% 
#   select(Genus_Species) %>% 
#   unique()
# 
# write.csv(splist, "Data 2018/pre_species.csv", row.names = F)
lf<-read.csv("pre_species_lifeform.csv")

plotslist<-pre%>%
  select(Trasect, Frame) %>% 
  rename(Transect=Trasect) %>%
  unique()

# pre_class_total<-pre %>% 
#   left_join(lf) %>% 
#   filter(Lifeform!="Drop")%>%
#   group_by(Trasect, Frame) %>% 
#   summarize(totcov=sum(Cover_class))
#   
# pre_class<-pre %>% 
#   left_join(lf) %>% 
#   filter(Lifeform!="Drop")%>%
#   group_by(Trasect, Frame, Lifeform) %>% 
#   summarize(cov=sum(Cover_class)) %>% 
#   left_join(pre_class_total) %>% 
#   mutate(relcov=cov/totcov) %>% 
#   mutate(class=ifelse(Lifeform=="Woody"&relcov>0.55, 1, 0))%>%
#   filter(class==1)%>%
#   select(Trasect, Frame, class) %>% 
#   unique()

cover_class_key <- data.frame(
  Cover_class = 1:7,
  abs_cover = c(0.5, 3, 15, 37.5, 62.5, 85, 97.5)
)

pre_class<-pre %>%  ## I changed this to actual covers because the cover classes are not linear so calculating relcov with them doesn't work
  left_join(lf) %>% 
  filter(Lifeform!="Drop")%>%
  rename(Transect=Trasect) %>%
  group_by(Transect, Frame, Lifeform) %>% 
  left_join(cover_class_key, by="Cover_class") %>%
  summarize(abs_cover=sum(abs_cover)) %>% 
  ungroup() %>%
  group_by(Transect, Frame) %>%
  mutate(relcov=abs_cover/sum(abs_cover)) %>% 
  mutate(class=ifelse(Lifeform=="Woody"&abs_cover>40, 1, 0))%>%
  filter(class==1)%>%
  select(Transect, Frame, class) %>% 
  unique()

ggplot(data=subset(pre_class, Lifeform=="Woody"), aes(x=Frame, y=abs_cover))+
  geom_point()+
  geom_line()+
  facet_wrap(~Transect)+
  geom_hline(yintercept = 45)

testing <- pre_class %>%
  spread(key=Lifeform, value=abs_cover) %>%
  replace(is.na(.),0)
  
  
# Data frame identifying pre-burn grass versus woody states
plotcat<-plotslist%>%
  left_join(pre_class) %>% 
  mutate(category=ifelse(is.na(class), "grassy", 'woody')) %>% 
  rename(transect=Transect,
         plot=Frame) %>% 
  select(-class)

#how do classes differ in bareground?
bg<-bareground%>%
  left_join(plotcat)

boxplot(bg$Abundance~bg$category)  
summary(aov(Abundance~category, data = bg))#no sig diff

}

###
### 2. Clean and combine species comp data from 2018, 2019, 2021
###
{

species_key <- read.csv("KNZ_species_list.csv") %>%##meghan changed this slightly.
    mutate(Genus_Species = paste(genus, species, sep="_")) %>%
    bind_rows(data.frame(code=551, gen="quincu", spec="lobata", genus="quincula", species="lobata", family="solanaceae",
                         growthform="p", lifeform="f",origin="n", photo="unkn", Genus_Species="quincula_lobata")) %>%
    bind_rows(data.frame(code=552, gen="penste", spec="spp", genus="penstemon", species="spp", family="acanthaceae", # we should probably figure this one out and fix this
                         growthform="unkn", lifeform="f",origin="unkn", photo="unkn", Genus_Species="penstemon_spp")) %>%
    bind_rows(data.frame(code=553, gen="solida", spec="spp", genus="solidago", species="spp", family="asteraceae", # we should probably figure this one out and fix this
                         growthform="unkn", lifeform="f",origin="unkn", photo="unkn", Genus_Species="solidago_spp"))
    
### 2018

#meghan changed the path slightly
spcomp_2018 <- read.csv("Data 2018\\Species_Comp_k20a_before.csv") %>%
  left_join(cover_class_key, by="Cover_class") %>%
  dplyr::select(-Cover) %>%
  rename(Transect = Trasect) %>%
  mutate(Genus_Species = tolower(Genus_Species)) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="helianthis_annuas", "helianthus_annuus")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="rubus_spp", "rubus_occidentalis")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="rosa_spp", "rosa_arkansana")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="triaga_betonicifolia", "tragia_betonicifolia")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="mimosa_nuttallii", "mimosa_quadrivalvis")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="quincula_lobata", "quincula_lobata")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="erechtites_hieraciifolius", "erechtites_hieraciifolia")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="tridus_flavus", "tridens_flavus")) %>%
  left_join(dplyr::select(species_key, Genus_Species, code), by="Genus_Species") %>%
  rename(sp_code = code, Plot=Frame, Genus_Species_clean=Genus_Species) %>%
  filter(!Genus_Species_clean %in% c("rock","bare_ground","litter_herbeceous","litter_woody")) %>%
  mutate(Year=2018) %>%
  dplyr::select(Year, Transect, Plot, sp_code, Genus_Species_clean, abs_cover)

## Check for missing species
spcomp_2018 %>% filter(is.na(sp_code))

# Adding species numbers and checking 2018 data
# test <- spcomp_2018 %>% 
#   dplyr::filter(!Genus_Species %in% c("rock","bare_ground","litter_herbeceous","litter_woody")) %>%
#   left_join(
#     dplyr::select(species_key, code, Genus_Species),
#     by="Genus_Species")
# 
# sp_not_matching <- test %>%
#   dplyr::select(Genus_Species, code) %>%
#   filter(is.na(code)) %>%
#   unique()
  
### 2019 data
#meghan changed the path slightly
missing_spcode_2019 <- read.csv("Data 2019//Woody removal plots_K20A_Konza2019_200702_partial.csv") %>%
  rename(Genus_Species=Species) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  dplyr::select(Num_ID, Genus_Species) %>% unique() %>% filter(Num_ID=="")

spcomp_2019 <-  read.csv("Data 2019//Woody removal plots_K20A_Konza2019_200702_partial.csv") %>%
  rename(Genus_Species=Species, Year = "ï..Year", Plot=Plot..Frame.) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="liatris", "liatris_punctata")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="symphoricarpos", "symphoricarpos_orbiculatus")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="koleria_pyramidata", "koeleria_macrantha")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="ulmus", "ulmus_americana")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="lespedeza_virginica", "lespedeza_cuneata")) %>% ### Not positive but this could be lespedeza cuneata, which is on KNZ species list
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="pennstemon_sp", "penstemon_spp")) %>%
  mutate(Num_ID = replace(Num_ID, Num_ID=="KW PIC", "")) %>%
  mutate(Num_ID = as.numeric(Num_ID)) %>%
  left_join(dplyr::select(species_key, Genus_Species, code), by="Genus_Species") %>%
  rename(sp_code = code) %>%
  mutate(Num_ID = ifelse(Num_ID=="", sp_code, Num_ID)) %>%
  filter(!is.na(Num_ID)) %>%
  dplyr::select(-sp_code, -Notes) %>%
  rename(sp_code = Num_ID) %>%
  left_join(
    dplyr::select(species_key, Genus_Species, code) %>%
    rename(Genus_Species_clean=Genus_Species), 
    by=c("sp_code"="code")) %>%
  filter(!is.na(Genus_Species_clean)) %>%
  group_by(Year, Transect, Plot, sp_code, Genus_Species_clean) %>%
  summarize(June=max(June), August=max(August)) %>%
  ungroup() %>%
  mutate(abs_cover = pmax(June, August)) %>%
  dplyr::select(Year, Transect, Plot, sp_code, Genus_Species_clean, abs_cover)
  

### 2021 data
#meghan changed the path slightly
missing_spcode_2021 <- read.csv("Data 2021//Woody removal plots_K20A_Konza2021.csv") %>%
  rename(Genus_Species=Species) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  dplyr::select(Num_ID, Genus_Species) %>% unique() %>% filter(Num_ID=="")

spcomp_2021 <-  read.csv("Data 2021//Woody removal plots_K20A_Konza2021.csv") %>%
  rename(Genus_Species=Species, Year = "ï..Year", Plot=Plot..Frame.) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="liatris", "liatris_punctata")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="symphoricarpos", "symphoricarpos_orbiculatus")) %>% #
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="koleria_pyramidata", "koeleria_macrantha")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="ulmus", "ulmus_americana")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="lespedeza_virginica", "lespedeza_cuneata")) %>% ### Not positive but this could be lespedeza cuneata, which is on KNZ species list
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="solidago_sp.", "solidago_spp")) %>%
  mutate(Num_ID = replace(Num_ID, Num_ID=="KW PIC", "")) %>%
  mutate(Num_ID = as.numeric(Num_ID)) %>%
  left_join(dplyr::select(species_key, Genus_Species, code), by="Genus_Species") %>%
  rename(sp_code = code) %>%
  mutate(Num_ID = ifelse(Num_ID=="", sp_code, Num_ID)) %>%
  filter(!is.na(Num_ID)) %>%
  dplyr::select(-sp_code, -Notes) %>%
  rename(sp_code = Num_ID) %>%
  left_join(
    dplyr::select(species_key, Genus_Species, code) %>%
      rename(Genus_Species_clean=Genus_Species), 
    by=c("sp_code"="code")) %>%
  filter(!is.na(Genus_Species_clean)) %>%
  group_by(Year, Transect, Plot, sp_code, Genus_Species_clean) %>%
  summarize(June=max(June), August=max(August)) %>%
  ungroup() %>%
  mutate(abs_cover = pmax(June, August)) %>%
  dplyr::select(Year, Transect, Plot, sp_code, Genus_Species_clean, abs_cover)

spcomp_all <- spcomp_2018 %>%
  bind_rows(spcomp_2019, spcomp_2021) %>%
  rename(transect=Transect, plot=Plot, year=Year, genus_species=Genus_Species_clean) %>%
  full_join(plotcat, by=c("transect", "plot"))

}

###
### 3. Prep data for multivariate viewing and plot NMDS trajectories 
###
{
spcomp_wide <- spcomp_all %>%
  dplyr::select(-sp_code) %>%
  spread(key=genus_species, value=abs_cover) %>%
  replace(is.na(.),0)

# Separate out spcomp and environmental columns (cols are species) #
sp_data <- spcomp_wide %>% dplyr::select(-year:-category)
env_data <- spcomp_wide %>% dplyr::select(year:category)

mds_all <- metaMDS(sp_data)
mds_scores <- data.frame(env_data, scores(mds_all, display="sites"))

# ddply(irrt.scores.all, .(Year), summarize, ### Calculates centroids
# 	NMDS1 = mean(NMDS1),
# 	NMDS2 = mean(NMDS2))

## Analyses run seperately by year ##
ggplot(filter(mds_scores, category=="grassy"), aes(x=NMDS1, y=NMDS2, col=factor(year))) +
  geom_path(inherit.aes=F, aes(NMDS1, NMDS2, group=plot)) +
  geom_point() +
  facet_wrap(~ transect) +
  scale_shape_manual(values=c(1,16)) +
  theme_bw() 

### Mean NMDS
mds_means <- mds_scores %>%
  group_by(year, transect, category) %>%
  summarize_at(vars(NMDS1, NMDS2), .funs=list(mean=mean, sterr=SE_function))

ggplot(mds_means, aes(x=NMDS1_mean, y=NMDS2_mean, col=factor(year), shape=category, size=category, label=transect,
                      xmin=NMDS1_mean-NMDS1_sterr,xmax=NMDS1_mean+NMDS1_sterr,
                      ymin=NMDS2_mean-NMDS2_sterr,ymax=NMDS2_mean+NMDS2_sterr)) +
  geom_errorbarh(height=0, size=0.5) +
  geom_errorbar(width=0, size=0.5) +
  geom_point() +
  geom_text(size=3, nudge_x=0.05, nudge_y=0.05, col="black") +
#  facet_wrap(~ transect) +
  scale_shape_manual(values=c(1,16)) +
  scale_size_manual(values=c(1,3)) +
  theme_bw() 

ggplot(mds_means, aes(x=NMDS1_mean, y=NMDS2_mean, col=factor(year), shape=category,
                      xmin=NMDS1_mean-NMDS1_sterr,xmax=NMDS1_mean+NMDS1_sterr,
                      ymin=NMDS2_mean-NMDS2_sterr,ymax=NMDS2_mean+NMDS2_sterr)) +
  geom_path(inherit.aes=F, aes(NMDS1_mean, NMDS2_mean, group=category)) +
  geom_errorbarh(height=0) +
  geom_errorbar(width=0) +
  geom_point() +
  facet_wrap(~ transect) +
  scale_shape_manual(values=c(1,16)) +
  scale_size_manual(values=c(1,3)) +
  theme_bw() 

plot(mds_all, "species")

vf <- envfit(mds_all, sp_data, perm = 999)
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) %>%
  rename("Genus_Species"="Species") %>%
  left_join(species_key, by="Genus_Species")

ggplot(spp.scrs, aes(x=NMDS1, y=NMDS2, col=lifeform)) +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Genus_Species),
            size = 3) +
  theme_bw()

spp.scrs_sub <- filter(spp.scrs, !(NMDS1 < .15 & NMDS1 > -.15 & NMDS2 < .15 & NMDS2 > -.15))

ggplot(spp.scrs_sub, aes(x=NMDS1, y=NMDS2, col=lifeform)) +
  geom_segment(data = spp.scrs_sub,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs_sub, aes(x = NMDS1, y = NMDS2, label = Genus_Species),
            size = 3) +
  theme_bw()

}

###
### Look at functional group covers through time
###
{
  
  lifeform_sums <- spcomp_all %>%
    left_join(species_key, by=c("genus_species" = "Genus_Species")) %>%
    group_by(year, transect, plot, category, lifeform) %>%
    summarize(abs_cover = sum(abs_cover))
  
  
  ggplot(filter(lifeform_sums, lifeform %in% c("f","g","w")), aes(x=factor(year), y=abs_cover, col=category)) +
    geom_jitter(position=position_dodge(width=.7)) +
    geom_violin(alpha=.5) +
    facet_wrap(~lifeform) +
    theme_bw()
  
  filter(lifeform_sums, year==2018 & category=="woody" & abs_cover<50)
  filter(lifeform_sums, year==2018 & transect==5 & plot==3)
  filter(spcomp_all, year==2018 & transect==5 & plot==3)
  
  ?geom_jitter
}


###using codyn to look at changes through time.

spcomp_all2<-spcomp_all %>% 
  mutate(plotid=paste(transect, plot, sep="_"),
         trtid=paste(transect, category, sep="_"))

###within a transect how much are woody vs grassy plots changing over time?
cent_change<-multivariate_change(spcomp_all2, time.var="year", abundance.var="abs_cover", replicate.var="plotid", treatment.var = 'trtid', species.var="sp_code", reference.time=2018)

cent_changeplot<-cent_change %>% 
  separate(trtid, into=c("transect", "category"), sep="_") %>% 
  filter(year2==2021) %>% 
  group_by(category) %>% 
  summarize(change=mean(composition_change), se=SE_function(composition_change))

ggplot(data=cent_changeplot, aes(x=category, y=change))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=change-se, ymax=change+se), width=0.1, position=position_dodge())

##across the whole watershed how is dispersion changing?
cent_change2<-multivariate_change(spcomp_all2, time.var="year", abundance.var="abs_cover", replicate.var="plotid", treatment.var = 'category', species.var="sp_code", reference.time=2018)


















