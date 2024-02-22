### Calculate and visualize shrub trajectories post-burn
###
### Author: Meghan Avolio, Kevin Wilcox (wilcoxkr@gmail.com)
### Last updated June 7, 2022


### Set up workstation
setwd("C:\\Users\\mavolio2\\Dropbox/Konza Research/Shrub Islands/")
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\Konza projects_other\\ShrubRecolonization\\")

library(tidyverse)
library(codyn)
library(vegan)
library(gridExtra)
library(ggrepel)
library(lme4)
library(lmerTest)
library(emmeans)

theme_set(theme_bw(12))

### Standard Error function
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

###read in 2018 data
pre<-read.csv("Data 2018/Species_Comp_k20a_before.csv") %>%
  mutate(Plot=Frame-1)
post <- read.csv("Data 2018/Woody removal plots_K20A_after.csv") %>%
  rename(transect=Transect,
         plot=Plot..Frame.)
bareground<-post%>%
  filter(Species=="Other- Bare Ground") # there is so much bareground I am not sure this is going to work

transects<-read.csv("shrub transect locations.csv")

lf<-read.csv("Data 2018/pre_species_lifeform.csv")


# splist<-pre %>% 
#   select(Genus_Species) %>% 
#   unique()
# 
# write.csv(splist, "Data 2018/pre_species.csv", row.names = F)

plotslist<-pre%>%
  select(Trasect, Plot) %>% 
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
  group_by(Transect, Plot, Lifeform) %>% 
  left_join(cover_class_key, by="Cover_class") %>%
  summarize(abs_cover=sum(abs_cover)) %>% 
  ungroup() %>%
  group_by(Transect, Plot) %>%
  mutate(relcov=abs_cover/sum(abs_cover)) %>% 
  mutate(class=ifelse(Lifeform=="Woody"&abs_cover>40, 1, 0))%>%
  filter(class==1)%>%
  select(Transect, Plot, class) %>% 
  unique()

ggplot(data=subset(pre_class, Lifeform=="Woody"), aes(x=Plot, y=abs_cover))+
  geom_point()+
  geom_line()+
  facet_wrap(~Transect)+
  geom_hline(yintercept = 45)

testing <- pre_class %>%
  spread(key=Lifeform, value=abs_cover) %>%
  replace(is.na(.),0)
  
  
# Data frame identifying pre-burn grass versus woody states - even when correcting for plot shift, this doesn't match what we have observed in 2022
transect<-transects %>% 
  rename(Plot=plot)

plotcat<-plotslist%>%
  left_join(pre_class) %>% 
  mutate(category=ifelse(is.na(class), "grassy", 'woody')) %>% 
  rename(transect=Transect) %>% 
  select(-class) %>% 
  left_join(transect)

#how do classes differ in bareground?
bg<-bareground%>%
  left_join(plotcat)

boxplot(bg$Abundance~bg$category)  
summary(aov(Abundance~category, data = bg))#no sig diff



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
  #pre 2018 data had transects at plots at different points, Plot 1 is actually plot 0 and this was also done after shrubs were cut down, so I am not sure how much we want to use this.
spcomp_2018 <- read.csv("Data 2018/Species_Comp_k20a_before.csv") %>%
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
  dplyr::select(Year, Transect, Plot, sp_code, Genus_Species_clean, abs_cover) %>% 
  mutate(newplot=Plot-1) %>% 
  select(-Plot) %>% 
  rename(Plot=newplot) %>% 
  filter(Plot!=0)
  

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
missing_spcode_2019 <- read.csv("Data 2019/Woody removal plots_K20A_Konza2019_200702_all.csv") %>%
  rename(Genus_Species=Species) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  dplyr::select(Num_ID, Genus_Species) %>% unique() %>% filter(Num_ID=="")

spcomp_2019 <-  read.csv("Data 2019/Woody removal plots_K20A_Konza2019_200702_all.csv") %>%
  rename(Genus_Species=Species, Plot=Plot..Frame.) %>%
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
  select(-sp_code) %>%
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
  dplyr::select(Transect, Plot, sp_code, Genus_Species_clean, abs_cover) %>% 
  mutate(abs_cover=as.numeric(abs_cover), Year=2019)
  

### 2021 data
#meghan changed the path slightly
missing_spcode_2021 <- read.csv("Data 2021//Woody removal plots_K20A_Konza2021.csv") %>%
  rename(Genus_Species=Species) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  dplyr::select(Num_ID, Genus_Species) %>% unique() %>% filter(Num_ID=="")

spcomp_2021 <-  read.csv("Data 2021//Woody removal plots_K20A_Konza2021.csv") %>%
  rename(Genus_Species=Species, Plot=Plot..Frame.) %>%
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

spcomp_2022 <-  read.csv("Data 2022//Shrub Removal Species Comp 2022.csv") %>%
  select(-X) %>% 
  rename(Genus_Species=Species, Plot=Plot..Frame.) %>%
  mutate(Genus_Species = tolower(sub(" ", "_", Genus_Species))) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="liatris", "liatris_punctata")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="symphoricarpos", "symphoricarpos_orbiculatus")) %>% #
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="koleria_pyramidata", "koeleria_macrantha")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="ulmus", "ulmus_americana")) %>%
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="lespedeza_virginica", "lespedeza_cuneata")) %>% ### Not positive but this could be lespedeza cuneata, which is on KNZ species list
  mutate(Genus_Species = replace(Genus_Species, Genus_Species=="solidago_sp.", "solidago_spp")) %>%
  mutate(Num_ID = as.numeric(Num_ID)) %>%
  left_join(dplyr::select(species_key, Genus_Species, code), by="Genus_Species") %>%
  rename(sp_code = code) %>%
  mutate(Num_ID = ifelse(Num_ID=="", sp_code, Num_ID)) %>%
  filter(!is.na(Num_ID)) %>%
  dplyr::select(-sp_code) %>%
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
  bind_rows(spcomp_2019, spcomp_2021, spcomp_2022) %>%
  rename(transect=Transect, plot=Plot, year=Year, genus_species=Genus_Species_clean) %>%
  left_join(transects) %>% 
  rename(shrubSP=Notes)

write.csv(spcomp_all, "species_compositon_allyears.csv" )

###
### 3. Prep data for multivariate viewing and plot NMDS trajectories 
##MA note - after looking at the data, I do not think we shoudl include 2018, they are just different and we used different methods to get cover data. I no longer think they are comparable.
###
##also are only doing this for the 7 transects that we have all years of data for

#running each plot separate
spcomp_wide <- spcomp_all %>%
  filter(year!=2018, X2022.data=="yes") %>% 
  filter(shrubSP %in% c('', 'Cornus', 'Rhus glabra')) %>% #this drops 3 plots
  dplyr::select(-sp_code, -X2022.data) %>%
  spread(key=genus_species, value=abs_cover) %>%
  replace(is.na(.),0)

#taking the average for each transect unique data info
spcomp_wideave <- spcomp_all %>%
  filter(year!=2018, X2022.data=="yes") %>%
  filter(shrubSP %in% c('', 'Cornus', 'Rhus glabra')) %>% 
  dplyr::select(-sp_code, -X2022.data) %>% 
  group_by(year, transect, cover, shrubSP, genus_species) %>% 
  summarise(mcov=mean(abs_cover)) %>% 
  pivot_wider(names_from = genus_species, values_from = mcov, values_fill = 0) 

# Separate out spcomp and environmental columns (cols are species) #
sp_data <- spcomp_wideave %>% ungroup %>% dplyr::select(-year:-shrubSP)
env_data <- spcomp_wideave %>% ungroup %>% dplyr::select(year:shrubSP)

#this model is working. best solution not repeated
mds_all <- metaMDS(sp_data,trymax = 1000)
mds_scores <- data.frame(env_data, scores(mds_all, display="sites"))

# ddply(irrt.scores.all, .(Year), summarize, ### Calculates centroids
# 	NMDS1 = mean(NMDS1),
# 	NMDS2 = mean(NMDS2))

## Analyses run seperately by year ##
# ggplot(filter(mds_scores, cover=="Sh"), aes(x=NMDS1, y=NMDS2, col=factor(year))) +
#   geom_path(inherit.aes=F, aes(NMDS1, NMDS2, group=plot)) +
#   geom_point() +
#   facet_wrap(~ transect) +
#   scale_shape_manual(values=c(1,16)) +
#   theme_bw() 

### Mean NMDS
mds_means <- mds_scores %>%
  group_by(year, transect, cover, shrubSP) %>%
  summarize_at(vars(NMDS1, NMDS2), .funs=list(mean=mean, sterr=SE_function))

#this is the figure i think I want
a<-ggplot(mds_means, aes(x=NMDS1_mean, y=NMDS2_mean, col=as.factor(year),label=transect, shape=shrubSP)) +
  #geom_errorbarh(aes(xmin=NMDS1_mean-NMDS1_sterr,xmax=NMDS1_mean+NMDS1_sterr), height=0, size=0.5) +
  #geom_errorbar(aes(ymin=NMDS2_mean-NMDS2_sterr,ymax=NMDS2_mean+NMDS2_sterr), width=0, size=0.5) +
  geom_path(aes(group=transect), color="black")+
  geom_point(size=3) +
  facet_grid(~cover)+
  geom_text(size=3, nudge_x=0.05, nudge_y=0.05, col="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(mds_means, aes(x=NMDS1_mean, y=NMDS2_mean, col=cover, shape=as.factor(year),
#                       xmin=NMDS1_mean-NMDS1_sterr,xmax=NMDS1_mean+NMDS1_sterr,
#                       ymin=NMDS2_mean-NMDS2_sterr,ymax=NMDS2_mean+NMDS2_sterr)) +
#   geom_path(inherit.aes=F, aes(NMDS1_mean, NMDS2_mean, group=cover)) +
#   geom_errorbarh(height=0) +
#   geom_errorbar(width=0) +
#   geom_point(size=4) +
#   facet_wrap(~ transect) +
#   theme_bw() 
# 
# plot(mds_all, "species")
# 
# vf <- envfit(mds_all, sp_data, perm = 999)
# spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
# spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) %>%
#   rename("Genus_Species"="Species") %>%
#   left_join(species_key, by="Genus_Species")
# 
# ggplot(spp.scrs, aes(x=NMDS1, y=NMDS2, col=lifeform)) +
#   geom_segment(data = spp.scrs,
#                aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
#                arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#   geom_text(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Genus_Species),
#             size = 3) +
#   theme_bw()
# 
# spp.scrs_sub <- filter(spp.scrs, !(NMDS1 < .15 & NMDS1 > -.15 & NMDS2 < .15 & NMDS2 > -.15))
# 
# ggplot(spp.scrs_sub, aes(x=NMDS1, y=NMDS2, col=lifeform)) +
#   geom_segment(data = spp.scrs_sub,
#                aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
#                arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#   geom_text(data = spp.scrs_sub, aes(x = NMDS1, y = NMDS2, label = Genus_Species),
#             size = 3) +
#   theme_bw()





###using codyn to look at changes through time.

spcomp_all2<-spcomp_all %>% 
  filter(year!=2018, X2022.data=="yes") %>% 
  filter(shrubSP %in% c('', 'Cornus', 'Rhus glabra')) %>% 
  group_by(year, transect, cover, shrubSP) %>%
  mutate(plotid=paste(transect, plot, sep="_"),
         trtid=paste(transect, cover, shrubSP, sep="_"))

###within a transect how much are woody vs grassy plots changing over time?
cent_change<-multivariate_change(spcomp_all2, time.var="year", abundance.var="abs_cover", replicate.var="plotid", treatment.var = 'trtid', species.var="sp_code")


cent_changestats<-cent_change %>% 
  separate(trtid, into=c("transect", "cover", 'shrubSP'), sep="_") %>% 
  mutate(type=ifelse(shrubSP!="", shrubSP, cover))

m1<-lmer(composition_change~type + (1|transect), data=cent_changestats)
anova(m1)
emmeans(m1, ~type)

cent_changeplot<-cent_change %>% 
  separate(trtid, into=c("transect", "cover", 'shrubSP'), sep="_") %>% 
  group_by(cover, shrubSP) %>% 
  summarize(change=mean(composition_change), se=SE_function(composition_change), ci=se*1.96) %>% 
  mutate(type=ifelse(shrubSP!="", shrubSP, cover))

b<-ggplot(data=cent_changeplot, aes(x=type, y=change))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=change-ci, ymax=change+ci), width=0.1, position=position_dodge())+
  theme_bw()+
  scale_x_discrete(limits=c('G', 'T', 'Cornus', 'Rhus glabra'))+
  annotate('text', x=1, y=0.4, size=5, label="B")+
  annotate('text', x=2, y=0.55, size=5, label="A")+
  annotate('text', x=3, y=0.7, size=5, label="A")+
  annotate('text', x=4, y=0.5, size=5, label="AB")
b
grid.arrange(a,b)

##across the whole watershed how is dispersion changing?
cent_change2<-multivariate_change(spcomp_all2, time.var="year", abundance.var="abs_cover", replicate.var="plotid", treatment.var = 'category', species.var="sp_code", reference.time=2018)

###what is going on in the shrub patch communities?
species_key<-species_key %>% 
  rename(genus_species=Genus_Species)
shrubislands<-spcomp_all %>% 
  filter(year!=2018, X2022.data=="yes", cover=="Sh", abs_cover>0) %>% 
  filter(shrubSP %in% c('', 'Cornus', 'Rhus glabra')) %>% 
  group_by(year, transect, plot) %>% 
  mutate(rank=rank(-abs_cover, ties.method = 'random')) %>% 
  left_join(species_key) %>% 
  mutate(key=paste(transect, plot, sep="_")) %>% 
  mutate(lifeform2=ifelse(lifeform=="s"|lifeform=='g', 'graminoid',
                   ifelse(lifeform=='w', 'Woody', 'Forb'))) %>% 
  mutate(lifeform3=ifelse(rank<2, paste(gen, spec, sep="_"), ""))

labels<-shrubislands %>% 
  ungroup() %>% 
  select(key, shrubSP) %>% 
  unique()


ggplot(data=shrubislands, aes(x=rank, y=abs_cover, color=lifeform2, label=lifeform3))+
  geom_line(color='black', aes(group=year))+
  geom_text_repel(color='black')+
  geom_point(size=2, aes(shape=as.factor(year)))+
  facet_wrap(~key)+
  geom_text(data=labels, aes(x=Inf, y=Inf, label=shrubSP), hjust=1.5, vjust=1.5, color='black')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###compare with grass communities
grassyplots<-spcomp_all %>% 
  filter(year!=2018, X2022.data=="yes", cover=="G", abs_cover>0) %>% 
  group_by(year, transect, plot) %>% 
  mutate(rank=rank(-abs_cover, ties.method = 'random')) %>% 
  left_join(species_key) %>%
  mutate(key=paste(transect, plot, sep="_")) %>% 
  filter(key %in% c('10_1', '2_1', '5_10', '5_9', '7_2', '7_3', '7_7', '8_1', '9_8', '9_9', '9_10')) %>% 
  mutate(lifeform2=ifelse(lifeform=="s"|lifeform=='g', 'graminoid',
                          ifelse(lifeform=='w', 'Woody', 'Forb'))) %>% 
  mutate(lifeform3=ifelse(rank<2, paste(gen, spec, sep="_"), ""))

ggplot(data=grassyplots, aes(x=rank, y=abs_cover, color=lifeform2, label=lifeform3))+
  geom_line(color='black', aes(group=year))+
  geom_text_repel(color='black')+
  geom_point(size=2, aes(shape=as.factor(year)))+
  facet_wrap(~key)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##geom_text_repel()##basic analyses

###
### Look at functional group covers through time
###
{
  
  lifeform_sums <- spcomp_all %>%
    left_join(species_key, by=c("genus_species" = "Genus_Species")) %>%
    group_by(year, transect, plot, cover, lifeform) %>%
    summarize(abs_cover = sum(abs_cover))
  
  
ggplot(filter(lifeform_sums, lifeform %in% c("f","g","w")), aes(x=factor(year), y=abs_cover, col=cover)) +
    geom_jitter(position=position_dodge(width=.7)) +
    geom_violin(alpha=.5) +
    facet_wrap(~lifeform) +
    theme_bw()
  
meancov<-lifeform_sums %>% 
  group_by(year, lifeform, cover) %>% 
  summarize(mcov=mean(abs_cover), ster=SE_function(abs_cover)) %>% 
  filter(lifeform!="s"&lifeform!="m")

ggplot(data=meancov, aes(x=year, y=mcov, color=lifeform))+
  geom_errorbar(aes(ymin=mcov-ster, ymax=mcov+ster), width=0.5)+
  geom_point()+
  geom_line()+
  facet_wrap(~cover)+
  scale_color_manual(name="Lifeform", labels=c('Forb', 'Grass', 'Woody'), values = c('goldenrod2', 'darkgreen','brown'))

covertype<-spcomp_all %>% 
  select(transect, plot, cover) %>% 
  unique()

richeven<-community_structure(spcomp_all2, abundance.var = "abs_cover", time.var = "year", replicate.var = "plotid") %>% 
  separate(plotid, into = c('transect', 'plot')) %>% 
  mutate(transect=as.integer(transect), plot=as.integer(plot)) %>% 
  left_join(covertype)

meanrich<-richeven %>% 
  group_by(year, cover) %>% 
  summarize(mrich=mean(richness), ster=SE_function(richness))

ggplot(data=meanrich, aes(x=year, y=mrich, color=cover))+
  geom_errorbar(aes(ymin=mrich-ster, ymax=mrich+ster), width=0.5)+
  geom_point()+
  geom_line()+
  scale_color_manual(name="Plot type", labels=c('Grassy', 'Shrubby', 'Transitional'), values = c('darkgreen', 'brown','goldenrod2'))

}

















