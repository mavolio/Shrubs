setwd("C://Users/mavolio2/Dropbox/Konza Research/Shrub Islands/")

library(tidyverse)
library(codyn)


##1. for each transect, designate a plot as grassland, transitional, or burned, based on 2018 data
post<-read.csv("Data 2018/Woody removal plots_K20A_after.csv") %>% 
  rename(transect=Transect,
         plot=Plot..Frame.)

bareground<-post%>%
  filter(Species=="Other- Bare Ground") # there is so much bareground I am not sure this is going to work

hist(bareground$Abundance)

pre<-read.csv("Data 2018/Species_Comp_k20a_before.csv")

# splist<-pre %>% 
#   select(Genus_Species) %>% 
#   unique()
# 
# write.csv(splist, "Data 2018/pre_species.csv", row.names = F)
lf<-read.csv("Data 2018/pre_species_lifeform.csv")

plotslist<-pre%>%
  select(Trasect, Frame) %>% 
  unique()

pre_class_total<-pre %>% 
  left_join(lf) %>% 
  filter(Lifeform!="Drop")%>%
  group_by(Trasect, Frame) %>% 
  summarize(totcov=sum(Cover_class))
  
pre_class<-pre %>% 
  left_join(lf) %>% 
  filter(Lifeform!="Drop")%>%
  group_by(Trasect, Frame, Lifeform) %>% 
  summarize(cov=sum(Cover_class)) %>% 
  left_join(pre_class_total) %>% 
  mutate(relcov=cov/totcov) %>% 
  mutate(class=ifelse(Lifeform=="Woody"&relcov>0.55, 1, 0))%>%
  filter(class==1)%>%
  select(Trasect, Frame, class) %>% 
  unique()

plotcat<-plotslist%>%
  left_join(pre_class) %>% 
  mutate(category=ifelse(is.na(class), "grassy", 'woody')) %>% 
  rename(transect=Trasect,
         plot=Frame) %>% 
  select(-class)

#how do classes differ in bareground?
bg<-bareground%>%
  left_join(plotcat)

boxplot(bg$Abundance~bg$category)  
summary(aov(Abundance~category, data = bg))#no sig diff

##are these plots still woody in subsequent years?
d2019<-read.csv("Data 2019/Woody removal plots_K20A_Konza2019_200702_partial.csv") %>% 
  rename(transect=Transect,
         plot=Plot..Frame.)


