#Author: Giacomo Bignardi
#Modified: 2023-10-10
#
#Description: compute image reliability
#Program: Power analysis ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)

#clean working environment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

Discovery = read_csv(sprintf("%s/%s/01_Germine_2015/01_CTD_Germine2015.csv", wdOA,wdOA_output))%>% mutate(category = recode(category, FA_TOT = "faces", SC = "scenes", AO= "abstracts"))
Validation  = read_csv(sprintf("%s/%s/02_Sutherland_2020/01_CTD_Sutherland2020.csv", wdOA,wdOA_output))%>% mutate(category = recode(category, FA = "faces", SC = "scenes"))

#DISCOVERY####
#select only repeated images and one twin per pair
Discovery1 = Discovery %>%
  dplyr::select(-c(FamId_2 )) %>% 
  filter (is_repeated !=0,
          SibId == 1)
Discovery2 = Discovery %>%
  dplyr::select(-c(FamId_2 )) %>% 
  filter (is_repeated !=0,
          SibId == 2)

#one twin per pair intra-image reliability
intra_item_FA1 = Discovery1 %>% 
  filter(novariance_FA == 0,intraR_FA >= .5, category == "faces") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item_2,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))
intra_item_SC1 = Discovery1 %>% 
  filter(novariance_SC == 0,intraR_SC >= .5, category == "scenes") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item_2,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))
intra_item_AO1 = Discovery1 %>% 
  filter(novariance_AO == 0,intraR_AO >= .5, category == "abstracts") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item_2,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))


#other twin per pair intra-image reliability
intra_item_FA2 = Discovery2 %>% 
  filter(novariance_FA == 0,intraR_FA >= .5, category == "faces") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item_2,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))
intra_item_SC2 = Discovery2 %>% 
  filter(novariance_SC == 0,intraR_SC >= .5, category == "scenes") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item_2,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))
intra_item_AO2 = Discovery2 %>% 
  filter(novariance_AO == 0,intraR_AO >= .5, category == "abstracts") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item_2,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))


#prepare for plotting
intra_item1 = rbind(intra_item_FA1,intra_item_SC1,intra_item_AO1)
intra_item2 = rbind(intra_item_FA2,intra_item_SC2,intra_item_AO2)

p1 = intra_item1 %>%
  ggplot(aes(as.factor(category),intra_item)) +
  geom_boxplot() +
  ylim(.5,1) +
  labs(x = "domain", y = expression(R[xx-Image]))+
  theme_classic()
p2 = intra_item2 %>%
  ggplot(aes(as.factor(category),intra_item)) +
  geom_boxplot() +
  labs(x = "domain", y = expression(R[xx-Image]))+
  ylim(.5,1) +
  theme_classic()

#VALIDATION####
#select only repeated images and one twin per pair
Validation1 = Validation %>%
  filter (is_repeated !=0,
          twinN == "A")
Validation2 = Validation %>%
  filter (is_repeated !=0,
          twinN == "A")
#one twin per pair intra-image reliability
intra_item_FA1_val = Validation1 %>% 
  filter(novariance_FA == 0,intraR_FA >= .5, category == "faces") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))
intra_item_SC1_val = Validation1 %>% 
  filter(novariance_SC == 0,intraR_SC >= .5, category == "scenes") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))

#other twin per pair intra-image reliability
intra_item_FA2_val = Validation2 %>% 
  filter(novariance_FA == 0,intraR_FA >= .5, category == "faces") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))
intra_item_SC2_val = Validation2 %>% 
  filter(novariance_SC == 0,intraR_SC >= .5, category == "scenes") %>% 
  pivot_wider(values_from = Value, names_from = is_repeated) %>% 
  group_by(Item,category) %>% 
  rename(first = `1`,second = `2`) %>% 
  summarise(intra_item = cor(first,second))

#prepare for plotting
intra_item1_val = rbind(intra_item_FA1_val,intra_item_SC1_val)
intra_item2_val = rbind(intra_item_FA2_val,intra_item_SC2_val)

p1_val = intra_item1_val %>%
  ggplot(aes(as.factor(category),intra_item)) +
  geom_boxplot() +
  ylim(.5,1) +
  labs(x = "domain", y = expression(R[xx-Image]))+
  theme_classic()
p2_val = intra_item2_val %>%
  ggplot(aes(as.factor(category),intra_item)) +
  geom_boxplot() +
  labs(x = "domain", y = expression(R[xx-Image]))+
  ylim(.5,1) +
  theme_classic()

pdf(sprintf("%s/%s/00_F1b_VPC.pdf",wdOA,wdOA_ImageOutput),
    width = 10,
    height = 3.5)
F1b
dev.off()

#SF2####
png(sprintf("%s/%s/supplementary/00_SF2_image_reliability.png",wdOA,wdOA_ImageOutput),
    width = 4*325,
    height =4*325,
    unit = "px",
    res = 300,
    bg = "white")
((p1|p2)/
    (p1_val|p2_val)) + plot_layout(guides="collect") +plot_annotation(tag_levels = "a")
dev.off()