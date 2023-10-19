#Author: Giacomo Bignardi
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Description: 
#1. create two dataframe, one for twin one and one for twin 2, and then convert to wide format
#2. convert wide to long
#3. create a list of exclusion criteria then create a set of dataframe to be used in following analysis
# criteria: 
# ratings are all NAs or SD in ratings is = 0)
# Rxx-intra <.5 -> exclusion criteria based on Vessel et al. (2018, Cognition) and Chen et al. (2022, Current Biology)
# create separate df to be used for different types of analysis
#4. CTD -> to be used for following CTD modeling (contains all participants -> exclusion criteria = T)
#5. MLM -> to be used to fit the multilevel model to do VCA (contains ~90% of the participants -> exclusion criteria = T)
#6. MM2 -> to be used to compute mean minus 2 scores (contains ~90% of the participants -> exclusion criteria = T)
#Program: 01_prepare_df ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)
library(patchwork)
#ggpubr

#clean working enviroment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#load dataFrames:
Germine2015  = read_csv(sprintf("%s/Germine_2015/germine_CurrentBiology_twinaesthetics_PUBLIC.csv",wdOA_Data))
#load functions:
source(sprintf("%s/%s/functions/SDcheck.R", wdOA,wdOA_scripts)) #function to check SD
source(sprintf("%s/%s/functions/IntraRater.R", wdOA,wdOA_scripts)) #function to calculate Rxx-intra

####1. WIDE####
#Select FamId SibID Sex Age Objects, Scenes, faces (repeated first exp.) 
#images repeated (first exposure)= 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47
Twin1 = Germine2015 %>% 
  dplyr::select(...1:Zygosity.twin1, 
                OBJECTS1.twin1:OBJECTS65.twin1, 
                SCENES1.twin1:SCENES65.twin1,
                Faces.MIT1.twin1:Faces.MIT65.twin1,
                Faces.GENHEAD1.twin1:Faces.GENHEAD65.twin1, 
                Faces.GLASGOW1.twin1:Faces.GLASGOW65.twin1, 
                Faces.OTHER1.twin1:Faces.OTHER65.twin1) %>%
  dplyr::rename(FamId = ...1, SibId = Twin_Num_of2.twin1,Zygosity = Zygosity.twin1,Sex = sex_x.twin1, Age = age.twin1)%>%
  mutate(OBJECTS1.twin1 = as.numeric(OBJECTS1.twin1),
         SCENES1.twin1 = as.numeric(SCENES1.twin1),
         Faces.MIT1.twin1 = as.numeric(Faces.MIT1.twin1),
         Faces.GENHEAD1.twin1 = as.numeric(Faces.GENHEAD1.twin1),
         Faces.GLASGOW1.twin1 = as.numeric(Faces.GLASGOW1.twin1),
         Faces.OTHER1.twin1 = as.numeric(Faces.OTHER1.twin1))
names(Twin1) = sub('.twin1','',names(Twin1))
#note the FamID == 94 NA is a MZ female with all NAs

#twin2
Twin2 = Germine2015 %>% 
  dplyr::select(...1,Twin_Num_of2.twin2:Zygosity.twin2, 
                OBJECTS1.twin2:OBJECTS65.twin2, 
                SCENES1.twin2:SCENES65.twin2,
                Faces.MIT1.twin2:Faces.MIT65.twin2,
                Faces.GENHEAD1.twin2:Faces.GENHEAD65.twin2, 
                Faces.GLASGOW1.twin2:Faces.GLASGOW65.twin2, 
                Faces.OTHER1.twin2:Faces.OTHER65.twin2) %>%
  dplyr::rename(FamId = ...1, SibId = Twin_Num_of2.twin2,Zygosity = Zygosity.twin2,Sex = sex_x.twin2, Age = age.twin2)%>%
  mutate(OBJECTS1.twin2 = as.numeric(OBJECTS1.twin2),
         SCENES1.twin2 = as.numeric(SCENES1.twin2),
         Faces.MIT1.twin2 = as.numeric(Faces.MIT1.twin2),
         Faces.GENHEAD1.twin2 = as.numeric(Faces.GENHEAD1.twin2),
         Faces.GLASGOW1.twin2 = as.numeric(Faces.GLASGOW1.twin2),
         Faces.OTHER1.twin2 = as.numeric(Faces.OTHER1.twin2))
names(Twin2) = sub('.twin2','',names(Twin2))

#final wide dataframe
TwinWide = rbind(Twin1,Twin2)
#note that one pair of twins reported different ages (fam 182. Age is estimated to be their average)
TwinWide[TwinWide$FamId == 182,]$Age = mean(TwinWide[TwinWide$FamId == 182,]$Age)

#inspect complete NA for at least one domain:
TwinWide[(rowSums(is.na(TwinWide)) > 0),]
#338_1 miss ratings Glasgow
#418_2 miss ratings Face Mit 
#94_1 miss every ratings

###Save wide####
write_csv(TwinWide,sprintf("%s/%s/01_Germine_2015/01_TwinWide_Germine2015.csv",wdOA,wdOA_output))
ncol(TwinWide)

####2. LONG####
#Create long dataframe for each domain
#####long:object####
Obj_TwinLong = TwinWide %>% 
  dplyr::select(c(c(FamId:Zygosity),c(OBJECTS1:OBJECTS65))) %>%
  gather(Item, Value, OBJECTS1:OBJECTS65)%>%
  mutate(category ="AO",
         #IS the image repeated in the database ? 0 no, 1 yes (first impression), 2, yes (second impression)
         is_repeated= ifelse(Item %in% sprintf("OBJECTS%s",c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47)),1,
                             ifelse(Item %in% sprintf("OBJECTS%s",rep(51:65)),2,0))
  )

#####long:sceneries####
Sce_TwinLong = TwinWide %>% 
  dplyr::select(c(c(FamId:Zygosity),c(SCENES1:SCENES65))) %>%
  gather(Item, Value, SCENES1:SCENES65)%>%
  mutate(category ="SC",
         #IS the image repeated in the database ? 0 no, 1 yes (first impression), 2, yes (second impression)
         is_repeated= ifelse(Item %in% sprintf("SCENES%s",c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47)),1,
                             ifelse(Item %in% sprintf("SCENES%s",rep(51:65)),2,0))
  )

#####long:faces####
FacTot_TwinLong = TwinWide %>% 
  dplyr::select(c(c(FamId:Zygosity),c(Faces.MIT1:Faces.OTHER65))) %>%
  gather(Item, Value, Faces.MIT1:Faces.OTHER65)%>%
  mutate(category ="FA_TOT",
         #IS the image repeated in the database ? 0 no, 1 yes (first impression), 2, yes (second impression)
         is_repeated= ifelse(Item %in% sprintf("Faces.OTHER%s",c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47)),
                             1,
                             ifelse(Item %in% sprintf("Faces.OTHER%s",rep(51:65)),2,
                                    ifelse(Item %in% sprintf("Faces.GLASGOW%s",c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47)),
                                           1,
                                           ifelse(Item %in% sprintf("Faces.GLASGOW%s",rep(51:65)),2,
                                                  ifelse(Item %in% sprintf("Faces.GENHEAD%s",c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47)),
                                                         1,
                                                         ifelse(Item %in% sprintf("Faces.GENHEAD%s",rep(51:65)),2,
                                                                ifelse(Item %in% sprintf("Faces.MIT%s",c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43, 47)),
                                                                       1,
                                                                       ifelse(Item %in% sprintf("Faces.MIT%s",rep(51:65)),2,
                                                                              0))))))))
  )


TwinLong = rbind(FacTot_TwinLong,Obj_TwinLong,Sce_TwinLong)
TwinLong = TwinLong%>%mutate(FamId_2 = sprintf("%s_%s",FamId,SibId))
nrow(TwinLong)

####3.REMOVE####
#create a list of exclusion for particpants with sd = 0
TwinLong%>%filter(category == "SC")%>%select(Item,Value, FamId_2,is_repeated)
TwinLong%>%filter(category == "AO")%>%select(Item,Value, FamId_2,is_repeated)

#exclude participants with no variance (check funtctions/SDcheck for details on sdCheck)
sdCheck_Fac = sdCheck(TwinLong%>%filter(category == "FA_TOT")%>%select(Item,Value, FamId_2,is_repeated)) 
sdCheck_Sce = sdCheck(TwinLong%>%filter(category == "SC")%>%select(Item,Value, FamId_2,is_repeated)) 
sdCheck_Obj = sdCheck(TwinLong%>%filter(category == "AO")%>%select(Item,Value, FamId_2,is_repeated)) 

#print frequency of participant with sd = 0 per Zygosity type (390 is the total n of images)
table(TwinLong[TwinLong$FamId_2 %in% sdCheck_Fac$IDs,][,"Zygosity"])/390
table(TwinLong[TwinLong$FamId_2 %in% sdCheck_Sce$IDs,][,c('FamId',"Zygosity")])/390
table(TwinLong[TwinLong$FamId_2 %in% sdCheck_Obj$IDs,][,c('FamId',"Zygosity")])/390

#remove participant 94_1, which has only NA
TwinLong = TwinLong%>%filter(FamId_2!="94_1")

#create a dummy variable to store participant with sd == 0
TwinLong = TwinLong%>%
  mutate(novariance_FA = ifelse(FamId_2 %in% sdCheck_Fac$IDs,1,0),
         novariance_SC = ifelse(FamId_2 %in% sdCheck_Sce$IDs,1,0),
         novariance_AO = ifelse(FamId_2 %in% sdCheck_Obj$IDs,1,0))

#Additional removal: id 418 fam member 2 and id 338 fam member 1 per the domain == FA have only NA
TwinLong%>%filter(is.na(Value))%>%select(FamId_2,category)
TwinLong[TwinLong$FamId_2=="418_2",]$novariance_FA = rep(1,length(TwinLong[TwinLong$FamId_2=="418_2",]$novariance_FA))
TwinLong[TwinLong$FamId_2=="338_1",]$novariance_FA = rep(1,length(TwinLong[TwinLong$FamId_2=="338_1",]$novariance_FA))

#descriptives
table(TwinLong$Zygosity,TwinLong$Sex)/780
TwinLong%>%rstatix::get_summary_stats(Age)

#####save long####
#note that this creates a very long (and heavy) dataframe
write_csv(TwinLong,sprintf("%s/%s/01_Germine_2015/01_TwinLong_Germine2015.csv",wdOA,wdOA_output))

####4.CTD####
#tidy naming of Twin long
TwinLong = TwinLong%>%
  mutate(
    #name of domain
    domain_2 = ifelse(grepl("MIT", Item), "MIT", 
                      ifelse(grepl("GENHEAD", Item),"GENHEAD", 
                             ifelse(grepl("GLASGOW", Item),"GLASGOW",
                                    ifelse(grepl("OTHER", Item),"OTHER",category)))),
    #name of domain (as per original df)
    string = ifelse(domain_2 == "GENHEAD", str_sub(Item,1,13),
                    ifelse(domain_2 == "GLASGOW", str_sub(Item,1,13),
                           ifelse(domain_2 == "MIT", str_sub(Item,1,9),
                                  ifelse(domain_2 == "OTHER", str_sub(Item,1,11),
                                         ifelse(domain_2 == "AO",str_sub(Item,1,7),str_sub(Item,1,6)))))),
    #rename repeated item to be equal to the original image name
    Item_2 = ifelse(grepl("51",Item),sprintf("%s3",string),
                    ifelse(grepl("52",Item),sprintf("%s6",string),
                           ifelse(grepl("53",Item),sprintf("%s9",string),
                                  ifelse(grepl("54",Item),sprintf("%s12",string),
                                         ifelse(grepl("55",Item),sprintf("%s15",string),
                                                ifelse(grepl("56",Item),sprintf("%s18",string),
                                                       ifelse(grepl("57",Item),sprintf("%s21",string),
                                                              ifelse(grepl("58",Item),sprintf("%s24",string),
                                                                     ifelse(grepl("59",Item),sprintf("%s27",string),
                                                                            ifelse(grepl("60",Item),sprintf("%s30",string),
                                                                                   ifelse(grepl("61",Item),sprintf("%s33",string),
                                                                                          ifelse(grepl("62",Item),sprintf("%s36",string),
                                                                                                 ifelse(grepl("63",Item),sprintf("%s39",string),
                                                                                                        ifelse(grepl("64",Item),sprintf("%s43",string),
                                                                                                               ifelse(grepl("65",Item),sprintf("%s47",string),Item)))))))))))))))
  )%>%
  select(-c(string,Item))

#####IntraRater####
#intra participant correlations: reliability (only on repeated images)
#Exclusion criteria is set as in Vessel et al. 2018 (intraR <.5)
intraR_Twin = TwinLong%>%filter(is_repeated ==1 | is_repeated==2)

#prepare dataframe for within reliability analysis                                                                                        
Fac_intraR_Twin = intraR_Twin%>%
  select(Item_2,Value,FamId_2,is_repeated,category,SibId,novariance_FA)%>%
  filter(category=="FA_TOT")%>%
  select(-c(SibId,novariance_FA,category))%>%
  mutate(Value= as.numeric(Value))

Sce_intraR_Twin = intraR_Twin%>%
  select(Item_2,Value,FamId_2,is_repeated,category,SibId,novariance_SC)%>%
  filter(category=="SC")%>% ## 4 removed
  select(-c(SibId,novariance_SC,category))%>%
  mutate(Value= as.numeric(Value))

Obj_intraR_Twin = intraR_Twin%>%
  select(Item_2,Value,FamId_2,is_repeated,category,SibId,novariance_AO)%>%
  filter(category=="AO")%>%
  select(-c(SibId,novariance_AO,category))%>%
  mutate(Value = as.numeric(Value))


#Calculate Rxx-intra rater correlation (reliability score to use for later exclusion)
Fac_intraR_Twin_output = intraRater(Fac_intraR_Twin)
Sce_intraR_Twin_output = intraRater(Sce_intraR_Twin)
Obj_intraR_Twin_output = intraRater(Obj_intraR_Twin)
#warnings refer to sd == 0, which are handled by the intraRater function

#number of excluded participants
sum(Fac_intraR_Twin_output[[1]]$intraR<.5, na.rm = T)
sum(Sce_intraR_Twin_output[[1]]$intraR<.5, na.rm = T)
sum(Obj_intraR_Twin_output[[1]]$intraR<.5, na.rm = T)

#_Figure:SF1_####
#To vizualize exclusion criteria:
p1_exclusion_fac = ggpubr::gghistogram(Fac_intraR_Twin_output[[1]], "intraR") + geom_vline(xintercept = .5) +xlab("intra-rater r faces (main)")
p1_exclusion_pla = ggpubr::gghistogram(Sce_intraR_Twin_output[[1]], "intraR") + geom_vline(xintercept = .5) +xlab("intra-rater r scenes (main)")
p1_exclusion_abs = ggpubr::gghistogram(Obj_intraR_Twin_output[[1]], "intraR") + geom_vline(xintercept = .5) +xlab("intra-rater r abstract (main)")
save(p1_exclusion_fac,p1_exclusion_pla,p1_exclusion_abs,file = sprintf("%s/%s/01_Germine_2015/01_sup_FS1_exclusion_Germine2015.Rdata",wdOA,wdOA_output))

#rename in order to merge:
names(Fac_intraR_Twin_output[[1]]) = c("FamId_2","intraR")
names(Sce_intraR_Twin_output[[1]]) = c("FamId_2","intraR")
names(Obj_intraR_Twin_output[[1]]) = c("FamId_2","intraR")

#save general df
TwinLong = merge(TwinLong,Fac_intraR_Twin_output[[1]][1:2],by = "FamId_2", all = T)
TwinLong = merge(TwinLong,Sce_intraR_Twin_output[[1]][1:2], by = "FamId_2", all = T)
TwinLong= merge(TwinLong,Obj_intraR_Twin_output[[1]][1:2], by = "FamId_2", all = T)
TwinLong = TwinLong%>%dplyr::rename(intraR_FA = "intraR.x",intraR_SC = "intraR.y",intraR_AO = "intraR")

####save CTD####
write_csv(TwinLong,sprintf("%s/%s/01_Germine_2015/01_CTD_Germine2015.csv",wdOA,wdOA_output))

##_Table 1_####
#all
TwinLong%>%
  #apply filtering
  select(Zygosity,Sex, SibId,FamId, FamId_2, Age)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair,Age)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(Zygosity,Sex), mean(Age), range(Age))

TwinLong%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2)%>%
  distinct()%>%
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(Zygosity,Sex))

#faces
TwinLong%>%
  #apply filtering
  filter(novariance_FA == 0)%>%
  filter(intraR_FA >.5)%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2, Age)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair,Age)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(Zygosity,Sex), mean(Age), range(Age))

TwinLong%>%
  filter(novariance_FA == 0)%>%
  filter(intraR_FA >.5)%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2)%>%
  distinct()%>%
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(Zygosity,Sex))

#scenes
TwinLong%>%
  #apply filtering
  filter(novariance_SC == 0)%>%
  filter(intraR_SC >.5)%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2,Age)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair,Age)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(Zygosity,Sex), mean(Age), range(Age))

TwinLong%>%
  filter(novariance_SC == 0)%>%
  filter(intraR_SC >.5)%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2)%>%
  distinct()%>%
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(Zygosity,Sex))

#abstracts
TwinLong%>%
  #apply filtering
  filter(novariance_AO == 0)%>%
  filter(intraR_AO >.5)%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2,Age)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair,Age)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(Zygosity,Sex), mean(Age), range(Age))

TwinLong%>%
  filter(novariance_AO == 0)%>%
  filter(intraR_AO >.5)%>%
  select(Zygosity,Sex, SibId,FamId, FamId_2)%>%
  distinct()%>%
  mutate(Pair = ifelse(duplicated(FamId)| duplicated(FamId, fromLast = TRUE),2,1))%>%
  select(Zygosity,Sex,FamId,Pair)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(Zygosity,Sex))

####descriptives: Intra-rater####
#compute once more intrarater but after exclusion criteria
intraRater(TwinLong%>%
             filter(is_repeated ==1 | is_repeated==2)%>%
             select(Item_2,Value,FamId_2,is_repeated,category,SibId,novariance_FA,intraR_FA)%>%
             filter(category=="FA_TOT")%>%  
             filter(intraR_FA>.5)%>%
             filter(novariance_FA==0)%>%
             select(-c(SibId,novariance_FA,category,intraR_FA))%>%
             mutate(Value= as.numeric(Value)))
intraRater(TwinLong%>%
             filter(is_repeated ==1 | is_repeated==2)%>%
             select(Item_2,Value,FamId_2,is_repeated,category,SibId,novariance_SC,intraR_SC)%>%
             filter(category=="SC")%>%  
             filter(intraR_SC>.5)%>%
             filter(novariance_SC==0)%>%
             select(-c(SibId,novariance_SC,category,intraR_SC))%>%
             mutate(Value= as.numeric(Value)))
intraRater(TwinLong%>%
             filter(is_repeated ==1 | is_repeated==2)%>%
             select(Item_2,Value,FamId_2,is_repeated,category,SibId,novariance_AO,intraR_AO)%>%
             filter(category=="AO")%>%  
             filter(intraR_AO>.5)%>%
             filter(novariance_AO==0)%>%
             select(-c(SibId,novariance_AO,category,intraR_AO))%>%
             mutate(Value= as.numeric(Value)))



####5.MLM####
#prepare df for MLM (VCA) for twin 1 and 2
#####twin1####
#create a dummy dataframe with repeated images containing NA for non-repeated images
MLM_Image_1 = TwinLong%>%
  filter(is_repeated ==0, SibId ==1)%>%
  dplyr::mutate(is_repeated = 1)%>%
  dplyr::select(Item_2,Value,FamId_2,Sex, Age,is_repeated,category,Zygosity,novariance_FA,novariance_SC, novariance_AO)
MLM_Image_1na = TwinLong%>%
  filter(is_repeated ==0,SibId ==1)%>%
  dplyr::mutate(Value= NA,is_repeated =2 )%>%
  dplyr::select(Item_2,Value,FamId_2,Sex, Age,is_repeated,category,Zygosity,novariance_FA,novariance_SC, novariance_AO)

MLM_Image_Twin1 = rbind(MLM_Image_1, MLM_Image_1na)

#merge the dataframe with the dataframe containing only repeated images
MLM_Twin1 = rbind(intraR_Twin%>%filter(SibId ==1)%>%select(Item_2,Value,FamId_2,is_repeated,Sex, Age,category,Zygosity, novariance_FA,novariance_SC, novariance_AO),MLM_Image_Twin1)

######faces####
#create a dataframe for MLM faces
Fac_MLM_Twin1 = MLM_Twin1%>%
  filter(category=="FA_TOT")%>%
  merge(Fac_intraR_Twin_output$intraRater[1:2], by = "FamId_2", all.x = T)

#sample size per zygosity before exclusion (Descriptives)
table(Fac_MLM_Twin1$Zygosity,Fac_MLM_Twin1$Sex)/400
#apply exclusion criteria
Fac_MLM_Twin1 = Fac_MLM_Twin1%>%
  filter(intraR>=.5, novariance_FA==0)
#tidylog: filter: removed 7,200 rows (2%), 302,000 rows remaining
#sample size per zygosity after exclusion (Descriptives)
table(Fac_MLM_Twin1$Zygosity,Fac_MLM_Twin1$Sex)/400

#entry as requested in the VCA function
Fac_MLM_Twin1 = Fac_MLM_Twin1%>%
  select(Item_2,Value,FamId_2,is_repeated)%>%
  mutate(FamId = substr(FamId_2,1,nchar(FamId_2)-2))%>%
  select(Item_2,Value,FamId,is_repeated) 

######sceneries####
#create a dataframe for MLM faces
Sce_MLM_Twin1 = MLM_Twin1%>%
  filter(category=="SC")%>%
  merge(Sce_intraR_Twin_output$intraRater[1:2], by = "FamId_2", all.x = T)

#sample size per zygosity before exclusion (Descriptives)
table(Sce_MLM_Twin1$Zygosity,Sce_MLM_Twin1$Sex)/100
Sce_MLM_Twin1 = Sce_MLM_Twin1%>%
  filter(intraR>=.5, novariance_SC==0)
#log: filter: removed 1,100 rows (1%), 76,200 rows remaining
#sample size per zygosity after exclusion (Descriptives)
table(Sce_MLM_Twin1$Zygosity,Sce_MLM_Twin1$Sex)/100

#entry as requested in the VCA toolbox
Sce_MLM_Twin1 = Sce_MLM_Twin1%>%
  select(Item_2,Value,FamId_2,is_repeated)%>%
  mutate(FamId = substr(FamId_2,1,nchar(FamId_2)-2))%>%
  select(Item_2,Value,FamId,is_repeated) 

######abstract####
#create a dataframe for MLM faces
Obj_MLM_Twin1 = MLM_Twin1%>%
  filter(category=="AO")%>%
  merge(Obj_intraR_Twin_output$intraRater[1:2], by = "FamId_2", all.x = T)

#sample size per zygosity before exclusion (Descriptives)
table(Obj_MLM_Twin1$Zygosity,Obj_MLM_Twin1$Sex)/100
Obj_MLM_Twin1 = Obj_MLM_Twin1%>%
  filter(intraR>=.5, novariance_AO==0)
#log: filter: removed 7,000 rows (9%), 70,300 rows remaining
#sample size per zygosity after exclusion (Descriptives)
table(Obj_MLM_Twin1$Zygosity,Obj_MLM_Twin1$Sex)/100

#entry as requested in the VCA toolbox
Obj_MLM_Twin1 = Obj_MLM_Twin1%>%
  select(Item_2,Value,FamId_2,is_repeated)%>%
  mutate(FamId = substr(FamId_2,1,nchar(FamId_2)-2))%>%
  select(Item_2,Value,FamId,is_repeated) 

#####twin2####
#create a dummy dataframe with repeated images containing NA for non-repeated images
MLM_Image_2a = TwinLong%>%
  filter(is_repeated ==0, SibId ==2)%>%
  dplyr::mutate(is_repeated = 1)%>%
  dplyr::select(Item_2,Value,FamId_2,Sex, Age,is_repeated,category,Zygosity,novariance_FA,novariance_SC, novariance_AO)
MLM_Image_2b = TwinLong%>%
  filter(is_repeated ==0,SibId ==2)%>%
  dplyr::mutate(Value= NA,is_repeated =2 )%>%
  dplyr::select(Item_2,Value,FamId_2,Sex, Age,is_repeated,category,Zygosity,novariance_FA,novariance_SC, novariance_AO)

MLM_Image_Twin2 = rbind(MLM_Image_2a, MLM_Image_2b)

#merge the dataframe with the dataframe containing only repeated images
MLM_Twin2 = rbind(intraR_Twin%>%filter(SibId ==2)%>%select(Item_2,Value,FamId_2,is_repeated,Sex, Age,category, Zygosity,novariance_FA,novariance_SC, novariance_AO),MLM_Image_Twin2)

######faces####
#create a dataframe for MLM faces
Fac_MLM_Twin2 = MLM_Twin2%>%
  filter(category=="FA_TOT")%>%
  merge(Fac_intraR_Twin_output$intraRater[1:2], by = "FamId_2", all.x = T)

#sample size per zygosity before exclusion (Descriptives)
table(Fac_MLM_Twin2$Zygosity,Fac_MLM_Twin2$Sex)/400
Fac_MLM_Twin2 = Fac_MLM_Twin2%>%
  filter(intraR>=.5, novariance_FA==0)#18 removed
#filter: removed 6,800 rows (2%), 302,800 rows remaining
#sample size per zygosity after exclusion (Descriptives)
table(Fac_MLM_Twin2$Zygosity,Fac_MLM_Twin2$Sex)/400

#entry as requested in the VCA toolbox
Fac_MLM_Twin2 = Fac_MLM_Twin2%>%
  select(Item_2,Value,FamId_2,is_repeated)%>%
  mutate(FamId = substr(FamId_2,1,nchar(FamId_2)-2))%>%
  select(Item_2,Value,FamId,is_repeated) 

######sceneries####
#create a dataframe for MLM faces
Sce_MLM_Twin2 = MLM_Twin2%>%
  filter(category=="SC")%>%
  merge(Sce_intraR_Twin_output$intraRater[1:2], by = "FamId_2", all.x = T)

#sample size per zygosity before exclusion (Descriptives)
table(Sce_MLM_Twin2$Zygosity,Sce_MLM_Twin2$Sex)/100
Sce_MLM_Twin2 = Sce_MLM_Twin2%>%
  filter(intraR>=.50, novariance_SC==0)
#filter: removed 1,100 rows (1%), 76,300 rows remaining
#sample size per zygosity after exclusion (Descriptives)
table(Sce_MLM_Twin2$Zygosity,Sce_MLM_Twin2$Sex)/100

#entry as requested in the VCA toolbox
Sce_MLM_Twin2 = Sce_MLM_Twin2%>%
  select(Item_2,Value,FamId_2,is_repeated)%>%
  mutate(FamId = substr(FamId_2,1,nchar(FamId_2)-2))%>%
  select(Item_2,Value,FamId,is_repeated) 

######abstract####
#create a dataframe for MLM faces
Obj_MLM_Twin2 = MLM_Twin2%>%
  filter(category=="AO")%>%
  merge(Obj_intraR_Twin_output$intraRater[1:2], by = "FamId_2", all.x = T)

#sample size per zygosity before exclusion (Descriptives)
table(Obj_MLM_Twin2$Zygosity,Obj_MLM_Twin2$Sex)/100
Obj_MLM_Twin2 = Obj_MLM_Twin2%>%
  filter(intraR>=.50, novariance_AO==0)
#filter: removed 6,300 rows (8%), 71,100 rows remaining
#sample size per zygosity after exclusion (Descriptives)
table(Obj_MLM_Twin2$Zygosity,Obj_MLM_Twin2$Sex)/100

#entry as requested in the VCA toolbox
Obj_MLM_Twin2 = Obj_MLM_Twin2%>%
  select(Item_2,Value,FamId_2,is_repeated)%>%
  mutate(FamId = substr(FamId_2,1,nchar(FamId_2)-2))%>%
  select(Item_2,Value,FamId,is_repeated) 

#save for MLM####
write_csv(Fac_MLM_Twin1,sprintf("%s/%s/01_Germine_2015/01_MLM_Fac_Twin1_Germine2015.csv",wdOA,wdOA_output))
write_csv(Sce_MLM_Twin1,sprintf("%s/%s/01_Germine_2015/01_MLM_Sce_Twin1_Germine2015.csv",wdOA,wdOA_output))
write_csv(Obj_MLM_Twin1,sprintf("%s/%s/01_Germine_2015/01_MLM_Obj_Twin1_Germine2015.csv",wdOA,wdOA_output))
write_csv(Fac_MLM_Twin2,sprintf("%s/%s/01_Germine_2015/01_MLM_Fac_Twin2_Germine2015.csv",wdOA,wdOA_output))
write_csv(Sce_MLM_Twin2,sprintf("%s/%s/01_Germine_2015/01_MLM_Sce_Twin2_Germine2015.csv",wdOA,wdOA_output))
write_csv(Obj_MLM_Twin2,sprintf("%s/%s/01_Germine_2015/01_MLM_Obj_Twin2_Germine2015.csv",wdOA,wdOA_output))

#6. MM2####
#entry for mm2 as required from the MM2 function
Fac_MLM_Twin1$FamId = sprintf("%s_1",Fac_MLM_Twin1$FamId)
Fac_MLM_Twin2$FamId = sprintf("%s_2",Fac_MLM_Twin2$FamId)
Sce_MLM_Twin1$FamId = sprintf("%s_1",Sce_MLM_Twin1$FamId)
Sce_MLM_Twin2$FamId = sprintf("%s_2",Sce_MLM_Twin2$FamId)
Obj_MLM_Twin1$FamId = sprintf("%s_1",Obj_MLM_Twin1$FamId)
Obj_MLM_Twin2$FamId = sprintf("%s_2",Obj_MLM_Twin2$FamId)

Fac_MLM_Twin = rbind(Fac_MLM_Twin1, Fac_MLM_Twin2)
Sce_MLM_Twin = rbind(Sce_MLM_Twin1, Sce_MLM_Twin2)
Obj_MLM_Twin = rbind(Obj_MLM_Twin1, Obj_MLM_Twin2)

Fac_MLM_Twin = Fac_MLM_Twin%>%mutate(domain = "FA_TOT")
Sce_MLM_Twin = Sce_MLM_Twin%>%mutate(domain = "SC")
Obj_MLM_Twin = Obj_MLM_Twin%>%mutate(domain = "AO")

MM2_Twin = rbind(Fac_MLM_Twin,Sce_MLM_Twin,Obj_MLM_Twin)

#save for mm2####
write_csv(MM2_Twin,sprintf("%s/%s/01_Germine_2015/01_MM2_Germine2015.csv",wdOA,wdOA_output))
