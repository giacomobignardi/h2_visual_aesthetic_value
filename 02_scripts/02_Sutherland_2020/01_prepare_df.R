#Author: Giacomo Bignardi
#Adapted from: Sutherland et al., 2020 PNAS
#Date: 28-04-2021
#Last modified: 19-09-2023
#
#Description: repeate procedure highlighted in 01_Germine_2015/01_prepare_df.R
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
Sutherland2020  = read_csv(sprintf("%s/%s/Sutherland_2020/01_Twin_dataScored.csv", wdOA,wdOA_Data))
#load functions:
source(sprintf("%s/%s/functions/SDcheck.R", wdOA,wdOA_scripts))
source(sprintf("%s/%s/functions/IntraRater.R", wdOA,wdOA_scripts))

#####1. WIDE####
#following avaiable code from Sutherland et al.
#import exclusion information based on participant TEXT ENTRY (e.g. technical difficulties):
ExcList1  = read_csv(sprintf("%s/%s/Sutherland_2020/0_Twin_ZygosityEstimated.csv", wdOA,wdOA_Data))%>%
  select(c("refcode1","matchzygosity","TRA_Q_match1","TRA_Q_match2"))%>%
  filter(matchzygosity == 0 | TRA_Q_match1 == 0 | TRA_Q_match1 == 0)%>%
  select(refcode1)
ExcList2 = read_csv(sprintf("%s/%s/Sutherland_2020/0_Twin_ExclusionList.csv", wdOA,wdOA_Data))%>%
  filter(Exclude == 1)%>%
  select(refCode)
# #load functions:
# source(sprintf("%s/%s/functions/xxx.R", wdOA,wdOA_scripts))
#Apply exclusion criteria:
Sutherland2020 = Sutherland2020%>%
  filter(!(refCode %in% ExcList1$refcode1))%>% #49 exclusions
  filter(!(refCode %in% ExcList2$refCode))%>%
  mutate(twinN = substring(refCode,7,7)) #create a column with birth order

#Select FamId SibID Sex Age Objects, Scenes, faces (repeated first exp.)
Twin1 = Sutherland2020 %>% 
  filter(twinN == "A",TwinType != "UN") %>% #remove additional uknown zygosity
  dplyr::select(refCode,
                TwinType,
                twinN,
                ZygosityStatus,
                Sex,
                AgeAtTestingCalculated,
                ageMatch,
                Anthony_Rimmer_11_oval_1:sc194_1)

Twin2 = Sutherland2020 %>% 
  filter(twinN == "B",TwinType != "UN") %>%
  dplyr::select(refCode,
                TwinType,
                twinN,
                ZygosityStatus,
                Sex,
                AgeAtTestingCalculated,ageMatch,
                Anthony_Rimmer_11_oval_1:sc194_1)

TwinWide = rbind(Twin1,Twin2)
###Save wide####
write_csv(TwinWide,sprintf("%s/%s/02_Sutherland_2020/01_TwinWide_Sutherland2020.csv",wdOA,wdOA_output))

####2. LONG####

#FROM SUTHERLAND ET AL. 2020-------------------------------------------------------------------- #
#following avaiable code from Sutherland et al. 2020 (https://osf.io/35zf8/?view_only=e76c6755dcea4be2adc5b075cae896e8)
#select ratitngs 
f1 = match("Anthony_Rimmer_11_oval_1",names(TwinWide))
allRatings = TwinWide[,f1:(f1+523)]
match("sc194_1",names(TwinWide))
#locate column number of first stimulus from each test
a1 = match("Anthony_Rimmer_11_oval_1",names(TwinWide))
d1 = match("Aaron_Mink_9_oval_1", names(TwinWide))
s1 = match("sc001_1", names(TwinWide))

#break up the ratings into the different tests
attrRatings = TwinWide[a1:(a1+149)]
domRatings = TwinWide[d1:(d1+149)]
sceneRatings = TwinWide[s1:(s1+73)]

#####faces####
#Convert to long ATTR Faces
Fac_TwinLong = TwinWide %>% 
  select(c(
    refCode,
    TwinType,
    twinN,
    ZygosityStatus,
    Sex,
    AgeAtTestingCalculated,
    colnames(attrRatings[,1]):colnames(attrRatings[,length(attrRatings)]))) %>%
  gather(Item, Value, Anthony_Rimmer_11_oval_1:`Google_1_Wayne.Waller_8_oval_1`)%>%
  mutate(category = "FA")%>%
  mutate(is_repeated = ifelse(Item %in% colnames(attrRatings[,1:50]),1,ifelse(Item %in% colnames(attrRatings[,51:100]),2,0)))

#####scenes####
#Convert to long ATTR Scenes
Sce_TwinLong = TwinWide %>% 
  select(c(
    refCode,
    TwinType,
    twinN,
    ZygosityStatus,
    Sex,
    AgeAtTestingCalculated,#Exclude,
    colnames(sceneRatings[,1]):colnames(sceneRatings[,length(sceneRatings)]))) %>%
  gather(Item, Value, sc001_1:sc194_1)%>%
  mutate(category = "SC")%>%
  mutate(is_repeated = ifelse(Item %in% colnames(sceneRatings[,1:24]),1,ifelse(Item %in% colnames(sceneRatings[,25:48]),2,0)))

#####control####
#Convert to long DOMIN Scenes (CONTROL)
Control_long = TwinWide %>% 
  select(c(
    refCode,
    TwinType,
    twinN,
    ZygosityStatus,
    Sex,
    AgeAtTestingCalculated,#Exclude,
    colnames(domRatings[,1]):colnames(domRatings[,length(domRatings)]))) %>%
  gather(Item, Value, Aaron_Mink_9_oval_1:`Google_1_Wilma.Dixon_15_oval_1`) %>%
  mutate(category = "DO")%>%
  mutate(is_repeated = ifelse(Item %in% colnames(domRatings[,1:50]),1,ifelse(Item %in% colnames(domRatings[,51:100]),2,0)))

#concatenate df together
TwinLong = rbind(Fac_TwinLong,Sce_TwinLong,Control_long)
#resolve problem for future csv loading --- expected: 1/0/T/F/TRUE/FALSE 
TwinLong = TwinLong%>%mutate(Sex = ifelse(Sex == 'M',1 ,2),
                             Item = substr(Item,1,nchar(Item)-2))
TwinLong = TwinLong%>%mutate(refCode_2 = substr(refCode, 1,nchar(refCode)-2))


####3.REMOVE####
#exclude participants with no variance (check funtctions/SDcheck for details on sdCheck)
sdCheck_Fac = sdCheck(TwinLong%>%select(Item,Value,refCode,category,is_repeated)%>%filter(category=="FA")%>%select(-category))
sdCheck_Sce = sdCheck(TwinLong%>%select(Item,Value,refCode,category,is_repeated)%>%filter(category=="SC")%>%select(-category))
sdCheck_Con = sdCheck(TwinLong%>%select(Item,Value,refCode,category,is_repeated)%>%filter(category=="DO")%>%select(-category))

#print frequency of participant with sd = 0 per Zygosity type (390 is the total n of images, 374 is the tot number of images per participants across domains and control)
table(TwinLong[TwinLong$refCode %in% sdCheck_Fac$IDs,][,c("TwinType")])/374
table(TwinLong[TwinLong$refCode %in% sdCheck_Sce$IDs,][,c("TwinType")])/374
table(TwinLong[TwinLong$refCode %in% sdCheck_Con$IDs,][,c("TwinType")])/374

#create a dummy variable to store if participant sd == 0
TwinLong = TwinLong%>%
  mutate(novariance_FA = ifelse(refCode %in% sdCheck_Fac$IDs,1,0),
         novariance_SC = ifelse(refCode %in% sdCheck_Sce$IDs,1,0),
         novariance_DO = ifelse(refCode %in% sdCheck_Con$IDs,1,0))

#descriptives
#n of pairs
TwinLong%>%
  #apply filtering
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(TwinType,Sex), mean(AgeAtTestingCalculated), range(AgeAtTestingCalculated), sd(AgeAtTestingCalculated))

#n of coplete pairs
TwinLong%>%
  #apply filtering
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  #tabulate descritives
  filter(Pair == 2)%>%
  summarise(table(TwinType,Sex))

#####save long####
write_csv(TwinLong,sprintf("%s/%s/02_Sutherland_2020/01_TwinLong_Sutherland2020.csv",wdOA,wdOA_output))

####4.CTD####
#####IntraRater####
#intra participant correlations: reliability (only on repeated images)
#Exclusion criteria is set as in Vessel et al. 2018 (intraR <.5)
intraR_Twin = TwinLong%>%filter(is_repeated ==1 | is_repeated==2)

#prepare dataframe for within reliability analysis                                                                                        
Fac_intraR_Twin = intraR_Twin%>%
  select(Item,Value,refCode,is_repeated,category,twinN,novariance_FA)%>%
  filter(category=="FA")%>%
  select(-c(twinN,novariance_FA,category))%>%
  mutate(Value= as.numeric(Value))

Sce_intraR_Twin = intraR_Twin%>%
  select(Item,Value,refCode,is_repeated,category,twinN,novariance_SC)%>%
  filter(category=="SC")%>% ## 4 removed
  select(-c(twinN,novariance_SC,category))%>%
  mutate(Value= as.numeric(Value))

Con_intraR_Twin = intraR_Twin%>%
  select(Item,Value,refCode,is_repeated,category,twinN,novariance_DO)%>%
  filter(category=="DO")%>%
  select(-c(twinN,novariance_DO,category))%>%
  mutate(Value = as.numeric(Value))


#Calculate interRater correlation (reliability score)
Fac_intraR_Twin_output = intraRater(Fac_intraR_Twin)
Sce_intraR_Twin_output = intraRater(Sce_intraR_Twin)
Con_intraR_Twin_output = intraRater(Con_intraR_Twin)

#number of excluded participants
sum(Fac_intraR_Twin_output[[1]]$intraR<.5, na.rm = T)
sum(Sce_intraR_Twin_output[[1]]$intraR<.5, na.rm = T)
# sum(Con_intraR_Twin_output[[1]]$intraR<.5, na.rm = T)

#_Figure:SF1_####
p1_exclusion_fac_val = ggpubr::gghistogram(Fac_intraR_Twin_output[[1]], "intraR") + geom_vline(xintercept = .5) +xlab("intra-rater r faces (val.)")
p1_exclusion_pla_val = ggpubr::gghistogram(Sce_intraR_Twin_output[[1]], "intraR") + geom_vline(xintercept = .5) +xlab("intra-rater r scenes (val.)")
# ggpubr::gghistogram(Con_intraR_Twin_output[[1]], "intraR") + geom_vline(xintercept = .5) +xlab("intra-rater r control (val.)")
save(p1_exclusion_fac_val,p1_exclusion_pla_val,file = sprintf("%s/%s/02_Sutherland_2020/01_sup_FS1_exclusion_Sutherland2020.Rdata",wdOA,wdOA_output))

#rename in order to merge:
names(Fac_intraR_Twin_output[[1]]) = c("refCode","intraR")
names(Sce_intraR_Twin_output[[1]]) = c("refCode","intraR")
names(Con_intraR_Twin_output[[1]]) = c("refCode","intraR")
#save general df
TwinLong = merge(TwinLong,Fac_intraR_Twin_output[[1]][1:2],by = "refCode", all = T)
TwinLong = merge(TwinLong,Sce_intraR_Twin_output[[1]][1:2], by = "refCode", all = T)
TwinLong= merge(TwinLong,Con_intraR_Twin_output[[1]][1:2], by = "refCode", all = T)
TwinLong = TwinLong%>%dplyr::rename(intraR_FA = "intraR.x",intraR_SC = "intraR.y",intraR_DO = "intraR")

####save CTD####
write_csv(TwinLong,sprintf("%s/%s/02_Sutherland_2020/01_CTD_Sutherland2020.csv",wdOA,wdOA_output))

##_Table 1_####
#all
TwinLong%>%
  #apply filtering
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(TwinType,Sex), mean(AgeAtTestingCalculated), range(AgeAtTestingCalculated))

TwinLong%>%
  #apply filtering
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(TwinType,Sex))

#faces
TwinLong%>%
  #apply filtering
  filter(novariance_FA == 0)%>%
  filter(intraR_FA >.5)%>%
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  #tabulate descritives
  summarise(table(TwinType,Sex), mean(AgeAtTestingCalculated), range(AgeAtTestingCalculated))

TwinLong%>%
  #apply filtering
  filter(novariance_FA == 0)%>%
  filter(intraR_FA >.5)%>%
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(TwinType,Sex))

#scenes
TwinLong%>%
  #apply filtering
  filter(novariance_SC == 0)%>%
  filter(intraR_SC >.5)%>%
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  summarise(table(TwinType,Sex))

TwinLong%>%
  #apply filtering
  filter(novariance_SC == 0)%>%
  filter(intraR_SC >.5)%>%
  select(TwinType,Sex, twinN,refCode_2, refCode, AgeAtTestingCalculated)%>%
  #select distinct families
  distinct()%>%
  #create a dummy to indicate whether pairs are complete
  mutate(Pair = ifelse(duplicated(refCode_2)| duplicated(refCode_2, fromLast = TRUE),2,1))%>%
  select(TwinType,Sex,refCode_2,Pair,AgeAtTestingCalculated)%>%
  distinct()%>%
  filter(Pair == 2)%>%
  summarise(table(TwinType,Sex))

####descriptives: Intra-rater####
#compute once more intrarater but after exclusion criteria while cheking that values are similar for first nd second-born twin separately
intraRater(TwinLong%>%
             filter(is_repeated ==1 | is_repeated==2)%>%
             select(Item,Value,refCode,is_repeated,category,twinN,novariance_FA,intraR_FA)%>%
             filter(category=="FA")%>%  
             filter(intraR_FA>.5)%>%
             filter(novariance_FA==0)%>%
             select(-c(twinN,novariance_FA,category,intraR_FA))%>%
             mutate(Value= as.numeric(Value)))
intraRater(TwinLong%>%
             filter(is_repeated ==1 | is_repeated==2)%>%
             select(Item,Value,refCode,is_repeated,category,twinN,novariance_SC,intraR_SC)%>%
             filter(category=="SC")%>%  
             filter(intraR_SC>.5)%>%
             filter(novariance_SC==0)%>%
             select(-c(twinN,novariance_SC,category,intraR_SC))%>%
             mutate(Value= as.numeric(Value)))


####MLM####
#prepare df for MLM (VCA) for twin 1 and 2
#####twin1####
#create a dummy dataframe with repeated images containing NA for non-repeated images
MLM_Image_1 = TwinLong%>%
  filter(is_repeated ==0, twinN == "A")%>%
  dplyr::mutate(is_repeated = 1)%>%
  dplyr::select(Item,Value,refCode,is_repeated,category,TwinType,novariance_FA,novariance_SC, novariance_DO)
MLM_Image_1na = TwinLong%>%
  filter(is_repeated ==0,twinN == "A")%>%
  dplyr::mutate(Value= NA,is_repeated =2 )%>%
  dplyr::select(Item,Value,refCode,is_repeated,category,TwinType,novariance_FA,novariance_SC, novariance_DO)

MLM_Image_Twin1 = rbind(MLM_Image_1, MLM_Image_1na)

#merge the dataframe with the dataframe containing only repeated images
MLM_Twin1 = rbind(intraR_Twin%>%filter(twinN == "A")%>%select(Item,Value,refCode,is_repeated,category,TwinType, novariance_FA,novariance_SC, novariance_DO),MLM_Image_Twin1)
######faces####
#create a dataframe for MLM faces
Fac_MLM_Twin1 = MLM_Twin1%>%
  filter(category=="FA")%>%
  merge(Fac_intraR_Twin_output$intraRater[1:2], by = "refCode", all.x = T)

#sample size per TwinType before exclusion
#table(Fac_MLM_Twin1$TwinType)/400
Fac_MLM_Twin1 = Fac_MLM_Twin1%>%
  filter(intraR>=.5, novariance_FA==0)
#log: filter: removed 5,800 rows (5%), 115,400 rows remaining
#sample size per TwinType after exclusion
# table(Fac_MLM_Twin1$TwinType)/400

#entry as requested in the VCA toolbox
Fac_MLM_Twin1 = Fac_MLM_Twin1%>%
  select(Item,Value,refCode,is_repeated)%>%
  mutate(FamId = substr(refCode,1,nchar(refCode)-2))%>%
  select(Item,Value,FamId,is_repeated) 
######sceneries####
#create a dataframe for MLM faces
Sce_MLM_Twin1 = MLM_Twin1%>%
  filter(category=="SC")%>%
  merge(Sce_intraR_Twin_output$intraRater[1:2], by = "refCode", all.x = T)

#sample size per TwinType before exclusion
# table(Sce_MLM_Twin1$TwinType)/100
Sce_MLM_Twin1 = Sce_MLM_Twin1%>%
  filter(intraR>=.5, novariance_SC==0)
#log: removed 400 rows (1%), 60,200 rows remaining
#sample size per TwinType after exclusion
# table(Sce_MLM_Twin1$TwinType)/100

#entry as requested in the VCA toolbox
Sce_MLM_Twin1 = Sce_MLM_Twin1%>%
  select(Item,Value,refCode,is_repeated)%>%
  mutate(FamId = substr(refCode,1,nchar(refCode)-2))%>%
  select(Item,Value,FamId,is_repeated) 

######control####
#create a dataframe for MLM faces
Con_MLM_Twin1 = MLM_Twin1%>%
  filter(category=="DO")%>%
  merge(Con_intraR_Twin_output$intraRater[1:2], by = "refCode", all.x = T)

#sample size per TwinType before exclusion
# table(Con_MLM_Twin1$TwinType)/100
Con_MLM_Twin1 = Con_MLM_Twin1%>%
  filter(novariance_DO==0) #NOTE, here we don't want to esclude participants with unreliable ratings,as we will model their sum score as a noise control
#log: filter: removed 1,000 rows (1%), 120,200 rows remaining
#sample size per TwinType after exclusion
# table(Con_MLM_Twin1$TwinType)/100

#entry as requested in the VCA toolbox
Con_MLM_Twin1 = Con_MLM_Twin1%>%
  select(Item,Value,refCode,is_repeated)%>%
  mutate(FamId = substr(refCode,1,nchar(refCode)-2))%>%
  select(Item,Value,FamId,is_repeated) 

#####twin2####

#create a dummy dataframe with repeated images containing NA for non-repeated images
MLM_Image_2a = TwinLong%>%
  filter(is_repeated ==0, twinN == "B")%>%
  dplyr::mutate(is_repeated = 1)%>%
  dplyr::select(Item,Value,refCode,is_repeated,category,TwinType,novariance_FA,novariance_SC, novariance_DO)
MLM_Image_2b = TwinLong%>%
  filter(is_repeated ==0,twinN == "B")%>%
  dplyr::mutate(Value= NA,is_repeated =2 )%>%
  dplyr::select(Item,Value,refCode,is_repeated,category,TwinType,novariance_FA,novariance_SC, novariance_DO)

MLM_Image_Twin2 = rbind(MLM_Image_2a, MLM_Image_2b)

#merge the dataframe with the dataframe containing only repeated images
MLM_Twin2 = rbind(intraR_Twin%>%filter(twinN == "B")%>%select(Item,Value,refCode,is_repeated,category, TwinType,novariance_FA,novariance_SC, novariance_DO),MLM_Image_Twin2)

######faces####
#create a dataframe for MLM faces
Fac_MLM_Twin2 = MLM_Twin2%>%
  filter(category=="FA")%>%
  merge(Fac_intraR_Twin_output$intraRater[1:2], by = "refCode", all.x = T)

#sample size per TwinType before exclusion
# table(Fac_MLM_Twin2$TwinType)/400
Fac_MLM_Twin2 = Fac_MLM_Twin2%>%
  filter(intraR>=.5, novariance_FA==0)#18 removed
#filter: removed 4,200 rows (3%), 123,200 rows remaining
#sample size per TwinType after exclusion
table(Fac_MLM_Twin2$TwinType)/400

#entry as requested in the VCA toolbox
Fac_MLM_Twin2 = Fac_MLM_Twin2%>%
  select(Item,Value,refCode,is_repeated)%>%
  mutate(FamId = substr(refCode,1,nchar(refCode)-2))%>%
  select(Item,Value,FamId,is_repeated) 


######sceneries####
#create a dataframe for MLM faces
Sce_MLM_Twin2 = MLM_Twin2%>%
  filter(category=="SC")%>%
  merge(Sce_intraR_Twin_output$intraRater[1:2], by = "refCode", all.x = T)

#sample size per TwinType before exclusion
# table(Sce_MLM_Twin2$TwinType)/100
Sce_MLM_Twin2 = Sce_MLM_Twin2%>%
  filter(intraR>=.50, novariance_SC==0)
#filter: removed 300 rows (<1%), 63,400 rows remaining
#sample size per TwinType after exclusion
# table(Sce_MLM_Twin2$TwinType)/100

#entry as requested in the VCA toolbox
Sce_MLM_Twin2 = Sce_MLM_Twin2%>%
  select(Item,Value,refCode,is_repeated)%>%
  mutate(FamId = substr(refCode,1,nchar(refCode)-2))%>%
  select(Item,Value,FamId,is_repeated) 

######control####
#create a dataframe for MLM faces
Con_MLM_Twin2 = MLM_Twin2%>%
  filter(category=="DO")%>%
  merge(Con_intraR_Twin_output$intraRater[1:2], by = "refCode", all.x = T)

#sample size per TwinType before exclusion
# table(Con_MLM_Twin2$TwinType)/100
Con_MLM_Twin2 = Con_MLM_Twin2%>%
  filter( novariance_DO==0)
#filter: removed 800 rows (1%), 126,600 rows remaining
#sample size per TwinType after exclusion
# table(Con_MLM_Twin2$TwinType)/100

#entry as requested in the VCA toolbox
Con_MLM_Twin2 = Con_MLM_Twin2%>%
  select(Item,Value,refCode,is_repeated)%>%
  mutate(FamId = substr(refCode,1,nchar(refCode)-2))%>%
  select(Item,Value,FamId,is_repeated) 


#####save for MLM####
write_csv(Fac_MLM_Twin1,sprintf("%s/%s/02_Sutherland_2020/01_MLM_Fac_Twin1_Sutherland2020.csv",wdOA,wdOA_output))
write_csv(Sce_MLM_Twin1,sprintf("%s/%s/02_Sutherland_2020/01_MLM_Sce_Twin1_Sutherland2020.csv",wdOA,wdOA_output))
write_csv(Fac_MLM_Twin2,sprintf("%s/%s/02_Sutherland_2020/01_MLM_Fac_Twin2_Sutherland2020.csv",wdOA,wdOA_output))
write_csv(Sce_MLM_Twin2,sprintf("%s/%s/02_Sutherland_2020/01_MLM_Sce_Twin2_Sutherland2020.csv",wdOA,wdOA_output))

####MM2####
#entry for mm2 as required from the MM2 function
Fac_MLM_Twin1$FamId = sprintf("%s-A",Fac_MLM_Twin1$FamId)
Fac_MLM_Twin2$FamId = sprintf("%s-B",Fac_MLM_Twin2$FamId)
Sce_MLM_Twin1$FamId = sprintf("%s-A",Sce_MLM_Twin1$FamId)
Sce_MLM_Twin2$FamId = sprintf("%s-B",Sce_MLM_Twin2$FamId)
Con_MLM_Twin1$FamId = sprintf("%s-A",Con_MLM_Twin1$FamId)
Con_MLM_Twin2$FamId = sprintf("%s-B",Con_MLM_Twin2$FamId)

Fac_MLM_Twin = rbind(Fac_MLM_Twin1, Fac_MLM_Twin2)
Sce_MLM_Twin = rbind(Sce_MLM_Twin1, Sce_MLM_Twin2)
Con_MLM_Twin = rbind(Con_MLM_Twin1, Con_MLM_Twin2)

Fac_MLM_Twin = Fac_MLM_Twin%>%mutate(domain = "FA")
Sce_MLM_Twin = Sce_MLM_Twin%>%mutate(domain = "SC")
Con_MLM_Twin = Con_MLM_Twin%>%mutate(domain = "DO")

MM2_Twin = rbind(Fac_MLM_Twin,Sce_MLM_Twin,Con_MLM_Twin)

#####save for mm2####
write_csv(MM2_Twin,sprintf("%s/%s/02_Sutherland_2020/01_MM2_Sutherland2020.csv",wdOA,wdOA_output))



