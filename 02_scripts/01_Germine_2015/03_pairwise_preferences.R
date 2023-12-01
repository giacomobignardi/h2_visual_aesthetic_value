#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 23-11-2023
#
#
#Description:
#analysis of pairwise preferences (gropued by pair class and domain)
#1. create Unrelated pseudo-Random pairs (UR)
#2. test for pairwise differences across visual domains and class pairs
#3. run post hoc analysis
#Program: 03_pairwise_preferences ------------------------------------------------------------------------------------------------------------------------------

#load packages
library(tidyverse)
library(tidylog)
library(readr)
library(emmeans)

#clean working enviroment 
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdNOA_ImageOutput = "05_Figures"

#load dataFrames:
BioMetric  = read_csv(sprintf("%s/%s/01_Germine_2015/01_CTD_Germine2015.csv", wdOA,wdOA_output))
BioMetric = BioMetric%>%select(-c(FamId_2 )) 
# #load functions:
# source(sprintf("%s/%s/functions/IntraRater.R", wdOA,wdOA_scripts))
####1. PSEUDO-RANDOM PAIRS####
#create a df with ratings per twin pairs in 3 columns. 
#The first 2 colums have ratings from twin 1i and twin 2i TWIN PAIRS
#The third column has the raintings of the twin2i-1 RANDOM PAIRS

#####faces####
Fac_Biometric_Randompairs_r = BioMetric%>%
  #drop no variance
  filter(category =="FA_TOT" & novariance_FA == 0,)%>%
  select(-c(novariance_FA, novariance_SC,novariance_AO,intraR_SC,intraR_AO))%>%
  filter(is_repeated ==1 | is_repeated==2)%>%
  ##average repeated measure
  pivot_wider(names_from = is_repeated, values_from = Value)%>%
  mutate(value = (`1` + `2`)/2)%>%
  select(-c(`1`,`2`))%>%
  #long to wide
  pivot_wider(names_from = SibId, values_from = c("value", "intraR_FA"))%>%
  dplyr::rename(`1` = "value_1",`2`="value_2")%>%
  #drop if one of the pair is missing
  drop_na(`1`)%>%
  drop_na(`2`)%>%
  #drop  if one of the two pair has a whitin score lesse than .5
  filter(intraR_FA_1 >= .5)%>%
  filter(intraR_FA_2 >= .5 )%>%
  select(-c(intraR_FA_1,intraR_FA_2)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex,Age,FamId,Item_2)%>%
  mutate(`3` = lag(`2`,60), age_2 = lag(Age,60), sex_2 = lag(Sex,60)  )

#remove mismatched sex (replace with NA the ratings)
Fac_Biometric_Randompairs_r[which(Fac_Biometric_Randompairs_r$Sex !=Fac_Biometric_Randompairs_r$sex_2),]$`3` = rep(NA,60)
Fac_Biometric_Randompairs_r[which(Fac_Biometric_Randompairs_r$Sex !=Fac_Biometric_Randompairs_r$sex_2),]$age_2 = rep(NA,60)

#create a df with the non repeated images
Fac_Biometric_Randompairs_nr = BioMetric%>%
  filter(category =="FA_TOT" & novariance_FA == 0)%>%
  filter(is_repeated ==0)%>%
  select(-c(is_repeated, novariance_FA, novariance_SC,novariance_AO,intraR_SC,intraR_AO))%>%
  pivot_wider(names_from = SibId, values_from = c("Value", "intraR_FA"))%>%
  dplyr::rename(`1` = "Value_1",`2`="Value_2")%>%
  #drop missing pair
  drop_na(`1`)%>%
  drop_na(`2`)%>%
  #drop if one of the two pair has a whitin score less than .5
  filter(intraR_FA_1 >= .5)%>%
  filter(intraR_FA_2 >= .5)%>%
  select(-c(intraR_FA_1,intraR_FA_2)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex, Age, FamId,Item_2)%>%
  mutate(`3` = lag(`2`,140),age_2 = lag(Age,140), sex_2 = lag(Sex,140))

#remove mismatched sex (replace NA in original ratings)
Fac_Biometric_Randompairs_nr[which(Fac_Biometric_Randompairs_nr$Sex !=Fac_Biometric_Randompairs_nr$sex_2),]$`3` = rep(NA,140)
Fac_Biometric_Randompairs_nr[which(Fac_Biometric_Randompairs_nr$Sex !=Fac_Biometric_Randompairs_nr$sex_2),]$age_2 = rep(NA,140)


#merge repeated and not images in one df
Fac_Biometric_Randompairs = rbind(Fac_Biometric_Randompairs_r,Fac_Biometric_Randompairs_nr)
#print descriptives of differences in age between age groups 
mean(Fac_Biometric_Randompairs$Age - Fac_Biometric_Randompairs$age_2, na.rm = T)
sd(Fac_Biometric_Randompairs$Age - Fac_Biometric_Randompairs$age_2, na.rm = T)
#is the age between the two groups significantly different?
t.test(
  Fac_Biometric_Randompairs%>%distinct(FamId, .keep_all = T)%>%select(Age),
  Fac_Biometric_Randompairs%>%distinct(FamId, .keep_all = T)%>%select(age_2)
)

Fac_Biometric_Randompairs = Fac_Biometric_Randompairs%>%select(FamId,Sex,Age,Zygosity,category,`1`,`2`,`3`)

#####sceneries####
Sce_Biometric_Randompairs_r = BioMetric%>%
  #drop no variance
  filter(category =="SC" & novariance_SC== 0)%>%
  select(-c(novariance_FA, novariance_SC,novariance_AO,intraR_FA,intraR_AO))%>%
  filter(is_repeated ==1 | is_repeated==2)%>%
  ##average repeated measure
  pivot_wider(names_from = is_repeated, values_from = Value)%>%
  mutate(value = (`1` + `2`)/2)%>%
  select(-c(`1`,`2`))%>%
  #long to wide
  pivot_wider(names_from = SibId, values_from = c("value", "intraR_SC"))%>%
  dplyr::rename(`1` = "value_1",`2`="value_2")%>%
  #drop if one of the pair is missing
  drop_na(`1`)%>%
  drop_na(`2`)%>%
  #drop  if one of the two pair has a whitin score lesse than .5
  filter(intraR_SC_1 >= .5)%>%
  filter(intraR_SC_2 >= .5 )%>%
  select(-c(intraR_SC_1,intraR_SC_2)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex,Age,FamId,Item_2)%>%
  mutate(`3` = lag(`2`,15), age_2 = lag(Age,15), sex_2 = lag(Sex,15)  )

#remove mismatched sex (replace with NA the ratings)
Sce_Biometric_Randompairs_r[which(Sce_Biometric_Randompairs_r$Sex !=Sce_Biometric_Randompairs_r$sex_2),]$`3` = rep(NA,15)
Sce_Biometric_Randompairs_r[which(Sce_Biometric_Randompairs_r$Sex !=Sce_Biometric_Randompairs_r$sex_2),]$age_2 = rep(NA,15)

#create a df with the non repeated images
Sce_Biometric_Randompairs_nr = BioMetric%>%
  filter(category =="SC" & novariance_SC== 0)%>%
  filter(is_repeated ==0)%>%
  select(-c(is_repeated, novariance_FA, novariance_SC, novariance_AO, intraR_FA, intraR_AO))%>%
  pivot_wider(names_from = SibId, values_from = c("Value", "intraR_SC"))%>%
  dplyr::rename(`1` = "Value_1",`2`="Value_2")%>%
  #drop missing pair
  drop_na(`1`)%>%
  drop_na(`2`)%>%
  #drop if one of the two pair has a whitin score less than .5
  filter(intraR_SC_1 >= .5)%>%
  filter(intraR_SC_2 >= .5)%>%
  select(-c(intraR_SC_1,intraR_SC_2)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex, Age, FamId,Item_2)%>%
  mutate(`3` = lag(`2`,35),age_2 = lag(35), sex_2 = lag(Sex,35))

#remove mismatched sex (replace NA in original ratings)
Sce_Biometric_Randompairs_nr[which(Sce_Biometric_Randompairs_nr$Sex !=Sce_Biometric_Randompairs_nr$sex_2),]$`3` = rep(NA,35)
Sce_Biometric_Randompairs_nr[which(Sce_Biometric_Randompairs_nr$Sex !=Sce_Biometric_Randompairs_nr$sex_2),]$age_2 = rep(NA,35)

#merge repeated and not images in one df
Sce_Biometric_Randompairs = rbind(Sce_Biometric_Randompairs_r,Sce_Biometric_Randompairs_nr)
#print descriptives of differences in age between age groups 
mean(Sce_Biometric_Randompairs$Age - Sce_Biometric_Randompairs$age_2, na.rm = T)
sd(Sce_Biometric_Randompairs$Age - Sce_Biometric_Randompairs$age_2, na.rm = T)
#is the age between the two groups significantly different?
t.test(
  Sce_Biometric_Randompairs%>%distinct(FamId, .keep_all = T)%>%select(Age),
  Sce_Biometric_Randompairs%>%distinct(FamId, .keep_all = T)%>%select(age_2)
)
Sce_Biometric_Randompairs = Sce_Biometric_Randompairs%>%select(FamId,Sex,Age,Zygosity,category,`1`,`2`,`3`)

#####abstracts####
Obj_Biometric_Randompairs_r = BioMetric%>%
  #drop no variance
  filter(category =="AO" & novariance_AO== 0)%>%
  select(-c(novariance_FA, novariance_SC,novariance_AO,intraR_FA,intraR_SC))%>%
  filter(is_repeated ==1 | is_repeated==2)%>%
  ##average repeated measure
  pivot_wider(names_from = is_repeated, values_from = Value)%>%
  mutate(value = (`1` + `2`)/2)%>%
  select(-c(`1`,`2`))%>%
  #long to wide
  pivot_wider(names_from = SibId, values_from = c("value", "intraR_AO"))%>%
  dplyr::rename(`1` = "value_1",`2`="value_2")%>%
  #drop if one of the pair is missing
  drop_na(`1`)%>%
  drop_na(`2`)%>%
  #drop  if one of the two pair has a whitin score lesse than .5
  filter(intraR_AO_1 >= .5)%>%
  filter(intraR_AO_2 >= .5 )%>%
  select(-c(intraR_AO_1,intraR_AO_2)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex,Age,FamId,Item_2)%>%
  mutate(`3` = lag(`2`,15), age_2 = lag(Age,15), sex_2 = lag(Sex,15)  )

#remove mismatched sex (replace with NA the ratings)
Obj_Biometric_Randompairs_r[which(Obj_Biometric_Randompairs_r$Sex !=Obj_Biometric_Randompairs_r$sex_2),]$`3` = rep(NA,15)
Obj_Biometric_Randompairs_r[which(Obj_Biometric_Randompairs_r$Sex !=Obj_Biometric_Randompairs_r$sex_2),]$age_2 = rep(NA,15)

#create a df with the non repeated images
Obj_Biometric_Randompairs_nr = BioMetric%>%
  filter(category =="AO" & novariance_AO== 0)%>%
  filter(is_repeated ==0)%>%
  select(-c(is_repeated, novariance_FA, novariance_SC, novariance_AO, intraR_FA, intraR_SC))%>%
  pivot_wider(names_from = SibId, values_from = c("Value", "intraR_AO"))%>%
  dplyr::rename(`1` = "Value_1",`2`="Value_2")%>%
  #drop missing pair
  drop_na(`1`)%>%
  drop_na(`2`)%>%
  #drop if one of the two pair has a whitin score less than .5
  filter(intraR_AO_1 >= .5)%>%
  filter(intraR_AO_2 >= .5)%>%
  select(-c(intraR_AO_1,intraR_AO_2)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex, Age, FamId,Item_2)%>%
  mutate(`3` = lag(`2`,35),age_2 = lag(35), sex_2 = lag(Sex,35))

#remove mismatched sex (replace NA in original ratings)
Obj_Biometric_Randompairs_nr[which(Obj_Biometric_Randompairs_nr$Sex !=Obj_Biometric_Randompairs_nr$sex_2),]$`3` = rep(NA,35)
Obj_Biometric_Randompairs_nr[which(Obj_Biometric_Randompairs_nr$Sex !=Obj_Biometric_Randompairs_nr$sex_2),]$age_2 = rep(NA,35)

#merge repeated and not images in one df
Obj_Biometric_Randompairs = rbind(Obj_Biometric_Randompairs_r,Obj_Biometric_Randompairs_nr)
#print descriptives of differences in age between age groups 
mean(Obj_Biometric_Randompairs$Age - Obj_Biometric_Randompairs$age_2, na.rm = T)
sd(Obj_Biometric_Randompairs$Age - Obj_Biometric_Randompairs$age_2, na.rm = T)
#is the age between the two groups significantly different?
t.test(
  Obj_Biometric_Randompairs%>%distinct(FamId, .keep_all = T)%>%select(Age),
  Obj_Biometric_Randompairs%>%distinct(FamId, .keep_all = T)%>%select(age_2)
)


Obj_Biometric_Randompairs = Obj_Biometric_Randompairs%>%select(FamId,Sex,Age,Zygosity,category,`1`,`2`,`3`)

####2. PAIRWISE PREFERENCES####
#calculate correlations within twin pairs and between random pairs
Fac_InterR_twin = Fac_Biometric_Randompairs%>%
  group_by(FamId, Zygosity,category,Sex,Age)%>%
  summarise(r_twin_PP = cor(`1`,`2`), z_twin_PP = psych::fisherz(cor( `1`,`2`)))%>%ungroup()
Fac_InterR_random = Fac_Biometric_Randompairs%>%
  group_by(FamId, Zygosity,category,Sex,Age)%>%
  summarise(r_twin_PP = cor(`1`,`3`), z_twin_PP = psych::fisherz(cor( `1`,`3`)))%>%
  mutate(Zygosity = "UR")%>%ungroup()#unrelated

Sce_InterR_twin = Sce_Biometric_Randompairs%>%
  group_by(FamId, Zygosity,category,Sex,Age)%>%
  summarise(r_twin_PP = cor(`1`,`2`), z_twin_PP = psych::fisherz(cor( `1`,`2`)))%>%ungroup()
Sce_InterR_random = Sce_Biometric_Randompairs%>%
  group_by(FamId, Zygosity,category,Sex,Age)%>%
  summarise(r_twin_PP = cor(`1`,`3`), z_twin_PP = psych::fisherz(cor( `1`,`3`)))%>%
  mutate(Zygosity = "UR")%>%ungroup()#unrelated

Obj_InterR_twin = Obj_Biometric_Randompairs%>%
  group_by(FamId, Zygosity,category,Sex,Age)%>%
  summarise(r_twin_PP = cor(`1`,`2`), z_twin_PP = psych::fisherz(cor( `1`,`2`)))%>%ungroup()
Obj_InterR_random = Obj_Biometric_Randompairs%>%
  group_by(FamId, Zygosity,category,Sex,Age)%>%
  summarise(r_twin_PP = cor(`1`,`3`), z_twin_PP = psych::fisherz(cor( `1`,`3`)))%>%
  mutate(Zygosity = "UR")%>%ungroup()#unrelated

#bind all data together
InterR_random = rbind(
  Fac_InterR_twin,
  Fac_InterR_random,
  Sce_InterR_twin,
  Sce_InterR_random,
  Obj_InterR_twin,
  Obj_InterR_random
                 )

InterR_random = InterR_random%>%
  mutate(category = factor(category, levels = c("AO","SC","FA_TOT")),
         Zygosity = factor(Zygosity, levels = c("UR","DZ","MZ")),
         FamId = ifelse(Zygosity == "UR", as.numeric(FamId) +10000, FamId))%>%
  dplyr::rename(Domain = "category", 
         Pairs = "Zygosity")%>%
  drop_na(r_twin_PP)#drop 6 pairs (mismatched sex and age)

levels(InterR_random$Domain) = c("abstract","scenes","faces")
levels(InterR_random$Pairs) = c("UR","DZ","MZ")

#####ANOVA####
#test for differences between domains (faces, sceneries and objects) and groups (MZ, DZ twin and random pairs)
######with outliers####
anova_interR_PairsXDomain_unballanced_out = car::Anova(lm(z_twin_PP~ + Pairs*Domain, InterR_random),
                                                   contrasts = list(
                                                     Pairs = "contr.sum",
                                                     Domain = "contr.sum"
                                                   ),
                                                   #choose a contrasts setting that sums to zero
                                                   #https://easystats.github.io/effectsize/articles/anovaES.html
                                                   type = 3)
effect_size_anova_out = effectsize::eta_squared(anova_interR_PairsXDomain_unballanced_out, alternative = "two.sided")

######without outliers####
#remove extreme outlier
extreme_outlier = InterR_random %>%
  group_by(Domain, Pairs) %>%
  rstatix::identify_outliers(z_twin_PP)%>%
  filter(is.extreme == "TRUE") 
table(extreme_outlier$Pairs)

#removed 7 extreme outliers
InterR_random_noOut = InterR_random%>%filter(!FamId %in%  c(extreme_outlier$FamId))
#check assumption
#Normality
#inspect violation of normality
ggpubr::ggqqplot(InterR_random_noOut, "z_twin_PP", ggtheme = theme_bw()) +
  facet_grid(Pairs ~ Domain)
#homogeneity of variance (doubt on faces)
InterR_random_noOut %>%
  group_by(Domain) %>%
  rstatix::levene_test(z_twin_PP ~ Pairs)
#not met for faces and scenes
#check directionality
InterR_random_noOut %>%
  group_by(Domain,Pairs) %>%
  summarise(sd(z_twin_PP, na.rm = T))
InterR_random_noOut = InterR_random_noOut%>%mutate(FamId = as.factor(FamId))
# #Homogeneity of covariance met
rstatix::box_m(as.matrix(as.numeric(unlist(InterR_random_noOut%>%drop_na()%>%select("z_twin_PP")))),
      as.matrix(as.numeric(unlist(InterR_random_noOut%>%drop_na()%>%select("Pairs")))))

#ANOVA model
PairsXDomain_mod = lm(z_twin_PP~ + Pairs*Domain, InterR_random_noOut)
#take into account unballanced desing (type III)
anova_interR_PairsXDomain_unballanced = car::Anova(PairsXDomain_mod,
                                                   contrasts = list(
                                                     Pairs = "contr.sum",
                                                     Domain = "contr.sum"
                                                   ),
                                                   #choose a contrasts setting that sums to zero
                                                   #https://easystats.github.io/effectsize/articles/anovaES.html
                                                   type = 3)

#ANOVA effect sizes
effect_size_anova = effectsize::eta_squared(anova_interR_PairsXDomain_unballanced, alternative = "two.sided")

#the follwing analysis are adapted from https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html and  https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/ 
####3. POST HOC####
#note that this section is pretty long but mainly because it is tiding the outputs from the emmeans comparison
#we test for 3 family-wise comparison (bonferroni) :
#1. pairwise difference in pairwise agreement across domain
#2. pairwise difference in pairwise agreement across pair classes
#3. pairwise difference in pairwise agreement across pair classes within domains
emmDomain = emmeans(PairsXDomain_mod, specs = pairwise ~ Domain)
#Estimates
emmDomain_estimates = emmDomain%>%
  as.data.frame()%>%
  rename(emmean_z = "emmean",
         lower.CL_z = "lower.CL",
         upper.CL_z = "upper.CL")%>%
  mutate(emmean_r = psych::fisherz2r(emmean_z), 
         lower.CL_r =  psych::fisherz2r(lower.CL_z),
         upper.CL_r = psych::fisherz2r(upper.CL_z))

emmDomain_estimates_csv = emmDomain_estimates%>%
  filter(contrast == ".")%>%
  select(-df)%>%
  mutate(emmean_z = round(emmean_z, 3),
         lower.CL_z =  round(lower.CL_z, 3),
         upper.CL_z = round(upper.CL_z, 3),
         emmean_r = round(emmean_r, 3),
         lower.CL_r = round(lower.CL_r, 3),
         upper.CL_r = round(upper.CL_r, 3)
  )%>%
  select(Domain,
         emmean_z,
         lower.CL_z ,
         upper.CL_z ,
         emmean_r ,
         lower.CL_r,
         upper.CL_r)

#Statistics
emmDomain_output = emmDomain$contrasts%>%
  rbind()%>% #control for bonferroni
  summary(infer = TRUE) %>%
  as.data.frame()%>%# calculate CI for changes
  rename(estimate_z = "estimate",
         lower.CL_z = "lower.CL",
         upper.CL_z = "upper.CL")%>%
  mutate(estimate_r = psych::fisherz2r(estimate_z), 
         lower.CL_r =  psych::fisherz2r(lower.CL_z),
         upper.CL_r = psych::fisherz2r(upper.CL_z))

#Effect sizes
emmDomain_effect = eff_size(emmDomain, 
                            sigma = sigma(PairsXDomain_mod), 
                            edf = df.residual(PairsXDomain_mod),
                            method = "identity")%>%
  as.data.frame()%>%
  rename(effect.size = "effect.size",
         e.s_lower.CL = "lower.CL",
         e.s_upper.CL= "upper.CL",
  )%>%
  mutate(contrast = substr(contrast, 2, nchar(as.character(contrast))-1))

#Prepare CSV
emmDomain_stats = merge(emmDomain_estimates%>%
                          select(-c(SE, df,lower.CL_z, upper.CL_z,lower.CL_r,  upper.CL_r)),emmDomain_output,by = "contrast")
emmDomain_csv = merge(emmDomain_stats,emmDomain_effect%>%select(-c(SE, df)),by = "contrast")%>%
  mutate(
    dif_emmean_z = round(emmean_z,3),
    dif_emmean_r  = round(emmean_r,3),
    lower.CL_z = round(lower.CL_z,3),
    upper.CL_z = round(upper.CL_z,3),
    lower.CL_r = round(lower.CL_r,3),  
    upper.CL_r = round(upper.CL_r,3),
    df = df,
    t.ratio =  round(t.ratio,3),
    effect.size = round(effect.size,3),
    e.s_lower.CL = round(e.s_lower.CL,3),
    e.s_upper.CL = round(e.s_upper.CL,3),
  )%>%
  select(
    contrast,
    dif_emmean_z,
    dif_emmean_r,
    lower.CL_z, 
    upper.CL_z, 
    lower.CL_r, 
    upper.CL_r, 
    df, 
    t.ratio, 
    p.value,  
    effect.size,
    e.s_lower.CL,
    e.s_upper.CL,
  )


####_Pair comparison####
emmPairs = emmeans(PairsXDomain_mod, specs = pairwise ~ Pairs, adjust = "bonferroni")

#Estimates
emmPairs_estimates = emmPairs%>%
  as.data.frame()%>%
  rename(emmean_z = "emmean",
         lower.CL_z = "lower.CL",
         upper.CL_z = "upper.CL")%>%
  mutate(emmean_r = psych::fisherz2r(emmean_z), 
         lower.CL_r =  psych::fisherz2r(lower.CL_z),
         upper.CL_r = psych::fisherz2r(upper.CL_z))

emmPairs_estimates_csv = emmPairs_estimates%>%
  filter(contrast == ".")%>%
  select(-df)%>%
  mutate(emmean_z = round(emmean_z, 3),
         lower.CL_z =  round(lower.CL_z, 3),
         upper.CL_z = round(upper.CL_z, 3),
         emmean_r = round(emmean_r, 3),
         lower.CL_r = round(lower.CL_r, 3),
         upper.CL_r = round(upper.CL_r, 3)
  )%>%
  select(Pairs,
         emmean_z,
         lower.CL_z ,
         upper.CL_z ,
         emmean_r ,
         lower.CL_r,
         upper.CL_r)

#Statistics
emmPairs_output = emmPairs$contrasts%>%
  rbind()%>% #control for bonferroni
  summary(infer = TRUE) %>%
  as.data.frame()%>%# calculate CI for changes
  rename(estimate_z = "estimate",
         lower.CL_z = "lower.CL",
         upper.CL_z = "upper.CL")%>%
  mutate(estimate_r = psych::fisherz2r(estimate_z), 
         lower.CL_r =  psych::fisherz2r(lower.CL_z),
         upper.CL_r = psych::fisherz2r(upper.CL_z))

#Effect sizeds
emmPairs_effect = eff_size(emmPairs, 
                           sigma = sigma(PairsXDomain_mod), 
                           edf = df.residual(PairsXDomain_mod),
                           method = "identity")%>%
  as.data.frame()%>%
  rename(effect.size = "effect.size",
         e.s_lower.CL = "lower.CL",
         e.s_upper.CL= "upper.CL",
  )%>%
  mutate(contrast = substr(contrast, 2, nchar(as.character(contrast))-1))

#Prepare CSV
emmPairs_stats = merge(emmPairs_estimates%>%
                         select(-c(SE, df,lower.CL_z, upper.CL_z,lower.CL_r,  upper.CL_r)),emmPairs_output,by = "contrast")
emmPairs_csv = merge(emmPairs_stats,emmPairs_effect%>%select(-c(SE, df)),by = "contrast")%>%
  mutate(
    dif_emmean_z = round(emmean_z,3),
    dif_emmean_r  = round(emmean_r,3),
    lower.CL_z = round(lower.CL_z,3),
    upper.CL_z = round(upper.CL_z,3),
    lower.CL_r = round(lower.CL_r,3),  
    upper.CL_r = round(upper.CL_r,3),
    df = df,
    t.ratio =  round(t.ratio,3),
    effect.size = round(effect.size,3),
    e.s_lower.CL = round(e.s_lower.CL,3),
    e.s_upper.CL = round(e.s_upper.CL,3),
  )%>%
  select(
    contrast,
    dif_emmean_z,
    dif_emmean_r,
    lower.CL_z, 
    upper.CL_z, 
    lower.CL_r, 
    upper.CL_r, 
    df, 
    t.ratio, 
    p.value,  
    effect.size,
    e.s_lower.CL,
    e.s_upper.CL,
  )

####_Pair|domain comparison####
emmPairs2Domain = emmeans(PairsXDomain_mod, specs = pairwise ~ Pairs|Domain, adjust = "bonferroni")

#Estimates
emmPairs2Domain_estimates = emmPairs2Domain%>%
  as.data.frame()%>%
  rename(emmean_z = "emmean",
         lower.CL_z = "lower.CL",
         upper.CL_z = "upper.CL")%>%
  mutate(emmean_r = psych::fisherz2r(emmean_z), 
         lower.CL_r =  psych::fisherz2r(lower.CL_z),
         upper.CL_r = psych::fisherz2r(upper.CL_z))

emmPairs2Domain_estimates_csv = emmPairs2Domain_estimates%>%
  filter(contrast == ".")%>%
  select(-df)%>%
  mutate(emmean_z = round(emmean_z, 3),
         lower.CL_z =  round(lower.CL_z, 3),
         upper.CL_z = round(upper.CL_z, 3),
         emmean_r = round(emmean_r, 3),
         lower.CL_r = round(lower.CL_r, 3),
         upper.CL_r = round(upper.CL_r, 3)
  )%>%
  select(Pairs,
         Domain,
         emmean_z,
         lower.CL_z ,
         upper.CL_z ,
         emmean_r ,
         lower.CL_r,
         upper.CL_r)

#Statistics
emmPairs2Domain_output = emmPairs2Domain$contrasts%>%
  rbind()%>% #control for bonferroni
  summary(infer = TRUE) %>% 
  as.data.frame()%>%# calculate CI for changes
  rename(estimate_z = "estimate",
         lower.CL_z = "lower.CL",
         upper.CL_z = "upper.CL")%>%
  mutate(estimate_r = psych::fisherz2r(estimate_z), 
         lower.CL_r =  psych::fisherz2r(lower.CL_z),
         upper.CL_r = psych::fisherz2r(upper.CL_z))

#Effect sizeds
emmPairs2Domain_effect = eff_size(emmPairs2Domain, 
                                  sigma = sigma(PairsXDomain_mod), 
                                  edf = df.residual(PairsXDomain_mod),
                                  method = "identity")%>%
  as.data.frame()%>%
  rename(effect.size = "effect.size",
         e.s_lower.CL = "lower.CL",
         e.s_upper.CL= "upper.CL",
  )%>%
  mutate(contrast = substr(contrast, 2, nchar(as.character(contrast))-1))

#REPORT:aesthetic agreement####
save(anova_interR_PairsXDomain_unballanced,
     anova_interR_PairsXDomain_unballanced_out,
     effect_size_anova,
     effect_size_anova_out,
     emmPairs_estimates_csv,
     emmDomain_estimates_csv,
     emmPairs2Domain_output,
     emmPairs2Domain_effect,file = sprintf("%s/%s/01_Germine_2015/03_report_aesthetic_agreement.Rdata",wdOA,wdOA_output))

#Prepare CSV
emmPairs2Domain_stats = merge(emmPairs2Domain_estimates%>%
                                select(-c(SE, df,lower.CL_z, upper.CL_z,lower.CL_r,  upper.CL_r)),emmPairs2Domain_output,by = c("contrast", "Domain"))
emmPairs2Domain_csv = 
  merge(emmPairs2Domain_stats,emmPairs2Domain_effect%>%select(-c(SE, df)),by = c("contrast","Domain"))%>%
  mutate(
    dif_emmean_z = round(emmean_z,3),
    dif_emmean_r  = round(emmean_r,3),
    lower.CL_z = round(lower.CL_z,3),
    upper.CL_z = round(upper.CL_z,3),
    lower.CL_r = round(lower.CL_r,3),  
    upper.CL_r = round(upper.CL_r,3),
    df = df,
    t.ratio =  round(t.ratio,3),
    effect.size = round(effect.size,3),
    e.s_lower.CL = round(e.s_lower.CL,3),
    e.s_upper.CL = round(e.s_upper.CL,3)
  )%>%
  select(
    contrast,
    Domain,
    dif_emmean_z,
    dif_emmean_r,
    lower.CL_z, 
    upper.CL_z, 
    lower.CL_r, 
    upper.CL_r, 
    df, 
    t.ratio, 
    p.value,  
    effect.size,
    e.s_lower.CL,
    e.s_upper.CL
  )

#####Supplementary pairwise csv####
write_csv(emmDomain_estimates_csv,sprintf("%s/%s/01_Germine_2015/03_sup_pairwise_domain_estimates_Germine2015.csv",wdOA,wdOA_output))
write_csv(emmDomain_csv,sprintf("%s/%s/01_Germine_2015/03_sup_pairwise_domain_statistics_Germine2015.csv",wdOA,wdOA_output))
write_csv(emmPairs_estimates_csv,sprintf("%s/%s/01_Germine_2015/03_sup_pairwise_pairs_estimates_Germine2015.csv",wdOA,wdOA_output))
write_csv(emmPairs_csv,sprintf("%s/%s/01_Germine_2015/03_sup_pairwise_pairs_statistics_Germine2015.csv",wdOA,wdOA_output))
write_csv(emmPairs2Domain_estimates_csv,sprintf("%s/%s/01_Germine_2015/03_sup_pairwise_pairs2domain_estimates_Germine2015.csv",wdOA,wdOA_output))
write_csv(emmPairs2Domain_csv,sprintf("%s/%s/01_Germine_2015/03_sup_pairwise_pairs2domain_statistics_Germine2015.csv",wdOA,wdOA_output))

#save pariswise agreeemnt
write_csv(InterR_random_noOut,sprintf("%s/%s/01_Germine_2015/03_pairwise_agreement_Germine2015.csv",wdOA,wdOA_output))

# ####Figure 2: before editorial####
# PairsXDomain_mod_r = lm(r_twin_PP ~ +Pairs * Domain, data = InterR_random_noOut)
# #final plot (r transformed)
# p_pairwise_r = emmip(PairsXDomain_mod_r, Pairs ~ Domain, CIs = TRUE) +
#   theme_bw(base_size = 14) + 
#   ylim(c(0.2,0.8))+
#   scale_color_viridis_d() + #note that final color has been changed after Review 01
#   labs(y = "pairwsie agreement",
#        x = "Visual domains")
# #save for later plotting
# save(p_pairwise_r,file =sprintf("%s/%s/01_Germine_2015/03_Fig_F2a_pairwise_agreement.Rdata",wdOA,wdOA_output))

####Figure 2: after editorial request####
#Sample mean and CI
InterR_random_noOut_mean = InterR_random_noOut %>%
  group_by(Domain,Pairs) %>% 
  summarise(mean = mean(z_twin_PP, na.rm = TRUE),
            sd = sd(z_twin_PP, na.rm = TRUE),
            n = n()) %>%
  #calculate lower and upped boundaries of 95% CI
  mutate(se = sd/sqrt(n),
         lower.ci = psych::fisherz2r(mean - qt(1 - (0.05 / 2), df = n - 1) * se),
         upper.ci = psych::fisherz2r(mean + qt(1 - (0.05 / 2), df = n - 1) * se),
         r_twin_PP = psych::fisherz2r(mean))%>% 
  select(Domain, Pairs,r_twin_PP, lower.ci,upper.ci) 

#visualzie individual datapoints and CI
set.seed(42)
p_pairwise_r_revised = ggplot(InterR_random_noOut_mean) + 
  geom_jitter(data = InterR_random_noOut,aes(Domain,r_twin_PP, color = Pairs),
              alpha = .05,
              width = 0.1,
              size =.50) +
  geom_line(aes(Domain, r_twin_PP, color = Pairs, group = Pairs))+
  geom_point(aes(Domain, r_twin_PP, color = Pairs), size = 1.5)+
  geom_errorbar(aes(x =Domain, ymin = lower.ci, ymax = upper.ci, color = Pairs), width = .1)+
  scale_color_viridis_d() + #note that final color has been changed after Review 01
  scale_y_continuous(lim= c(-0.5,1), breaks = seq(-0.5,1,by = .25))+
  labs(y = "pairwise agreement",
       x = "Visual domains")+
  theme_light()

p_pairwise_r_revised_zoom = ggplot(InterR_random_noOut_mean) + 
  geom_line(aes(Domain, r_twin_PP, color = Pairs, group = Pairs))+
  geom_point(aes(Domain, r_twin_PP, color = Pairs), size = 1.5)+
  geom_errorbar(aes(x =Domain, ymin = lower.ci, ymax = upper.ci, color = Pairs), width = .1)+
  scale_color_viridis_d() + #note that final color has been changed after Review 01
  scale_y_continuous(lim= c(.25,.85), breaks = seq(.25,.85,by = .1))+
  labs(y = "pairwise agreement",
       x = "Visual domains")+
  theme_light()

#save for later plotting
save(p_pairwise_r_revised,p_pairwise_r_revised_zoom,file =sprintf("%s/%s/01_Germine_2015/03_Fig_F2a_pairwise_agreement_revised.Rdata",wdOA,wdOA_output))

##Supplementary File 5 (Editorial Request)####
#save source data
write_csv(InterR_random_noOut_mean %>% 
            mutate(pairwise.r = round(r_twin_PP,3),lower.ci = round(lower.ci,3),upper.ci = round(upper.ci,3))%>% 
            select(Domain,Pairs, pairwise.r,lower.ci,upper.ci),
          sprintf("%s/%s/01_Germine_2015/03_Fig_F2a_summary_pairwise_agreement_revised_sourceData.csv",wdOA,wdOA_output))




#Supplementary Fig 3####
#intra rater dist
p1 = ggplot(InterR_random_noOut, aes( y = r_twin_PP, fill = Pairs)) +
  geom_density(alpha = .5) + 
  scale_fill_viridis_d() +
  facet_wrap(~Domain) +
  labs(y = "inter-rater r") + 
  theme_classic(base_size = 14)
#intra rater z transform dist
p2 = ggplot(InterR_random_noOut, aes( y = z_twin_PP, fill = Pairs)) +
  geom_density(alpha = .5) + 
  scale_fill_viridis_d() +
  labs(y = "inter-rater z") + 
  facet_wrap(~Domain) +
  theme_classic(base_size = 14)
#FIg2a for z values
p_pairwise_z = emmip(PairsXDomain_mod, Pairs ~ Domain) +
  theme_bw(base_size = 16) + 
  scale_color_viridis_d() +
  labs(y = "pairwsie agreement",
       x = "Visual domains")
#marginal mean effects domains
p3 = plot(emmDomain, comparisons = TRUE) + 
  theme_bw(base_size = 14) +
  labs(x = "Marginal means inter-rater Z", y = "Domains")
p4 = plot(emmPairs, comparisons = TRUE) + 
  theme_bw(base_size = 14) +
  labs(x = "Marginal means inter-rater Z", y = "Domains")
#marginal mean effects pairs|domains
p5 = plot(emmPairs2Domain, comparisons = TRUE) + 
  theme_bw(base_size = 14) +
  labs(x = "Marginal means inter-rater Z", y = "Pairs")

#save for later plotting
save(p1,p2,p3,p4,p5,p_pairwise_z, file =sprintf("%s/%s/01_Germine_2015/03_sup_FS6_pairwise_agreement.Rdata",wdOA,wdOA_output))
