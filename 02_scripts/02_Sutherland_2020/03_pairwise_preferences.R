#Author: Giacomo Bignardi
#Adapted from: NA
#Date: 28-04-2021
#Last modified: 18-09-2023
#
#
#Program: pairwise_preferences ------------------------------------------------------------------------------------------------------------------------------
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
BioMetric  = read_csv(sprintf("%s/%s/02_Sutherland_2020/01_CTD_Sutherland2020.csv", wdOA,wdOA_output))
####1. PSEUDO-RANDOM PAIRS####
#create a df with ratings per twin pairs in 3 columns. 
#The first 2 colums have ratings from twin 1i and twin 2i TWIN PAIRS
#The third column has the raintings of the twin2i-1 RANDOM PAIRS

#####faces####
Fac_Biometric_Randompairs_r = BioMetric%>%
  #drop no variance
  filter(category =="FA" & novariance_FA == 0,)%>%
  mutate(refCode = substr(refCode,1, nchar(refCode)-2))%>%
  select(-c(ZygosityStatus,novariance_FA, novariance_SC,novariance_DO,intraR_SC,intraR_DO))%>%
  filter(is_repeated ==1 | is_repeated==2)%>%
  ##averAgeAtTestingCalculated repeated measure
  pivot_wider(names_from = is_repeated, values_from = Value)%>%
  mutate(value = (`1` + `2`)/2)%>%
  select(-c(`1`,`2`))%>%
  #long to wide
  pivot_wider(names_from = twinN, values_from = c("value", "intraR_FA"))%>%
  dplyr::rename(A = "value_A",B="value_B")%>%
  #drop if one of the pair is missing
  drop_na(A)%>%
  drop_na(B)%>%
  #drop  if one of the two pair has a whitin score lesse than .5
  filter(intraR_FA_A >= .5)%>%
  filter(intraR_FA_B >= .5 )%>%
  select(-c(intraR_FA_A,intraR_FA_B)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex,AgeAtTestingCalculated,refCode,Item)%>%
  mutate(C = lag(B,50), AgeAtTestingCalculated_2 = lag(AgeAtTestingCalculated,50), sex_2 = lag(Sex,50)  )

#remove mismatched sex (replace with NA the ratings)
Fac_Biometric_Randompairs_r[which(Fac_Biometric_Randompairs_r$Sex !=Fac_Biometric_Randompairs_r$sex_2),]$C = rep(NA,50)
Fac_Biometric_Randompairs_r[which(Fac_Biometric_Randompairs_r$Sex !=Fac_Biometric_Randompairs_r$sex_2),]$AgeAtTestingCalculated_2 = rep(NA,50)

#create a df with the non repeated imAgeAtTestingCalculateds
Fac_Biometric_Randompairs_nr = BioMetric%>%
  filter(category =="FA" & novariance_FA == 0)%>%
  mutate(refCode = substr(refCode,1, nchar(refCode)-2))%>%
  filter(is_repeated ==0)%>%
  select(-c(ZygosityStatus, is_repeated, novariance_FA, novariance_SC,novariance_DO,intraR_SC,intraR_DO))%>%
  pivot_wider(names_from = twinN, values_from = c("Value", "intraR_FA"))%>%
  dplyr::rename(A = "Value_A",B="Value_B")%>%
  #drop missing pair
  drop_na(A)%>%
  drop_na(B)%>%
  #drop if one of the two pair has a whitin score less than .5
  filter(intraR_FA_A >= .5)%>%
  filter(intraR_FA_B >= .5)%>%
  select(-c(intraR_FA_A,intraR_FA_B)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex, AgeAtTestingCalculated, refCode,Item)%>%
  mutate(C = lag(B,50),AgeAtTestingCalculated_2 = lag(AgeAtTestingCalculated,50), sex_2 = lag(Sex,50))

#remove mismatched sex (replace NA in original ratings)
Fac_Biometric_Randompairs_nr[which(Fac_Biometric_Randompairs_nr$Sex !=Fac_Biometric_Randompairs_nr$sex_2),]$C = rep(NA,50)
Fac_Biometric_Randompairs_nr[which(Fac_Biometric_Randompairs_nr$Sex !=Fac_Biometric_Randompairs_nr$sex_2),]$AgeAtTestingCalculated_2 = rep(NA,50)


#merge repeated and not imAgeAtTestingCalculateds in one df
Fac_Biometric_Randompairs = rbind(Fac_Biometric_Randompairs_r,Fac_Biometric_Randompairs_nr)
#print descriptives of differences in AgeAtTestingCalculated between AgeAtTestingCalculated groups 
mean(Fac_Biometric_Randompairs$AgeAtTestingCalculated - Fac_Biometric_Randompairs$AgeAtTestingCalculated_2, na.rm = T)
sd(Fac_Biometric_Randompairs$AgeAtTestingCalculated - Fac_Biometric_Randompairs$AgeAtTestingCalculated_2, na.rm = T)
#is the AgeAtTestingCalculated between the two groups significantly different?
t.test(Fac_Biometric_Randompairs$AgeAtTestingCalculated,Fac_Biometric_Randompairs$AgeAtTestingCalculated_2)

Fac_Biometric_Randompairs = Fac_Biometric_Randompairs%>%select(refCode,Sex,AgeAtTestingCalculated,TwinType,category,Item,A,B,C)

#####sceneries####
Sce_Biometric_Randompairs_r = BioMetric%>%
  #drop no variance
  filter(category =="SC" & novariance_SC == 0,)%>%
  mutate(refCode = substr(refCode,1, nchar(refCode)-2))%>%
  select(-c(ZygosityStatus,novariance_SC, novariance_FA,novariance_DO,intraR_FA,intraR_DO))%>%
  filter(is_repeated ==1 | is_repeated==2)%>%
  ##averAgeAtTestingCalculated repeated measure
  pivot_wider(names_from = is_repeated, values_from = Value)%>%
  mutate(value = (`1` + `2`)/2)%>%
  select(-c(`1`,`2`))%>%
  #long to wide
  pivot_wider(names_from = twinN, values_from = c("value", "intraR_SC"))%>%
  dplyr::rename(A = "value_A",B="value_B")%>%
  #drop if one of the pair is missing
  drop_na(A)%>%
  drop_na(B)%>%
  #drop  if one of the two pair has a whitin score lesse than .5
  filter(intraR_SC_A >= .5)%>%
  filter(intraR_SC_B >= .5 )%>%
  select(-c(intraR_SC_A,intraR_SC_B)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex,AgeAtTestingCalculated,refCode,Item)%>%
  mutate(C = lag(B,24), AgeAtTestingCalculated_2 = lag(AgeAtTestingCalculated,24), sex_2 = lag(Sex,24)  )

#remove mismatched sex (replace with NA the ratings)
Sce_Biometric_Randompairs_r[which(Sce_Biometric_Randompairs_r$Sex !=Sce_Biometric_Randompairs_r$sex_2),]$C = rep(NA,24)
Sce_Biometric_Randompairs_r[which(Sce_Biometric_Randompairs_r$Sex !=Sce_Biometric_Randompairs_r$sex_2),]$AgeAtTestingCalculated_2 = rep(NA,24)

#create a df with the non repeated imAgeAtTestingCalculateds
Sce_Biometric_Randompairs_nr = BioMetric%>%
  filter(category =="SC" & novariance_SC == 0)%>%
  mutate(refCode = substr(refCode,1, nchar(refCode)-2))%>%
  filter(is_repeated ==0)%>%
  select(-c(ZygosityStatus, is_repeated, novariance_FA, novariance_SC,novariance_DO,intraR_FA,intraR_DO))%>%
  pivot_wider(names_from = twinN, values_from = c("Value", "intraR_SC"))%>%
  dplyr::rename(A = "Value_A",B="Value_B")%>%
  #drop missing pair
  drop_na(A)%>%
  drop_na(B)%>%
  #drop if one of the two pair has a whitin score less than .5
  filter(intraR_SC_A >= .5)%>%
  filter(intraR_SC_B >= .5)%>%
  select(-c(intraR_SC_A,intraR_SC_B)) %>%
  #shift ratings of twin 2 (i) to twin 2 (i+1). This will create a "random" pair (unfamiliar cluster)
  arrange(Sex, AgeAtTestingCalculated, refCode,Item)%>%
  mutate(C = lag(B,26),AgeAtTestingCalculated_2 = lag(AgeAtTestingCalculated,26), sex_2 = lag(Sex,26))

#remove mismatched sex (replace NA in original ratings)
Sce_Biometric_Randompairs_nr[which(Sce_Biometric_Randompairs_nr$Sex !=Sce_Biometric_Randompairs_nr$sex_2),]$C = rep(NA,26)
Sce_Biometric_Randompairs_nr[which(Sce_Biometric_Randompairs_nr$Sex !=Sce_Biometric_Randompairs_nr$sex_2),]$AgeAtTestingCalculated_2 = rep(NA,26)


#merge repeated and not imAgeAtTestingCalculateds in one df
Sce_Biometric_Randompairs = rbind(Sce_Biometric_Randompairs_r,Sce_Biometric_Randompairs_nr)
#print descriptives of differences in AgeAtTestingCalculated between AgeAtTestingCalculated groups 
mean(Sce_Biometric_Randompairs$AgeAtTestingCalculated - Sce_Biometric_Randompairs$AgeAtTestingCalculated_2, na.rm = T)
sd(Sce_Biometric_Randompairs$AgeAtTestingCalculated - Sce_Biometric_Randompairs$AgeAtTestingCalculated_2, na.rm = T)
#is the AgeAtTestingCalculated between the two groups significantly different?
t.test(Sce_Biometric_Randompairs$AgeAtTestingCalculated,Sce_Biometric_Randompairs$AgeAtTestingCalculated_2)

Sce_Biometric_Randompairs = Sce_Biometric_Randompairs%>%select(refCode,Sex,AgeAtTestingCalculated,TwinType,category,Item,A,B,C)

####2. PAIRWISE PREFERENCES####
#calculate correlations between twin pairs and between random pairs
Fac_InterR_twin = Fac_Biometric_Randompairs%>%
  group_by(refCode, TwinType,category)%>%
  summarise(r_twin_PP = cor(A,B), z_twin_PP = psych::fisherz(cor( A,B)))
Fac_InterR_random = Fac_Biometric_Randompairs%>%
  group_by(refCode, TwinType,category)%>%
  summarise(r_twin_PP = cor(A,C), z_twin_PP = psych::fisherz(cor( A,C)))%>%
  mutate(TwinType = "UR")#unrelated

Sce_InterR_twin = Sce_Biometric_Randompairs%>%
  group_by(refCode, TwinType,category)%>%
  summarise(r_twin_PP = cor(A,B), z_twin_PP = psych::fisherz(cor( A,B)))
Sce_InterR_random = Sce_Biometric_Randompairs%>%
  group_by(refCode, TwinType,category)%>%
  summarise(r_twin_PP = cor(A,C), z_twin_PP = psych::fisherz(cor( A,C)))%>%
  mutate(TwinType = "UR")#unrelated


InterR_random = rbind(
  Fac_InterR_twin,
  Fac_InterR_random,
  Sce_InterR_twin,
  Sce_InterR_random
                 )

InterR_random = InterR_random%>%
  mutate(category = factor(category, levels = c("SC", "FA")),
         TwinType = factor(TwinType, levels = c("UR","DZ","MZ")),
         refCode = ifelse(TwinType == "UR", paste0(refCode,"_UR"),  paste0(refCode,"_TW")))%>%
  dplyr::rename(Domain = "category", 
         Pairs = "TwinType")%>%
  ungroup(refCode,Pairs)%>%
  drop_na(r_twin_PP)#drop 6 pairs (mismatched sex and AgeAtTestingCalculated)

levels(InterR_random$Domain) = c("scenes val.", "faces val.")
levels(InterR_random$Pairs) = c("UR","DZ","MZ")

#####ANOVA####
#test for differences between domains (faces, sceneries and objects) and groups (MZ, DZ twin and random pairs)
#test if anova is the same without outliers

#####with outliers####
anova_interR_PairsXDomain_unballanced_val_out = car::Anova(lm(z_twin_PP~ + Pairs*Domain, InterR_random),
                                                       contrasts = list(
                                                         Pairs = "contr.sum",
                                                         Domain = "contr.sum"
                                                       ),
                                                       #choose a contrasts setting that sums to zero
                                                       #https://easystats.github.io/effectsize/articles/anovaES.html
                                                       type = 3)
effect_size_anova_val_out = effectsize::eta_squared(anova_interR_PairsXDomain_unballanced_val_out, alternative = "two.sided")

#remove extreme outlier
extreme_outlier = InterR_random %>%
  group_by(Domain, Pairs) %>%
  rstatix::identify_outliers(z_twin_PP)%>%
  filter(is.extreme == "TRUE") 

table(extreme_outlier$Pairs)
#removed 3 extreme outliers
InterR_random_noOut = InterR_random%>%filter(!refCode %in%  c(extreme_outlier$refCode))

####_Assumption: Normality####
#inspect violation of normality
ggpubr::ggqqplot(InterR_random_noOut, "z_twin_PP", ggtheme = theme_bw()) +
  facet_grid(Pairs ~ Domain)

#homogeneity of variance (doubt on faces)
InterR_random_noOut %>%
  group_by(Domain) %>%
  rstatix::levene_test(z_twin_PP ~ Pairs)


InterR_random_noOut = InterR_random_noOut%>%mutate(refCode = as.factor(refCode))
# #Homogeneity of covariance met
rstatix::box_m(as.matrix(as.numeric(unlist(InterR_random_noOut%>%drop_na()%>%select("z_twin_PP")))),
      as.matrix(as.numeric(unlist(InterR_random_noOut%>%drop_na()%>%select("Pairs")))))


PairsXDomain_mod = lm(z_twin_PP~ Pairs * Domain, InterR_random_noOut)
#take into account unballanced desing (type III)
anova_interR_PairsXDomain_unballanced_val = car::Anova(PairsXDomain_mod,
                                                        contrasts = list(
                                                          Pairs = "contr.sum",
                                                          Domain = "contr.sum"
                                                        ),
                                                        #choose a contrasts setting that sums to zero
                                                        #https://easystats.github.io/effectsize/articles/anovaES.html
                                                        type = 3)
effect_size_anova_val = effectsize::eta_squared(anova_interR_PairsXDomain_unballanced_val, alternative = "two.sided")

#the follwing analysis are adapted from https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html and  https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/ 
####3. POST HOC####
#note that this section is pretty long but mainly because it is tyding the outputs from the emmeans comparison
#we test for 3 family wise comparison and bonferroni :
#1. pairwise difference in pairwise agreement across domain
#2. pairwise difference in pairwise agreement across pair classes
#3. pairwise difference in pairwise agreement across pair classes within domainsemmDomain = emmeans(PairsXDomain_mod, specs = pairwise ~ Domain, adjust = "Bonferroni")
#Estimates
emmDomain = emmeans(PairsXDomain_mod, specs = pairwise ~ Domain)
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

#Effect sizeds
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

emmPairs_estimates_csv_val = emmPairs_estimates%>%
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
emmPairs2Domain = emmeans(PairsXDomain_mod, specs = pairwise ~ Pairs|Domain,adjust = "bonferroni")

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
emmPairs2Domain_output_val = emmPairs2Domain$contrasts%>%
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
emmPairs2Domain_effect_val = eff_size(emmPairs2Domain, 
                                  sigma = sigma(PairsXDomain_mod), 
                                  edf = df.residual(PairsXDomain_mod),
                                  method = "identity")%>%
  as.data.frame()%>%
  rename(effect.size = "effect.size",
         e.s_lower.CL = "lower.CL",
         e.s_upper.CL= "upper.CL",
  )%>%
  mutate(contrast = substr(contrast, 2, nchar(as.character(contrast))-1))

##REPORT:aesthetic agreement####
save(anova_interR_PairsXDomain_unballanced_val, 
     anova_interR_PairsXDomain_unballanced_val_out,
     effect_size_anova_val,
     effect_size_anova_val_out,
     emmPairs_estimates_csv_val,
     emmPairs2Domain_output_val,
     emmPairs2Domain_effect_val,file = sprintf("%s/%s/02_Sutherland_2020/03_report_aesthetic_agreement_val.Rdata",wdOA,wdOA_output))

#Prepare CSV
emmPairs2Domain_stats = merge(emmPairs2Domain_estimates%>%
                                select(-c(SE, df,lower.CL_z, upper.CL_z,lower.CL_r,  upper.CL_r)),emmPairs2Domain_output_val,by = c("contrast", "Domain"))
emmPairs2Domain_csv = 
  merge(emmPairs2Domain_stats,emmPairs2Domain_effect_val%>%select(-c(SE, df)),by = c("contrast","Domain"))%>%
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
write_csv(emmDomain_estimates_csv,sprintf("%s/%s/02_Sutherland_2020/03_sup_pairwise_domain_estimates_Sutherland_2020.csv",wdOA,wdOA_output))
write_csv(emmDomain_csv,sprintf("%s/%s/02_Sutherland_2020/03_sup_pairwise_domain_statistics_Sutherland_2020.csv",wdOA,wdOA_output))
write_csv(emmPairs_estimates_csv_val,sprintf("%s/%s/02_Sutherland_2020/03_sup_pairwise_pairs_estimates_Sutherland_2020.csv",wdOA,wdOA_output))
write_csv(emmPairs_csv,sprintf("%s/%s/02_Sutherland_2020/03_sup_pairwise_pairs_statistics_Sutherland_2020.csv",wdOA,wdOA_output))
write_csv(emmPairs2Domain_estimates_csv,sprintf("%s/%s/02_Sutherland_2020/03_sup_pairwise_pairs2domain_estimates_Sutherland_2020.csv",wdOA,wdOA_output))
write_csv(emmPairs2Domain_csv,sprintf("%s/%s/02_Sutherland_2020/03_sup_pairwise_pairs2domain_statistics_Sutherland_2020.csv",wdOA,wdOA_output))

#save pariswise agreeemnt
write_csv(InterR_random_noOut,sprintf("%s/%s/02_Sutherland_2020/03_pairwise_agreement_Sutherland2020.csv",wdOA,wdOA_output))

# ####Figure 2: before editorial####
# PairsXDomain_mod_r = lm(r_twin_PP ~ +Pairs * Domain, data = InterR_random_noOut)
# #final plot (r transformed)
# p_pairwise_r_val = emmip(PairsXDomain_mod_r, Pairs ~ Domain) +
#   theme_bw(base_size = 14) + 
#   ylim(c(0.2,0.8))+
#   scale_color_viridis_d() +
#   labs(y = "pairwsie agreement",
#        x = "Visual domains")
# 
# #save for later plotting
# save(p_pairwise_r_val,file =sprintf("%s/%s/02_Sutherland_2020/03_Figure2b_pairwise_agreement.Rdata",wdOA,wdOA_output))

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
p_pairwise_r_val_revised = ggplot(InterR_random_noOut_mean) + 
  geom_jitter(data = InterR_random_noOut,aes(Domain,r_twin_PP, color = Pairs),
              alpha = .05,
              width = 0.1,
              size =.50) +
  geom_line(aes(Domain, r_twin_PP, color = Pairs, group = Pairs))+
  geom_point(aes(Domain, r_twin_PP, color = Pairs), size = 1.5)+
  geom_errorbar(aes(x =Domain, ymin = lower.ci, ymax = upper.ci, color = Pairs), width = .1)+
  scale_color_viridis_d() + #note that final color has been changed after Review 01
  scale_y_continuous(lim= c(0,1), breaks = seq(0,1,by = .1))+
  labs(y = "pairwsie agreement",
       x = "Visual domains")+
  theme_light()

p_pairwise_r_val_revised_zoom = ggplot(InterR_random_noOut_mean) + 
  geom_line(aes(Domain, r_twin_PP, color = Pairs, group = Pairs))+
  geom_point(aes(Domain, r_twin_PP, color = Pairs), size = 1.5)+
  geom_errorbar(aes(x =Domain, ymin = lower.ci, ymax = upper.ci, color = Pairs), width = .1)+
  scale_color_viridis_d() + #note that final color has been changed after Review 01
  scale_y_continuous(lim= c(.55,.85), breaks = seq(.55,.85,by = .1))+
  labs(y = "pairwsie agreement",
       x = "Visual domains")+
  theme_light()

##save for later plotting
save(p_pairwise_r_val_revised,p_pairwise_r_val_revised_zoom,file =sprintf("%s/%s/02_Sutherland_2020/03_Figure2b_pairwise_agreement_revised.Rdata",wdOA,wdOA_output))

##Supplementary File 6 (Editorial Request)####
#save source data
write_csv(InterR_random_noOut_mean %>% 
            mutate(pairwise.r = round(r_twin_PP,3),lower.ci = round(lower.ci,3),upper.ci = round(upper.ci,3))%>% 
                     select(Domain,Pairs, pairwise.r,lower.ci,upper.ci)
          ,sprintf("%s/%s/02_Sutherland_2020/03_Fig_F2b_summary_pairwise_agreement_revised_sourceData.csv",wdOA,wdOA_output))

##Supplementary File 7 (Additional Editorial Request)####
#save source data
write_csv(InterR_random_noOut %>% 
            mutate(pairwise.r = round(r_twin_PP,3))%>% 
            select(Domain,Pairs, pairwise.r),
          sprintf("%s/%s/02_Sutherland_2020/03_Fig_F2_individual_pairwise_agreement_revised_sourceData.csv",wdOA,wdOA_output))


#Supplementary Fig 3####
#intra rater dist
p1_val = ggplot(InterR_random_noOut, aes( y = r_twin_PP, fill = Pairs)) +
  geom_density(alpha = .5) + 
  scale_fill_viridis_d() +
  facet_wrap(~Domain) +
  labs(y = "inter-rater r") + 
  theme_classic(base_size = 14)
#intra rater z transform dist
p2_val = ggplot(InterR_random_noOut, aes( y = z_twin_PP, fill = Pairs)) +
  geom_density(alpha = .5) + 
  scale_fill_viridis_d() +
  labs(y = "inter-rater z") + 
  facet_wrap(~Domain) +
  theme_classic(base_size = 14)
#FIg2a for z values
p_pairwise_z_val = emmip(PairsXDomain_mod, Pairs ~ Domain) +
  theme_bw(base_size = 16) + 
  scale_color_viridis_d() +
  labs(y = "pairwsie agreement",
       x = "Visual domains")
#marginal mean effects domains
p3_val = plot(emmDomain, comparisons = TRUE) + 
  theme_bw(base_size = 14) +
  labs(x = "Marginal means inter-rater Z", y = "Domains")
p4_val = plot(emmPairs, comparisons = TRUE) + 
  theme_bw(base_size = 14) +
  labs(x = "Marginal means inter-rater Z", y = "Domains")
#marginal mean effects pairs|domains
p5_val = plot(emmPairs2Domain, comparisons = TRUE) + 
  theme_bw(base_size = 14) +
  labs(x = "Marginal means inter-rater Z", y = "Pairs")

#save for later plotting
save(p1_val,p2_val,p3_val,p4_val,p5_val,p_pairwise_z_val, file =sprintf("%s/%s/02_Sutherland_2020/03_sup_FS6_pairwise_agreement_val.Rdata",wdOA,wdOA_output))
