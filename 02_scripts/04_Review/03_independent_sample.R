#Author: Giacomo Bignardi
#Adapted from: Nichola Burton (Sutherland et al., 2022, PNAS)
#last modified:2023-10-10
#
#Description: sensitivity analysis: test-retest
# ------------------------------------------------------------------------------
#load packages
library(tidyverse)
library(psych)
library(tidylog)
library(readr)

#clean working enviroment 
rm(list = ls())

#clean working environment 
rm(list = ls())

#clear environment
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

#load functions:
source(sprintf("%s/%s/functions/SDcheck.R", wdOA,wdOA_scripts))
source(sprintf("%s/%s/functions/IntraRater.R", wdOA,wdOA_scripts))
source(sprintf("%s/%s/functions/MM2.R", wdOA,wdOA_scripts))

#load data:
online1 = read_csv(sprintf("%s/%s/Sutherland_2020/04_Review/01_Online_ExclusionsSession1.csv", wdOA,wdOA_Data))
online2 = read_csv(sprintf("%s/%s/Sutherland_2020/04_Review/04_Online_Session2Exclusions.csv", wdOA,wdOA_Data))
#apply exclusion criteria as in Sutherland et al.
online2 = online2  %>% filter(
  Exclude_RT == 0,
  Exclude_Var == 0,
  Exclude_SC == 0)

#retain only participants included at Session 2:
online2_id = online2 %>% pull(Code)
online1 = online1 %>% filter(Code %in% online2_id)

#EXTRACT RATINGS DATA
#locate the column number of the first rating stimulus and subset all ratings
f1 = match("Anthony_Rimmer_11_oval_1",names(online1))
allRatings1 = online1[,c(1:2,f1:(f1+523))]
allRatings2 = online2[,c(2:3,f1:(f1+523))]

#locate column number of first stimulus from each test
a1 = match("Anthony_Rimmer_11_oval_1",names(allRatings1))
d1 = match("Aaron_Mink_9_oval_1", names(allRatings1))
s1 = match("sc001_1", names(allRatings1))

a2 = match("Anthony_Rimmer_11_oval_1",names(allRatings2))
d2 = match("Aaron_Mink_9_oval_1", names(allRatings2))
s2 = match("sc001_1", names(allRatings2))

#break up the ratings into the different tests
attrRatings1 = allRatings1[c(1:2,a1:(a1+149))]
domRatings1 = allRatings1[c(1:2,d1:(d1+149))]
sceneRatings1 = allRatings1[c(1:2,s1:(s1+73))]

attrRatings1_l = attrRatings1 %>% 
  pivot_longer(3:ncol(attrRatings1),names_to = "item", values_to = "value") %>% 
  mutate(is_repeated = substr(item,nchar(item),nchar(item)))
domRatings1_l = domRatings1 %>% 
  pivot_longer(3:ncol(domRatings1),names_to = "item", values_to = "value") %>% 
  mutate(is_repeated = substr(item,nchar(item),nchar(item)))
sceneRatings1_l = sceneRatings1 %>% 
  pivot_longer(3:ncol(sceneRatings1),names_to = "item", values_to = "value") %>% 
  mutate(is_repeated = substr(item,nchar(item),nchar(item)))

ratings1 = rbind(attrRatings1_l %>% mutate(domain = "fac"),domRatings1_l%>% mutate(domain = "dom"),sceneRatings1_l%>% mutate(domain = "sce")) %>% 
  arrange(Code) %>% 
  mutate(Day = 1)

attrRatings2 = allRatings2[c(1:2,a2:(a2+149))]
domRatings2 = allRatings2[c(1:2,d2:(d2+149))]
sceneRatings2 = allRatings2[c(1:2,s2:(s2+73))]


attrRatings2_l = attrRatings2 %>% 
  pivot_longer(3:ncol(attrRatings2),names_to = "item", values_to = "value") %>% 
  mutate(is_repeated = substr(item,nchar(item),nchar(item)))
domRatings2_l = domRatings2 %>% 
  pivot_longer(3:ncol(domRatings2),names_to = "item", values_to = "value") %>% 
  mutate(is_repeated = substr(item,nchar(item),nchar(item)))
sceneRatings2_l = sceneRatings2 %>% 
  pivot_longer(3:ncol(sceneRatings2),names_to = "item", values_to = "value") %>% 
  mutate(is_repeated = substr(item,nchar(item),nchar(item)))

ratings2 = rbind(attrRatings2_l %>% mutate(domain = "fac"),domRatings2_l%>% mutate(domain = "dom"),sceneRatings2_l%>% mutate(domain = "sce")) %>% 
  arrange(Code) %>% 
  mutate(Day = 2)

ratings = rbind(ratings1,ratings2)

#compute intra rarter for further exclusion
intraRel_fac = intraRater(ratings %>% filter(domain == "fac") %>% select(item,value,Code,Day))
intraRel_sce = intraRater(ratings %>% filter(domain == "sce") %>% select(item,value,Code,Day))

intraRel_fac_summ = intraRel_fac$intraRater %>% mutate(domain ="fac")
intraRel_sce_summ = intraRel_sce$intraRater %>% mutate(domain ="sce")

ratings_final =
  merge(ratings %>% rename(Subj = Code), rbind(intraRel_fac_summ,intraRel_sce_summ), by = c("Subj","domain")) %>% 
  filter(novariance ==0)

#demographics
dem = online1 %>% filter(Code %in% (ratings_final %>% distinct(Subj) %>% pull(Subj))) 
table(dem$Sex)
dem %>% summarise(mean(How.old.are.you..in.years..),sd(How.old.are.you..in.years..), range(How.old.are.you..in.years..))

#AESTHETIC AGREEMENT####
ratings_final %>% distinct(Subj)
#average within rater across repeated images
ratings_final = ratings_final %>% 
  mutate(item = substr(item,1,nchar(item)-1)) %>% 
  pivot_wider(names_from = is_repeated, values_from = value) %>% 
  mutate(`2` = ifelse(is.na(`2`),`1`,`2`)) %>% 
  mutate(value =  (`1`+`2`)/2) %>% 
  select(-c(`1`,`2`))
#create random pairs
n = nrow(ratings_final %>% distinct(Subj))
pairs = n / 2
pair_id = c(rep(1:pairs), rep(1:pairs))
set.seed(42)  
pair_id_rnd = sample(pair_id)
pair_ids =ratings_final %>% 
  distinct(Subj) %>% 
  mutate(pair_id = pair_id_rnd) %>% 
  mutate(pair = ifelse(duplicated(pair_id), 2, 1))

ratings_final %>% filter(intraR <=.5) %>% distinct(Subj)
#final dataframe for aesthetic agreement analysis
paired_ratings = merge(ratings_final,pair_ids, by = "Subj") %>% 
  as.data.frame() 

paired_ratings_wid1 = paired_ratings %>%
  filter(Day == 1,intraR >=.5) %>% #for simplicity take only non repeated images
  select(-c(Subj,
            How.old.are.you..in.years..,intraR)) %>% 
  pivot_wider(names_from = "pair", values_from = "value")

#calculate aesthetic agreement at day 1
AA1 = paired_ratings_wid1%>%
  group_by(pair_id,domain)%>%
  summarise(r_PP = cor(`1`,`2`), z_PP = psych::fisherz(cor( `1`,`2`)))%>%ungroup()

#get 
sum_AA1 = AA1 %>% 
  group_by(domain) %>% 
  rstatix::get_summary_stats(z_PP) %>% 
  mutate(mean_r = psych::fisherz2r(mean),
         sd_r = psych::fisherz2r(sd))

sum_AA1 %>% select(domain,mean_r, sd_r)

#calculate aesthetic agreement at day 2
paired_ratings_wid2 = paired_ratings %>%
  filter(Day == 2,intraR >=.5) %>% 
  select(-c(Subj,
            How.old.are.you..in.years..,intraR))%>% #for simplicity take only non repeated images
  pivot_wider(names_from = "pair", values_from = "value")


AA2 = paired_ratings_wid2%>%
  group_by(pair_id,domain)%>%
  summarise(r_PP = cor(`1`,`2`), z_PP = psych::fisherz(cor( `1`,`2`)))%>%ungroup()

#get 
sum_AA2 = AA2 %>% 
  group_by(domain) %>% 
  rstatix::get_summary_stats(z_PP) %>% 
  mutate(mean_r = psych::fisherz2r(mean),
         sd_r = psych::fisherz2r(sd))

sum_AA2 %>% select(domain,mean_r, sd_r)

AA = merge(AA1,AA2, by = c("pair_id", "domain")) %>% as.data.frame()

#calculate test-retest reliability
AA_fac = AA %>% filter(domain == "fac")
AA_sce = AA %>% filter(domain == "sce")

##between day####
AA_cor_faces = cor.test(AA_fac$z_PP.x, AA_fac$z_PP.y)
AA_cor_scenes = cor.test(AA_sce$z_PP.x, AA_sce$z_PP.y)
#calculate test-retest reliability (ICC(2,1))
AA_ICC_faces = ICC(AA_fac %>% select(z_PP.x,z_PP.y))
AA_ICC_scenes = ICC(AA_sce %>% select(z_PP.x,z_PP.y))


#Taste-typicality####
#faces
fac_mat_ratings1 = ratings_final %>% filter(Day == 1, domain == "fac",intraR >=.5) %>% 
  select(item,Subj,value) %>% 
  pivot_wider(names_from= Subj, values_from = value)
fac_mat_ratings2 = ratings_final %>% filter(Day == 2, domain == "fac",intraR >=.5) %>% 
  select(item,Subj,value) %>% 
  pivot_wider(names_from= Subj, values_from = value)

tt_fac1 = MM2(fac_mat_ratings1)
tt_fac2 = MM2(fac_mat_ratings2)

tt_fac = merge(tt_fac1$summary,tt_fac2$summary, by = "Sub") %>% as.data.frame()

#scenes
sce_mat_ratings1 = ratings_final %>% filter(Day == 1, domain == "sce",intraR >=.5) %>% 
  select(item,Subj,value) %>% 
  pivot_wider(names_from= Subj, values_from = value)
sce_mat_ratings2 = ratings_final %>% filter(Day == 2, domain == "sce",intraR >=.5) %>% 
  select(item,Subj,value) %>% 
  pivot_wider(names_from= Subj, values_from = value)

tt_sce1 = MM2(sce_mat_ratings1)
tt_sce2 = MM2(sce_mat_ratings2)

tt_sce = merge(tt_sce1$summary,tt_sce2$summary, by = "Sub") %>% as.data.frame()

##between day####
tt_cor_faces = cor.test(as.numeric(tt_fac$mm2_z.x), as.numeric(tt_fac$mm2_z.y)) 
tt_ICC_faces = ICC(tt_fac %>% select(mm2_z.x,mm2_z.y) %>% mutate(mm2_z.x = as.numeric(mm2_z.x), mm2_z.y = as.numeric(mm2_z.y)))

tt_cor_scenes = cor.test(as.numeric(tt_sce$mm2_z.x), as.numeric(tt_sce$mm2_z.y)) 
tt_ICC_scenes = ICC(tt_sce %>% select(mm2_z.x,mm2_z.y) %>% mutate(mm2_z.x = as.numeric(mm2_z.x), mm2_z.y = as.numeric(mm2_z.y)))


#Evaluation-bias####
eb = ratings_final %>% filter(intraR >=.5) %>% 
  group_by(Subj,Day, domain) %>% 
  summarise(eval_bias = mean(value)) %>% 
  pivot_wider(names_from = Day, values_from = eval_bias) %>% 
  ungroup()

eb_fac = eb %>% filter( domain == "fac")
eb_sce = eb %>% filter( domain == "sce")

#replication phenotypic correlations
tt_cor = merge(tt_fac,tt_sce, by = "Sub")

#phenotypic correlation (taste-typicality faces scenes)
##between day####
eb_cor_faces   = cor.test(eb_fac$`1`,eb_fac$`2`)
eb_cor_scenes  = cor.test(eb_sce$`1`,eb_sce$`2`)
eb_ICC_faces   =  ICC(eb %>% filter( domain == "fac") %>% select(`1`,`2`))
eb_ICC_scenes  = ICC(eb %>% filter( domain == "sce") %>% select(`1`,`2`))

#Z-FISHER####
#test for difference between taste-typicality and aesthetic agreement between days correlations 
AA_fac_cor = AA_fac %>% rename(PAz_1 = z_PP.x,PAz_2 = z_PP.y) %>% select(pair_id,PAz_1,PAz_2)
AA_sce_cor = AA_sce %>% rename(PAz_1 = z_PP.x,PAz_2 = z_PP.y) %>% select(pair_id,PAz_1,PAz_2)

tt_fac_cor = tt_fac %>% rename(tt_1 = mm2_z.x,tt_2 = mm2_z.y) %>% select(Sub,tt_1,tt_2) %>% mutate(tt_1 = as.numeric(tt_1), tt_2 = as.numeric(tt_2))
tt_sce_cor = tt_sce %>% rename(tt_1 = mm2_z.x,tt_2 = mm2_z.y) %>% select(Sub,tt_1,tt_2) %>% mutate(tt_1 = as.numeric(tt_1), tt_2 = as.numeric(tt_2))


cocor_fac = cocor::cocor(~  tt_1 + tt_2 |PAz_1 + PAz_2 ,
                         data = list(tt_fac_cor, AA_fac_cor))

cocor_sce = cocor::cocor(~tt_1 + tt_2 | PAz_1 + PAz_2 ,
                         data = list(tt_sce_cor, AA_sce_cor))

##NOTE MANUAL REPORT OF COCOR####
cocor_fac
cocor_fac_est = data.frame(delta_r = 0.1419, tt_sample = 75, PA_sample= 36, fisher_z = 1.6072, p = 0.1080)

cocor_sce
cocor_sce_est = data.frame(delta_r = 0.0832, tt_sample = 78, PA_sample= 39, fisher_z = 0.9794, p = 0.3274)

save(eb_cor_faces,tt_cor_faces,eb_cor_scenes,tt_cor_scenes, AA_cor_faces,AA_cor_scenes,
     eb_ICC_faces,tt_ICC_scenes,eb_ICC_faces,eb_ICC_scenes, AA_ICC_faces,AA_ICC_scenes,
     cocor_fac_est,cocor_sce_est, file = sprintf("%s/%s/04_review/03_report_betweenday.RData",wdOA,wdOA_output))#save output


#S9####
#rename eb
eb_fac_cor = eb_fac %>% rename(eb_1_fac = `1` ,eb_2_fac = `2`, Sub = Subj ) %>% select(-domain)
eb_sce_cor = eb_sce %>% rename(eb_1_sce = `1` ,eb_2_sce = `2`, Sub = Subj ) %>% select(-domain)

tt_fac_cor = tt_fac_cor %>% rename(tt_1_fac = tt_1,tt_2_fac  = tt_2) 
tt_sce_cor = tt_sce_cor %>% rename(tt_1_sce = tt_1,tt_2_sce = tt_2) 

pheno_cor_tt = merge(tt_fac_cor,tt_sce_cor, by = "Sub")
pheno_tt_1 = cor.test(pheno_cor_tt$tt_1_fac,pheno_cor_tt$tt_1_sce)
pheno_tt_2 = cor.test(pheno_cor_tt$tt_2_fac,pheno_cor_tt$tt_2_sce)

pheno_cor_eb = merge(eb_fac_cor,eb_sce_cor, by = "Sub")
pheno_eb_1 = cor.test(pheno_cor_eb$eb_1_fac,pheno_cor_eb$eb_1_sce)
pheno_eb_2 =cor.test(pheno_cor_eb$eb_2_fac,pheno_cor_eb$eb_2_sce)

#SAVE####
save(pheno_tt_1,pheno_tt_2,
     pheno_eb_1,pheno_eb_2,
     file = sprintf("%s/%s/04_review/03_report_pheno_cor_replication.RData",wdOA,wdOA_output))#save output


