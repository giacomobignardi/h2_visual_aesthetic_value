#Author: Giacomo Bignardi
#Date: 15/04/2021
#
#
#
#Create a function to check whether one participant's distibution of ratings has sd = 0
intraRater = function(df){
  
  require(tidyverse)
  #it require(psych)
  #it require(Rmisc)
  
  names(df) = c("Obj","Rating","Subj","Block")
  
  N_subj = length(unlist(unique(df[,3])))
  N_items = length(unlist(unique(df[,1])))
  
  df_intra = df%>%spread(Block,Rating)
  
  #compute: intrarater correlations
  intrarater_df_long = df_intra%>%
    dplyr::group_by(Subj)%>%
    dplyr::summarise(cor = cor(`1`,`2`), use = "complete.obs")%>%
    dplyr::select("Subj", "cor")
  
  intrarater_df_long = as.data.frame(intrarater_df_long)
  
  names(intrarater_df_long) = c("Subj", "intraR")
  
  intrarater_df_long$intraR =  as.numeric(intrarater_df_long$intraR)
  
  #calculate intra rater correlation coefficient:
  #deal with 1 for Fisher Z transformation.Fisher Z transformation needed to obtain a distribution from which one can estimate CI 
  intrarater_df_long$intraR = ifelse(intrarater_df_long$intraR >= 0.99999999999, 0.99999,  intrarater_df_long$intraR)
  intrarater_df_long$novariance = ifelse(is.na(intrarater_df_long$intraR), 1,  0)
  sumy_intrarater =  psych::fisherz2r(summary(psych::fisherz(intrarater_df_long[intrarater_df_long$novariance==0,]$intraR)))
  sd = psych::fisherz2r(sd(psych::fisherz(intrarater_df_long[intrarater_df_long$novariance==0,]$intraR)))
  #psych::fisherz2r(Rmisc::CI(psych::fisherz(intrarater_df_long[intrarater_df_long$novariance==0,]$intraR),ci = 0.95))
    
  N_subj_novariance = nrow(intrarater_df_long[intrarater_df_long$novariance==0,])
  
  #Final output
  list(intraRater = intrarater_df_long,
       intraRater_Summary =sumy_intrarater,
       intraRater_sd = sd,
       TotalN = N_subj,
       FinalN = N_subj_novariance,
       N_items = N_items)
} 

#[[1]]= df with:
##Subj = subID
##intraR = intrarater correlations per subject
#[[2]]= average intrarater correlation
#[[3]]= fsample used to calculate intrarater correlation
