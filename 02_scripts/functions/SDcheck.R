#Author: Giacomo Bignardi
#Date: 15/04/2021
#
#
#
#Create a function to check whether one participant's distibution of ratings has sd = 0

sdCheck = function(df){
  
  require(tidyverse)
  names(df) = c("Obj","Rating","FamId","Block")
  df = df%>%mutate(Rating = as.numeric(Rating), Block = as.factor(Block))
  #number of unique participants
  N = length(unlist(unique(df[,3])))
  
  sds = df%>%group_by(FamId, Block)%>%summarise(sd = sd(Rating))
  
  N_removed = unique(sds[sds$sd == 0,]$FamId) #save for later
  
  N_final= N - length(N_removed) #new N
  
  list(InitalSample = N,
       finalSample = N_final,
       IDs = N_removed)
}

#InitalSample = n of participants
#finalSample = n of participants after removing NA and SD = 0
#IDs = participant IDs (and NA if at least on NA for sd has been computed)
