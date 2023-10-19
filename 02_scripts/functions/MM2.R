#Author: Giacomo Bignardi
#Date: 15/04/2021
#
#
#
#Create a function to compute mean minus 2

#laod MM1 function
MM2 = function(X){
  
  require(tidyverse)
  
  #drop stim name
  X = X[,-1]
  #save subj id
  id = colnames(X)
  
  #Compute means minus two
  means = c()
  for(i in 1:ncol(X)){ #iterate for all the rows (items)
    mean_i = rowMeans(X%>%select(!starts_with(sprintf("%s_",substr(id[i],1,nchar(id[i])-2)))))
    means = cbind(means,mean_i)
  }
  
  #Compute correlations between means minus two and participants ratings
  mm2 = c()
  mm2_low = c()
  mm2_high = c()
  for(j in 1:ncol(X)){
    cor = cor(pull(X[,j]),means[,j]) #correlations between participants ratings and their MM2
    mm2 = c(mm2,as.numeric(cor))
  }
  
  #transform correlations into z scores (r to z Fisher)
  zScores = psych::fisherz(mm2)
  #calculate mean and CI
  mean_mm2 = mean(zScores, na.rm = T)
  sd_mm2 = sd(zScores, na.rm = T)
  
  #save list (convert Z to r again)
  MM2List = list(MM1 = cbind(mm2 = psych::fisherz2r(mean_mm2),sd = sd_mm2), 
                 summary = as.data.frame(cbind(Sub = id,mm2_r = mm2, mm2_z = zScores),stringsAsFactors = F))
  names(MM2List) = c("MM2", "summary")
  MM2List
}

