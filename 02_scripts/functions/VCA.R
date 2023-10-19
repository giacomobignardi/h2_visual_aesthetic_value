#Author: Giacomo Bignardi
#Date: 15-04-2021
# NOTE: and updated version of this function can be found at https://github.com/giacomobignardi/empirical-aesthetics-VCA/tree/main/01_VCA/02_scripts/R/functions 
#
#
#Create a function to compute variance component from a multi level model

#laod MM1 function
VCA = function(mlm) {
  #Function to calcucalte Variance Component from a MLM (adpated from Sutherland et al, 2020 Nichola Burton and Martinez et al. 2020)
  #extract summary variance
  sum_mlm = summary(mlm)
  
  #calculate total variance
  total_variance = 
    (sum_mlm$varcor$`Sub:Obj`[1])+
    (sum_mlm$varcor$`Sub`[1])+
    (sum_mlm$varcor$`Obj`[1])+
    (sum_mlm$varcor$`Block`[1])+
    (sum_mlm$varcor$`Block:Obj`[1])+
    (sum_mlm$varcor$`Block:Sub`[1])+
    (sum_mlm$sigma^2)
  
  #calculate variance components
  VCSub =  (sum_mlm$varcor$`Sub`[1]) / total_variance
  VCSubXObj =  (sum_mlm$varcor$`Sub:Obj`[1]) / total_variance
  VCObj =  (sum_mlm$varcor$`Obj`[1]) / total_variance
  VCBlock = (sum_mlm$varcor$`Block`[1]) / total_variance
  VCBlockXSub = (sum_mlm$varcor$`Block:Sub`[1]) /total_variance
  VCBlockXObj = (sum_mlm$varcor$`Block:Obj`[1]) / total_variance
  VCResidual =  (sum_mlm$sigma^2) /  total_variance
  

   c(VCObj, VCSub, VCBlock, VCSubXObj, VCBlockXSub, VCBlockXObj, VCResidual)
  
}

