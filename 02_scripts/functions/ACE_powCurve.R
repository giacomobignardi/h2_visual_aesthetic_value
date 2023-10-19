#Author: Giacomo Bignardi
#Inspiref from: Verhulst, 2017; Behav Gen
#Date: 2023-05-29
#
#Description: Step 4 of power analysis for ACE model in Verhulst, 2017; Behav Gen
#--------------------------------------------------------------------------------
#Multiply WNCP per number of family and compute power
powCurveACE = function(quant,maxN = 5000, pcrit = .10, df = 1){
  criticalX = qchisq(1- pcrit, df)
  WncpA_seq = quant$WncpA* (seq(1,maxN,1))
  WncpC_seq = quant$WncpC* (seq(1,maxN,1))
  #compute β 
  BetaA_seq = pchisq(criticalX, df = df, WncpA_seq)
  BetaC_seq = pchisq(criticalX, df = df, WncpC_seq)
  #compute statistical power
  PowerA_seq =  1 - BetaA_seq 
  PowerC_seq =  1 - BetaC_seq 
  
  data.frame(family_n = (seq(1,maxN,1)),
             A = PowerA_seq,
             C = PowerC_seq) %>% 
    pivot_longer(names_to = "component", values_to = "power", c(A:C))
}
