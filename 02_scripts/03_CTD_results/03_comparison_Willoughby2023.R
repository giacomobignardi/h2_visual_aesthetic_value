# Author: Giacomo Bignardi 
# Date: 18 12 2021
# Last modified: 18-09-2023
#
# Description: Comparison of ACE estimates as obtained in this study and meta-analytic ACE estimates for other human complex traits (adapted from http://match.ctglab.nl/#/multiple/reported_ace)
# Adapted to: https://www.nature.com/articles/s43586-022-00191-x after review
# Load Libraries & Options
# load packages
library(readr)
library(tidyverse)
library(patchwork)

#clear wd
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_Data = "01_input"
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"
wdOA_ImageOutput = "05_Figures"

####DATA####
Bignardi = read_csv(sprintf("%s/%s/03_CTD_results/SF2_ACE_results.csv",wdOA,wdOA_output))
Bignardi_summary = Bignardi%>%
  filter(comparison == "AE" | comparison == "ACE" )%>%
  filter(!(domain == "abstract" & domain == "abstracts" & facet == "taste_typicality"))%>%
  mutate(A = as.numeric(SA),
         C_D = as.numeric(SC_SD),
         E= as.numeric(SE))%>%
  dplyr::select(comparison ,domain, facet,sample, A,C_D,E)%>%
  #select full ACE model or AE model if full model was ADE
  mutate( plot = ifelse((domain == "faces" & facet == "taste_typicality" & sample == "Sutherland et al. 2020") |
                (domain == "abstracts"& facet == "evaluation_bias"& sample == "Germine et al. 2015") |
                comparison == "ACE" , 1, 0)) %>% 
  filter(plot == 1) %>% 
  filter(!(domain == "abstract"& facet == "taste_typicality")) %>% 
  group_by(facet)%>%
  summarise(A = round(mean(A),2),
            C = round(mean(C_D),2),
            E= round(mean(E),2))%>%
  rename(category = "facet") %>% 
  mutate(category = factor(category, levels = c("taste_typicality","evaluation_bias"), labels = c("Taste-typicality","Evaluation-bias")))


Bignardi_absPref = Bignardi%>%
  filter(comparison == "CE")%>%
  filter(domain == "abstract" &  facet == "taste_typicality")%>%
  mutate(A = as.numeric(SA),
         C = as.numeric(SC_SD),
         E= as.numeric(SE))%>%
  dplyr::select(domain, facet,sample, A,C,E)%>%
  group_by(facet)%>%
  summarise(A = round(mean(A),2),
            C = round(mean(C),2),
            E= round(mean(E),2))%>%
  rename(category = "facet")%>%
  mutate(category = factor(category, levels = c("taste_typicality"), labels =c("Taste-typicality abstract")))
  
#compare results with results obtained for other domains (Polderman et al. 2015)
category = c(
  "Heart functions",
  "Weight maintenance",
  "Height",
  "Blood pressure",
  "Immunological system functions",
  "Bipolar affective disorder",
  "Depressive episode",
  "Eating disorders",
  "Anxiety disorders (non-phobic)",
  "Schizophrenia",
  "Alcohol-related disorders",
  "Higher-level cognitive functions",
  "Education",
  "Memory functions",
  "Brain structure (0 to 11 years)",
  "Brain structure (18 to 64 years)",
  "Temperament and personality",
  "Societal attitudes",
  "Religion and spirituality",
  "Social relationships (male)",
  "Social relationships (female)"
  
)

# Table 2 https://www.nature.com/articles/s43586-022-00191-x/tables/2

A = c(.29,
      .63,
      .80,
      .47,
      .42,
      .68,
      .34,
      .40,
      .39,
      .77,
      .42,
      .63,
      .51,
      .48,
      .83,
      .49,
      .45,
      .31,
      .40,
      .29,
      .34)
C = c(.16,
      .10,
      .10,
      .11,
      .08,
      .15,
      .11,
      .04,
      .05,
      .01,
      .12,
      .11,
      .28,
      .08,
      .09,
      .10,
      .12,
      .20,
      .20,
      .11,
      .12
      )
E= (1-(A+C))

Willoughby_2023 = data.frame(category, A,C,E)
comparison1 = rbind(Willoughby_2023, Bignardi_summary,Bignardi_absPref)
#filter for cognitive traits
comparison1 = comparison1 %>% filter(!(category %in% c("Heart functions", "Weight maintenance","Height", "Blood pressure", "Immunological system functions")))
comparison1 = comparison1%>%
  pivot_longer(A:E, names_to = "estimates", values_to = "value")%>%
  mutate(estimates = factor(estimates,
                            levels = rev(c("A","C","E"))))

ordered_category = comparison1 %>% arrange(estimates,value) %>% filter(estimates == "A") %>% pull(category)
comparison1 = comparison1 %>% mutate(category = factor(category, levels = c(ordered_category)))
colors5prof <- rev(c("#914fa3","#41af76","#697F8D" ))
h2_comparison  = ggplot()+
  geom_bar(data = comparison1, aes(x=value, y=category, fill=estimates), stat="identity", alpha =.8)+
  labs( y= "Estimates",
       x = "",
       fill = "Etiological Components")+
  theme_classic(base_size = 12)+
  geom_text(data = comparison1%>%filter(estimates == "A" ),aes(y = category, x =  value, label = value),
            hjust=1.1,
            color="white",
            size=3.5)+
  theme(axis.ticks.y = element_line(linewidth =.5),
        legend.position = "bottom")+
  scale_x_continuous(breaks = c(0, 1))+
  geom_vline(xintercept=0)+
  #Change the order of the levels in the breaks
  scale_fill_manual(values = colors5prof, breaks = levels(comparison1$estimates))

#SAVE####
save(h2_comparison,file = sprintf("%s/%s/03_CTD_results/03_CTD_comparison.Rdata",wdOA, wdOA_output))