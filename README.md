## Genetic effects on variability in visual aesthetic evaluations are partially shared across visual domains

[![DOI](https://zenodo.org/badge/707185913.svg)](https://zenodo.org/doi/10.5281/zenodo.10251278)

This repository contains the data, processed data, and the code to replicate the analysis made in:

> Bignardi, G., Smit, D., Vessel, E. A., Trupp, M. D., Ticini, L. F., Fisher, S. E., & Polderman, T. J. (2023). Shared and Distinct Genetic Sources of Variability in Visual Aesthetic Value. https://doi.org/10.31234/osf.io/79nbq 

## Repository

The following directories contain:
+ 01_input: input data to reproduce the analysis (see 04_metadata). 
+ 02_scripts: includes scripts for the analysis
+ 03_outputs: includes outputs generated via 02_scripts
+ 04_metadata: includes instructions to download original data from [Germine et al. (2015)](https://doi.org/10.1016/j.cub.2015.08.048) and [Sutherland et al. (2020)](https://doi.org/10.1073/pnas.1920131117).
+ 05_figures: includes unedited Figures

Please note that the raw data supporting the findings of this study are fully available at https://osf.io/c3hz6/ and https://osf.io/35zf8/?view_only=e76c6755dcea4be2adc5b075cae896e8. 
These data were originally collected by [Germine et al. (2015)](https://doi.org/10.1016/j.cub.2015.08.048) and [Sutherland et al. (2020)](https://doi.org/10.1073/pnas.1920131117). This study is entirely built upon these data alone. 


## Manuscript reproducibility

Code to reproduce specific results or figures included in the subsections of the main and supplementary material can be found under the specific directories:

### Inter-individual differences in visual aesthetic evaluations.
+ Text: /02_scripts/03_CTD_results/00_vizualize.R
+ Table 1:  
+ /02_scripts/01_Germine_2015/01_prepare_df.R  
+ /02_scripts/02_Sutherland_2020/01_prepare_df.R

### Measures of inter-individual differences in aesthetic evaluation.
+ Text: /03_outputs/report/01_PCA.Rmd
+ Fig.1b: /02_scripts/03_CTD_results/00_vizualize.R
+ Fig.1d: /02_scripts /01_Germine_2015/05_phenotyping_PCA.R

### Monozygotic twins show higher pairwise aesthetic agreement than dizygotic twins and unrelated pairs.
+ Text: /03_outputs/report/02_aesthetic_agreement.Rmd
+ Fig.2: /02_scripts/03_CTD_results/00_vizualize.R

### Univariate CTD estimates genetic effects for all dimensions of aesthetic value except taste-typicality for abstract images.
+ Text: /03_outputs/report/03_univariate_CTD.Rmd
+ Table 2: /02_scripts/02_CTD_outputs/03_comparison_Willoughby2023.R
+ Fig.3a-d: /02_scripts/03_CTD_results/00_vizualize.R
+ Fig.3g: /02_scripts/03_CTD_results/03_comparison_Willoughby2023.R

### Multivariate CTD model shows shared and distinct genetic influences on inter-individual differences in aesthetic value across visual domains.
+ Text: 03_outputs/report/04_phenotypic_cor.Rmd
+ Text: 03_outputs/report/05_multivariate_CTD.Rmd
+ Fig.4a-d: /02_scripts/03_CTD_results/00_vizualize.R
+ Fig.5c&d: /02_scripts/03_CTD_results/04_multivariate_CTD.R

### Sensitivity analyses discount contributions of potential confounding effects.
+ Text:  06_outputs/report/06_sensitivity_analysis.Rmd

###  S1: Intra-rater test-retest reliability and exclusion criteria
+ Supp.Fig.1: /02_scripts/03_CTD_results/00_vizualize.R

###  S2: Intra-image test-retest reliability 
+ Supp.Fig.2: /02_scripts/0_Review/01_Image_reliability.R

###  S3: Multilevel Modelling and Variance Partitioning Components
+ Supp.Fig.3: /02_scripts/03_CTD_results/00_vizualize.R

###  S4: Taste-typically comparison with pairwise aesthetic agreement
+ Supp.Fig.4: /02_scripts/03_CTD_results/00_vizualize.R

###  S5: Principal Component Analysis 
+ Text:  06_outputs/report/S05_PCA.Rmd
+ Supp.Fig.5a-c: /02_scripts/03_CTD_results/00_vizualize.R

###  S6: Monozygotic twins show higher pairwise aesthetic agreement than dizygotic twins and unrelated pairs: results with outliers and from the validation sample 
+ Text:  06_outputs/report/S06_aesthetic_agreement.Rmd
+ Supp.Fig.6: /02_scripts/03_CTD_results/00_vizualize.R
+ Supp.Fig.7: /02_scripts/03_CTD_results/00_vizualize.R

### S7. Within-day test-retest reliabilities for taste-typicality and evaluation-bias
+ Text:  06_outputs/report/S07_withinday_reliability.Rmd

### S8. Between-day reliabilities for pairwise aesthetic agreement, taste-typicality, and evaluation-bias
+ Text:  06_outputs/report/S08_betweenday_reliability.Rmd

### S9. Phenotypic correlations replicate in a third non-overlapping sample
+ Text:  06_outputs/report/S09_phenotypic_cor_replication.Rmd

### S10: Multivariate Classical Twin Design Sensitivity Analyses 
+ Supp.Fig.8: /02_scripts/03_CTD_results/04_multivariate_CTD.R
+ Supp.Fig.9: /02_scripts/03_CTD_results/04_multivariate_CTD.R

### S11: Power analysis 
+ Supp.Fig.10: /02_scripts/04_Review/04_CTD_power_C.R
+ Supp.Fig.11: /02_scripts/04_Review/04_CTD_power_C.R

### S12: Sample overlap
+ Supp.Fig.12: NA
+ Output: 02_scripts/04_Review/S05_sample_overlap.R

### Table S1: Phenotypic Twin correlations before and after controlling for confounding effects
Table: /02_scripts/03_CTD_results/01_SAT_twin_correlations.R

### Table S2: Univariate modelling of genetic and environmental contributions to inter-individual differences in aesthetic evaluation after controlling for confounding effects
Table: /02_scripts/03_CTD_results/02_CTD_outputs.R

### Supplementary Data 1
+ /03_outputs/processedData/03_CTD_results/02_SF1_SAT_results.csv; script: /02_scripts/03_CTD_results/02_CTD_outputs.R

### Supplementary Data 2
+ /03_outputs/processedData/03_CTD_results/02_SF2_ACE_results.csv; script: /02_scripts/03_CTD_results/02_CTD_outputs.R

### Supplementary Data 3
+ /03_outputs/processedData/03_CTD_results/02_SF3_ACE_ci.csv; script: /02_scripts/03_CTD_results/02_CTD_outputs.R

### Supplementary Data 4
+ /05_figures/00_F1b_VPC_sourceData.csv; script: /02_scripts/03_CTD_results/00_vizualize.R

### Supplementary Data 5
+ /03_outputs/processedData/01_Germine_2015/03_Fig_F2_individual_pairwise_agreement_revised_sourceData.csv; script: /02_scripts/01_Germine_2015/03_pairwise_preferences.R

### Supplementary Data 6
+ /03_outputs/processedData/01_Germine_2015/03_Fig_F2a_summary_pairwise_agreement_revised_sourceData.csv; script: /02_scripts/01_Germine_2015/03_pairwise_preferences.R

### Supplementary Data 7
+ /03_outputs/processedData/02_Sutherland_2020/03_Fig_F2_individual_pairwise_agreement_revised_sourceData.csv; script: /02_scripts/02_Sutherland_2020/03_pairwise_preferences.R

### Supplementary Data 8
+ /03_outputs/processedData/02_Sutherland_2020/03_Fig_F2b_summary_pairwise_agreement_revised_sourceData.csv; script: /02_scripts/02_Sutherland_2020/03_pairwise_preferences.R

+ *Note 1*: Please note that after revision, the supplementary material had to be reformatted. Sections S1, S2, etc.... were split into Supplementary Notes and Supplementary Figures.

+ *Note 2*: Please note that after revision, the order of Supplementary Data has changed. This readme is up to date with the final and official order for the Supplementary Data.

## License 

The license for the current repository is:

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
