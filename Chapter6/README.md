# Chapter 6 - Machine Learning for Variant Prioritisation
Jupyter+R: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD)<br />
RStudio: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD?urlpath=urlpath%3Drstudio)
(binder will only work when the repo is made public)


# Contents
[1. Introduction](#introduction)<br />
[2. Step-by-step pipeline guide](#step-by-step-pipeline-guide)<br />
[3. References](#references)<br />

# Introduction

In this chapter, I develop a variant prioritisation framework, which I test on BP GWAS data, investigating whether variant prioritisation that is specific to BP variants from the GWAS by Evangelou et al.[1] is possible using positive-unlabelled learning.


## Step-by-step Pipeline Guide

![overview](https://i.imgur.com/A8zOPtX.png)

**Overview of the Machine Learning Framework** 

Variant prioritisation uses all significantly associated variants BP GWAS variants from Evangelou et al. [1] (n=37,436 in LD > 0.8 r2 with a sentinel SNP) as input. These variants were annotated in the genomic analysis package Hail using seven resources (six inside of Hail: VEP, ClinVar, CADD, DANN, ENCODE, and dbNSFP, and one externally: the UCSC Genome Browser). The ClinVar annotation is then used to identify pathogenic variants to act as positive-labelled examples in machine learning. Positively labelled variants were identified in two groups (those pathogenic for cardiovascular diseases and those pathogenic for any disease). The second stage then cleans the collected annotations before they enter machine learning as features, with two cleaning approaches tested to explore keeping as many annotations as possible. Machine learning then tests positive-unlabelled machine learning approaches with different methods (bagging versus Elkanoto), models (XGBoost versus LightGBM) and open-source packages versus from scratch coding of bagging. Finally, the top performing method is used for variant prioritisation of all associated variants to then assess the scoring in comparison to the gene-level prioritisation framework developed in chapter 4.

# Code

All steps of the framework can be ran using the ```Hail_BP_Variant_Prioritisation.ipynb``` script [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter6/Hail_BP_Variant_Prioritisation.ipynb). This code includes all data pre-processing steps, training gene labelling, and machine learning application.

Additional re-runs of the framework without Hail were ran using the notebooks for prioritising all cardiovascular positively labelled variants [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter6/BPGWAS%20Variant%20Prioritisation%20-%20all%20betas%20FS.ipynb) versus using positively labelled variants that have ClinVar annotations to any pathogenic variants regardless of disease found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter6/All_Pathogenic_GWAS_Variant_Prioritisation%20-%20all%20betas.ipynb)



# References
1. Evangelou, E., Warren, H.R., Mosen-Ansorena, D. et al. Genetic analysis of over 1 million people identifies 535 new loci associated with blood pressure traits. Nat Genet 50, 1412–1425 (2018). https://doi.org/10.1038/s41588-018-0205-x
