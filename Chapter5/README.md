# Chapter 5 - Re-application to Blood Lipid Trait GWAS
Jupyter+R: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD)<br />
RStudio: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD?urlpath=urlpath%3Drstudio)


# Contents
[1. Introduction](#introduction)<br />
[2. Installation](#installation)<br />
[3. Step-by-step pipeline guide](#step-by-step-pipeline-guide)<br />
[4. References](#references)<br />


# Introduction
In this chapter, I investigate using blood lipid GWAS data in a re-application of the ML framework developed in chapters 2-4 and explore output prioritised genes and their potential roles in blood lipid biology. Utilising one study that recently genotyped 1.6 million individuals for five blood lipid traits[1]: high-density lipoprotein (HDL), non-high-density lipoprotein (nonHDL) low-density lipoprotein (LDL), total cholesterol (TC), and triglycerides (TG). This GWAS was followed by two further functional genomic analyses[2, 3], prioritising blood lipid trait genes that could serve as opportune examples in training data for ML prioritisation, which is explored here. 

## Step-by-step Pipeline Guide

![overview](https://i.imgur.com/AEN6Sun.png)

**Overview of the Machine Learning Framework** 

GWAS variants from Graham et al. were annotated to genes and evaluated by machine learning. Genes were scored according to the their likely associations to blood lipid traits based on previously mostly curated genes and a group of curated least likely influential genes.

their BP-drug relationship (genes interacting with BP drugs scored at 1), publication significance to BP (genes significantly named in BP publications scored at 0.75) and protein-protein interactions with known BP genes (genes least likely to affect BP with no interactions with BP genes scored at 0.1). Associated genes alongside insignificant least likely BP genes (p-value <0.15, not in 500kb+/- loci, no linkage disequilibrium) were annotated. Biological/functional data was collected from 20 databases. Further gene filtering identified genes that could be use as training data in the machine learning stage. Machine learning was applied benchmarking 14 models using 8 selected features. The top performing trained model was then used to score genes not in the training data, with the genes and their corresponding scores being then sorted into their loci to select the best genes per locus for gene enrichment analysis

### 1. Data preprocessing and integration:

Feature collection followed the same code and steps applied in chapters 2-4.

#### Least likely BP gene filtering

Genes identified as least likely to affect blood lipid traits by p-value and linkage disequilibrium across all five blood lipid trait GWAS'.



#### Training data curation
The final integrated file was then divided into training data and genes to be predicted by the trained model. This was done by identifying three gene groups that would make up the training data:
1. Gold standard genes were defined by causing Mendelian dyslipidaemias.
2. Silver standard genes were defined by genes that had mouse model knockouts with lipid phenotypes.
3. Least likely genes (genes fitlered by distance to blood lipid trait loci, high GWAS p-value, and no linkage disequilibrum r2 measured across all five blood lipid traits' GWAS').



### 2. Machine learning:

#### Feature pre-processing

Features in the training data were removed if there were >25% missing or had >0.9 correlation in the ```EDA-training-data.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter5/EDA-training-data.ipynb)

Feature imputation was performed using random forest imputation followed by feature selection using BorutaShap in the ```FS_figure.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter5/FS%20figure.ipynb)

#### Machine learning benchmarking
 
Regression classification was used with all models benchmarked using repeated nested 5-fold nested cross validation, with each model undergoing Bayesian hyperparamter tuning (code for benchmarking can be found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/Selected-Features-ML.ipynb)).  Median model performances across the repeated folds were then taken to assess model performance. The top performaning model using all features had its best model parameters applied in BorutaShap feature selection to select features.

After feature selection, the models were then benchmarked against each other in repeated nested cross-validation again, with the top performing model used for further analysis. Code for the final model benchmarking that selected the top performing model can be found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter5/FSML_pipeline.py) with the output in ``` fsBenchmarking.sh.o2551209 ```. Due to the size of the training data and the computational time to benchmarking the majority of the models, these results were not written via the ipynb jupyter notebook.

The models benchmarked (extreme gradient boosting, catboost, lightgbm gradient boosting, random forest, extratrees, decision tree, k-nearest neighbors, LASSO and elasticnet logistic regrssion) were further compared with voting and stacking regressors using all other models and a bagging regressor using extreme gradient boosting within the ```Meta-estimators benchmarking.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter5/Meta-estimators%20benchmarking.ipynb).

The trained top performing model was given non-training data genes to predict their likelihood of affecting BP within the ```XGB-Predictions.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter5/XGB-Predictions.ipynb).

The top performing model was further anlysed using SHAP in the ```XGB-SHAP-interpretation.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter5/XGB-SHAP-interpretation.ipynb).
 

# References
1. Graham SE, Clarke SL, Wu K-HH, Kanoni S, Zajac GJM, Ramdas S, et al. The power of genetic diversity in genome-wide association studies of lipids. Nature. 2021;600(7890):675-9.
2. Ramdas S, Judd J, Graham SE, Kanoni S, Wang Y, Surakka I, et al. A multi-layer functional genomic analysis to understand noncoding genetic variation in lipids. bioRxiv. 2021:2021.12.07.470215.
3. Kanoni S, Graham SE, Wang Y, Surakka I, Ramdas S, Zhu X, et al. Implicating genes, pleiotropy and sexual dimorphism at blood lipid loci through multi-ancestry meta-analysis. medRxiv. 2021:2021.12.15.21267852

