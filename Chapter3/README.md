# Chapter 3 - Multiclass Classification to Prioritise Genes Post-GWAS
Jupyter+R: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD)<br />
RStudio: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD?urlpath=urlpath%3Drstudio)
(binder will only work when the repo is made public)

# Contents
[1. Introduction](#introduction)<br />
[2. Step-by-step pipeline guide](#step-by-step-pipeline-guide)<br />


# Introduction

In this chapter, I apply ML benchmarking to the training data curated in chapter 2, posing gene prioritisation as a multiclass supervised learning problem. This approach allows for a thorough assessment of model performance across several metrics, with various methods (curation of training data, models, and class balancing approaches) tested to optimise the ML framework.


## Step-by-step Pipeline Guide

![overview](https://i.imgur.com/UjIebAt.png)

**Overview of the Multiclass Machine Learning Framework** 

GWAS variants from Evangelou et al. were annotated to genes and evaluated by machine learning. Genes were labelled according to the their likely associations to BP based on their BP-drug relationship (genes interacting with BP drugs labelled most likely), publication significance to BP (genes significantly named in BP publications labelled probable), experimental annotation to BP in Ingenuity Pathway Analysis (genes possibly influencing BP labelled possible) and protein-protein interactions with known BP genes (genes least likely to affect BP with no interactions with BP genes labelled least likely). Associated genes alongside insignificant least likely BP genes (p-value <0.15, not in 500kb+/- loci, no linkage disequilibrium) were annotated with biological/functional data was collected from 20 databases. Further gene filtering identified genes that could be use as training data in the machine learning stage. Machine learning was applied benchmarking 14 models using selected features. The top performing trained model was then used to label genes not in the training data. Furthermore, additonal testing was applied by running ML on training data that included all four labels (n=377) or removing the possible gene group (n=293).


#### Training data curation
The final integrated file was then divided into training data and genes to be predicted by the trained model. This was done by identifying  gene groups that would make up the training data.
Firstly a 3-label training dataset (n=293) was defined by:
1. Known BP genes (genes with BP drugs and BP mechanisms)
2. Genes that are probable to affect BP (genes with significance in BP publications in text-mining or had drug interactions known to have a BP side effect)
3. Least likely genes (genes fitlered by distance to BP loci, high GWAS p-value, no linkage disequilibrum r2 measured, and no close PPIs with known BP genes)

Secondly a 4-label training dataset (n=377) was defined by:
1. Known BP genes (genes with BP drugs and BP mechanisms)
2. Genes that are probable to affect BP (genes with significance in BP publications in text-mining or had drug interactions known to have a BP side effect)
3. Genes annotated as having BP study in Ingenuity Pathway Analysis
4. Least likely genes (genes fitlered by distance to BP loci, high GWAS p-value, no linkage disequilibrum r2 measured, and no close PPIs with known BP genes)


### 2. Machine learning:

#### Feature pre-processing

Features in the training data were removed if there were >25% missing or had >0.9 correlation in the ```EDA-training-data.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/EDA-training-data.ipynb)

Feature imputation was performed using random forest imputation followed by feature selection using BorutaShap in the ```All-Features-ML.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/All-Features-ML.ipynb)

#### Machine learning benchmarking
 
For initial machine learning benchmarking, 3-label and 4-label datasets were compared as well as different correlation thresholds before feature selection. 3-label performances can be found [here](https://github.com/hlnicholls/PhD-Thesis/tree/main/Chapter3/3%20label) and 4-label performances can be found [here](https://github.com/hlnicholls/PhD-Thesis/tree/main/Chapter3/4%20label)

Multiclass classification was used with all models benchmarked using repeated nested 5-fold nested cross validation, with each model undergoing Bayesian hyperparamter tuning (code for benchmarking can be found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter3/3%20label/correlation09/Model_benchmarking.ipynb)). Median model performances across the repeated folds were then taken to assess model performance. The top performaning model using all features had its best model parameters applied in BorutaShap feature selection to select features.

After feature selection, the models were then benchmarked against each other in repeated nested cross-validation again, with the top performing model used for further analysis. Final model benchmarking that selected the top performing model can be found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter3/3%20label/correlation09/Feature_Imputation_and_Selection_BorutaShap.ipynb).

The models benchmarked (extreme gradient boosting, catboost, lightgbm gradient boosting, random forest, extratrees, decision tree, k-nearest neighbors, neural network, logistic regrssion) were further compared with voting and stacking using all other models and a bagging model using catboost within the ```Meta-estimators benchmarking.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter3/3%20label/correlation09/Meta-estimators.ipynb).

Re-sampling methods were then explored for the 3-label dataset only, exploring Synthetic Minority Oversampling TEchnique (SMOTE) oversampling and class-weighting - both tested methods can be found [here](https://github.com/hlnicholls/PhD-Thesis/tree/main/Chapter3/3%20label/correlation09).

The trained top performing model on class-weighting was given non-training data genes to predict their likelihood of affecting BP within the ```Classweighting_CB_BPGWASPredict.ipynb.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter3/3%20label/correlation09/Classweighting_CB_BPGWASPredict.ipynb).

The top performing model was further anlysed using SHAP in the ```SHAP_classweighting_CB.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter3/3%20label/correlation09/SHAP_classweighting_CB.ipynb).
 