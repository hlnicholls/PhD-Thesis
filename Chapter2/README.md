# Chapter 2 Exploratory Data Analysis
Jupyter+R: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD)<br />
RStudio: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/PhD-Thesis/HEAD?urlpath=urlpath%3Drstudio)


# Contents
[1. Introduction](#introduction)<br />
[2. Installation](#installation)<br />
[3. Data preprocessing and integration](#Data-preprocessing-and-integration)<br />
[4. Exploratory Data Analysis](#Exploratory-Data-Analysis)<br />
[5. References](#references)<br />
[6. Contact](#contact)<br />

# Introduction

This repository explores the integration of a range of multi-omic databases and their potential as features to be used in machine learning, alongside the curation of a subset of training genes from an individual genome-wide association study (GWAS). This exploratory data analysis provides an overview of the training data that is then used in chapters 3 and 4. The data is curated for machine learning to prioritise genes that are most likely to influence blood pressure (BP) traits, with causal genes being defined as those that are likely to contribute to BP - as recognised by machine learning outputs that are trained on curated genes with known BP roles.

# Installation

Data pre-processing was applied in R and machine learning was applied in Python via the Anaconda environment. The setup of the Python Anaconda environment can be installed by running in the Anaconda command prompt:

```conda env create -f MLGWAS_environment.yml```

Another option available is to install required R and Python packages by running ```install.R``` and ```pip install -r requirements.txt```

Gene annotation was applied in unix using ```bedtools``` and ```ANNOVAR```. Installation for bedtools is described here (https://bedtools.readthedocs.io/en/latest/content/installation.html) and installation for ANNOVAR is described here (https://annovar.openbioinformatics.org/en/latest/user-guide/download/)



**Overview of the Machine Learning Framework** 

GWAS variants from Evangelou et al.[1] were annotated to genes and evaluated by machine learning. Genes were scored according to the their likely associations to BP based on their BP-drug relationship (genes interacting with BP drugs scored at 1), publication significance to BP (genes significantly named in BP publications scored at 0.75) and protein-protein interactions with known BP genes (genes least likely to affect BP with no interactions with BP genes scored at 0.1). Associated genes alongside insignificant least likely BP genes (p-value <0.15, not in 500kb+/- loci, no linkage disequilibrium) were annotated. Biological/functional data was collected from 20 databases. Further gene filtering identified genes that could be use as training data in the machine learning stage. Machine learning was applied benchmarking 14 models using 8 selected features. The top performing trained model was then used to score genes not in the training data, with the genes and their corresponding scores being then sorted into their loci to select the best genes per locus for gene enrichment analysis

# Data preprocessing and integration:

An R markdown can be found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter4/Data%20preprocessing/Data-Preprocessing.md) detailing the complete data-processing code steps outlined below.

#### Whole GWAS data preprocessing

Whole GWAS data was taken in order to annotate variant-level data and then take the most significant value per gene. The whole GWAS data was also used to identify variants with high p-values (>0.15) and no linkage disequilibrium measured with BP loci, provide an initial variants that could be identified as least likely to affect BP for further filtering. The further filtering involved selecting only variants not within 500kb+/- BP loci and not closely connected with known BP genes in protein-protein interaction (PPI) networks.

#### Variant-to-Gene annotation

Variants are then annotated to genes using bedtools within unix using a bash script ```Bedtools_Gene_Annotation.bash``` found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter4/Data%20preprocessing/Bedtools_Gene_Annotation.bash)

Bedtools gene annotation involved using the hg19/GRCh37 reference genome from Ensembl (release 92, Homo_sapiens.GRCh37.87) to annotate variants to their closest genes with a cut off of 50kb distance accepted to annotate a variant to that gene.

The annotated genes were further fitered by gene type, selecting only: protein-coding, processed transcripts, pseudogenes, and antisense genes - to ensure genes would no be heavily missing feature data on machine learning.

#### Variant feature annotation and filtering
Variants were annotated in unix using ANNVOAR (2018Apr16). 

The ANNOVAR annotation was completed using ```ANNOVAR_annotation.bash``` found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter4/Data%20preprocessing/ANNOVAR_annotation.bash)

Once ANNOVAR variant annotation was complete, the whole GWAS genes were filtered to only selected associated genes and genes found to be insignificant/least likely to affect BP (p-value >0.15, no linkage disequilibrium, not within 500kb+/- BP loci) - filtering out other genes in the whole GWAS data.

#### Least likely BP gene filtering

Genes identified as least likely to affect BP by p-value and linkage disequilibrium were further filtered by PPI distance to known BP genes, creating a final list of least likely genes to be scored at 0.1 on training.

#### Variant-level features
All variant-level annotations per gene had the gene's most significant variant-level value selected to represent that gene, excluding for beta values in the GWAS summary statistics, which had the maximum absolute value taken across all three BP phenotypes (systolic, diastolic and pulse pressure) for each gene.

#### Gene-level features
GTEx foldchange data was the only gene-level annotation that required further processing before merging into a final file, requiring all GTEx tissues to be merged into one file. All other collected feaures required no further processing and could be directly merged into a final file.

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




# Exploratory Data Analysis

The two training datasets containing either 3 labels or 4 labels both underwent EDA. Features were first removed if they were > 25% missing and then correlation thresholds of 0.85, 0.9 and 0.99 were each tested to remove highly correlating features. The 0.9 correlation threshold was selected due to having the higher ML performance - with EDA using the 0.9 threshold and the 3-label dataset explored in the ```EDA-training-data-3l.ipynb``` script found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter2/BP%20GWAS%20EDA/EDA-training-data-3l.ipynb)

EDA for the 4-label dataset was also explored [here](https://github.com/hlnicholls/PhD-Thesis/tree/main/Chapter2/BP%20GWAS%20EDA/4%20label) with again 0.9 correlation providing the best ML to be the selected feature cleaning method.

The feature selection, which in itself required embedded ML using BorutaShap, is further explored in chapter 3 and can can be found [here](https://github.com/hlnicholls/PhD-Thesis/blob/main/Chapter3/3%20label/correlation09/Feature_Imputation_and_Selection_BorutaShap.ipynb). However, the selected features were explored on EDA presented here in chapter 2 in the ```EDA-training-data-3l.ipynb``` and ```EDA-training-data-4l.ipynb```scripts.

# References
1. Evangelou, E., Warren, H.R., Mosen-Ansorena, D. et al. Genetic analysis of over 1 million people identifies 535 new loci associated with blood pressure traits. Nat Genet 50, 1412–1425 (2018). https://doi.org/10.1038/s41588-018-0205-x

# Contact

To report any issues or for further information, please contact: 

- Hannah Nicholls: hannahnicholls@qmul.ac.uk
