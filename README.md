# Cortico-limbic structural deviations from normative references in different neuroimmunological disorders

## Table of contents
1. [Introduction](#Introduction)  
2. [Scripts](#Scripts)  
   2.1 [Required R packages and installation](#Check-and-install-packagesR)  
   2.2 [Normative Modeling](#Normative-model-fitR)  
   2.3 [Normative curve estimation and peak age determination](#PeakAge-and-MedianTrajectoriesR)  
   2.4 [Statistical analyses of deviation scores across diseases](#Stastical-analysis-deviations-across-diseasesR)  
   2.5 [Disease Classification (ROC Analysis)]  
       2.5.1 [disease-vs-HC ROC analysis using centile-score-based model](#ROC-disease-vs-HC-QuantDataR)  
       2.5.2 [disease-vs-HC ROC analysis using Raw-data-based model)](#ROC-disease-vs-HC-RawdataR)  
       2.5.3 [One-vs-RestDisease ROC analysis using centile-score-based model)](#ROC-One-vs-RestDisease-QuantDataR)  
       2.5.4 [One-vs-RestDisease ROC analysis using Raw-data-based model)](#ROC-One-vs-RestDisease-RawDataR)  
   2.6 [Clinical scores association analysis](#Partial-correlationR)  
   2.7 [Clinical score prediction](#Clinical-score-predictionR)  
   2.8 [Prognosis risk stratification](#Prognosis-risk-stratificationR)  
3. [License](#License)

   ---

## 1. Introduction

The codes on Chinese normative reference construction of cortico-limbic system and their downstream clinical applications. The main scripts could be found in file “Scripts”. For other detailed scripts on figure plots and comparison analyses, readers could contact the authors: huangchuxin@csu.edu.cn or liuyaou@bjtth.org. Please note that some source functions (#Source-codes) are from Bethlehem, R.A.I., Seidlitz, J., White, S.R. et al. Brain charts for the human lifespan. Nature 604, 525–533 (2022). https://doi.org/10.1038/s41586-022-04554-y.

This GitHub repository contains the main codes required to replicate the Chinese population-specific normative references for the cortico-limbic system and their clinical applications in neuroimmunological diseases.The original datasets are not included in this repository. Due to data privacy and usage agreements, we do not have authorization to redistribute all of the datasets used in the study. However, some datasets can be made available upon reasonable request. Please contact the corresponding author, Prof. Yaou Liu (liuyaou@bjtth.org), for inquiries regarding data access.

---

## 2. Scripts

This repository contains a workflow for cortico-limbic structural normative modeling and clinical applications in neuroimmunological disorders.

### 2.1 Required R packages and installation  

`Check-and-install-packages.R`

Automatically installs all required R packages for normative modeling and statistical analyses.

---

### 2.2 Normative Modeling  

`Normative-model-fit.R`

Fits GAMLSS models for cortico-limbic regions.

---

### 2.3 Normative curve estimation and peak age determination

`PeakAge-and-MedianTrajectories.R`

Derives normative lifespan trajectories and peak ages.

---

### 2.4 Statistical analyses of deviation scores across diseases]

`Stastical-analysis-deviations-across-diseases.R`

Performs group comparisons, effect size estimation, and multiple-comparison correction.

---
 
### 2.5 Disease Classification (ROC Analysis)

`ROC-disease-vs-HC-QuantData.R`
`ROC-disease-vs-HC-Rawdata.R`
`ROC-One-vs-RestDisease-QuantData.R`
`ROC-One-vs-RestDisease-RawData.R`
 
ROC analysis for disease-vs-HC and One-vs-RestDisease classification tasks.The performance was evaluated using two distinct frameworks: one based on centile scores derived from normative modeling, and another based on raw input data.

---
 
### 2.6 Clinical scores association analysis
 
`Partial-correlation.R`

Performs partial correlation analysis to evaluate the association between centile scores and clinical metrics, adjusting for age and sex.

---

### 2.7 Clinical score prediction

`Clinical-score-prediction.R`

Predicts clinical performance using SVR with Elastic Net feature selection, validated via 10-fold cross-validation.
  
### 2.8 Prognosis risk stratification

`Prognosis-risk-stratification.R`

Evaluates the prognostic value of deviation scores via Cox proportional hazards modeling, with risk stratification and survival analysis.

---

## 3. License

**MIT License**

Copyright (c) 2026 huangcx-lab

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
