# Prognostic score-based methods for estimating center effects based on survival probability: Application to post-kidney transplant survival

Youjin Lee, Peter P. Reese, and Douglas E. Schaubel

## Overview

In evaluating the performance of different facilities or centers on survival outcomes, the standardized mortality ratio (SMR), which compares the observed to expected mortality has been widely used, particularly in the evaluation of kidney transplant centers. Despite its utility, the SMR may exaggerate center effects in settings where survival probability is relatively high. An example is one-year graft survival among U.S. kidney transplant recipients. We propose a novel approach to estimate center effects in terms of differences in survival probability (i.e., each center versus a reference population).  An essential component of the method is a prognostic score weighting technique, which permits accurately evaluating centers without necessarily specifying a correct survival model. Advantages of our approach over existing facility-profiling methods include a metric based on survival probability (greater clinical relevance than ratios of counts/rates); direct standardization (valid to compare between centers, unlike indirect standardization based methods, such as the SMR); and less reliance on correct model specification (since the assumed model is used to generate risk classes as opposed to fitted-value based `expected' counts). We establish the asymptotic properties of the proposed weighted estimator and evaluate its finite-sample performance under a diverse set of simulation settings. The method is then applied to evaluate U.S. kidney transplant centers with respect to graft survival probability.

## Data

For the real analysis, we used the data obtained from the Scientific Registry of Transplant Recipients (SRTR). Data access can be requested [here](https://www.srtr.org/requesting-srtr-data/data-requests). 

We provide the sample data `Data/sample.csv`, which has the same structure as the real data presented in the paper. The data is provided for illustrative purpose only, not for reproducing the results. 

## Code

- `Code/aux.R` : provides auxiliary functions that estimate the variance of the weighted cumulative hazards and derive prognostic scores from a stratified additive hazards model.

- `Code/exploratory_figure.R`: This code can be used to reproduce Figure 1 in the main manuscript. 

- `Code/main_sim.R`: generates the simulated data and estimates the center effects using prognostic scores. This code can be used to reproduce Tables 1 and 2 in the main manuscript and Table S1 in the Supporting Information. 

- `Code/main_realdata.R`: reads the real data example and estimates the center effects using the proposed methods. This code can be used to reproduce Figures 2-3 in the main manuscript and Figure S1 in the Supporting Information.

- `Code/supp_realdata.R`: explores the impact of varied number of risk classes on the real data application results and compares the estimated excess three-year survival probabilities to the one-year survival probabilities. This code can be used to reproduce Figures S2-S4 and Table S6 in the Supporting Information. 

- `Code/SectionS3.2.R`: investigates the sensitivity to the prognostic score model misspecification. This code can be used to reproduce Tables S2 and S3 in the Supporting Information. 

- `Code/SectionS3.3.R`: performs the sensitivity analysis against the no effect modification assumption. This code can be used to reproduce Table S4 in the Supporting Information. 

- `Code/SectionS3.4.R`: investigates the sensitivity to the violation of positivity assumption. This code can be used to reproduce Table S5 in the Supporting Information. 


## Instructions for the use of the sample data

```{r}
dat = read.csv("Data/sample.csv", header = TRUE, sep = ",")
dim(dat) 
names(dat)
summary(dat$folltime) # distribution of the observed survival times
table(dat$GF) # failure indicator
dat$REC_CTR_ID = as.character(dat$REC_CTR_ID)
table(dat$REC_CTR_ID)
```

Run `Code/supp_realdata.R` file with this `dat`. 

