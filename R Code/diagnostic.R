# 0. Data, libraries and functions----
## 0.1 Loading data----
#### imputed Data
load(file = "./Data/imputed_dataset.Rda")

#### Models
load(file = "./Data/models.Rda")

## 0.2 Loading libraries----
source("./R Code/libraries.R")

## 0.3 Loading graphic parameters----
source("./R Code/graphic_parameters.R")

# 1. Diagnostics----
## 1.1 Proportional hazard assumption----
### 1.1.1 Test----
#### Testing the proportional hazard assumptions is consistent.
#### All co-variables should have a non-significant test.
cox.zph(cox_aic)
cox.zph(cox_bic)

#### In both cases, except for the monocytes,
#### all covariates are significant (in some cases even strongly significant),
#### so that the proportional hazard assumption may not apply.

### 1.1.2 Schoenfeld Residuals over time for seperate covariates----
#### The Schoenfeld residuals should scatter around 0 over time
#### (covariates are independent of time when proportional hazards are assumed)
#### and should have a uniform variability.

#pdf(width = pdf_w_h[1]*2, height = pdf_w_h[2]*2,file = "./../Plots/schoenfeldres_AICmodel.pdf")
ggcoxzph(cox.zph(cox_aic), point.col = gray3[1])
#dev.off()

#pdf(width = pdf_w_h[1], height = pdf_w_h[2],file = "./../Plots/schoenfeldres_BICmodel.pdf")
ggcoxzph(cox.zph(cox_bic), point.col = gray3[1])
#dev.off()

#### In contrast to the tests from 1.1.1, the coefficients in both models are well distributed,
#### lying at 0 on the one hand and scattering fairly evenly (except for the edges) on the other.

