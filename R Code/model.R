# 0. Import data, libraries and functions----
## 0.1 Loading data----
load(file = "./../Data/prep.RData")

## 0.2 Loading functions----
source("./functions.R")

## 0.3 libraries----
source("./libraries.R")

# 1. Imputation of the data----
# simultaneous imputation and rephrasing from list to dataframe
set.seed(12345)
df_imp <- do.call("rbind", imputedata(df, m = 20))

# 2. Model
## 2.1 Calculate model without covariates and saturated----
cox_without <- coxph(Surv(time, as.numeric(status)-1) ~ 1, data = df_imp)
formula <- "Surv(time, as.numeric(status)-1) ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein +
ldh + mg"
cox_saturated <- coxph(as.formula(formula), data = df_imp, ties = "efron", singular.ok = F)

## 2.2 Calculating step AIC----
cox_aic <- step(cox_without, direction="forward",
                scope=list(lower=cox_without, upper=cox_saturated),
                k=2*length(df_imp)/mean(!is.na(df[,!c(names(df) %in% c("status", "time"))])))

## 2.3 Calculating step BIC----
cox_bic <- step(cox_without, direction="forward",
                scope = list(lower = cox_without, upper = cox_saturated),
                k=log(nrow(df_imp)) * length(df_imp) / mean(!is.na(df[,!c(names(df) %in% c("status", "time"))])))

# 3. Saving ----
## 3.1 Models----
save(cox_aic, cox_bic, file = "./../Data/models.Rda")

## 3.2 Imputed dataset----
save(df_imp, file = "./../Data/imputed dataset.Rda")
