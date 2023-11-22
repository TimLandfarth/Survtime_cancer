# Survival time for new cancer therapy
The code and pseudo data are provided here for the reproducibility of the survival time analysis of a new cancer therapy.

For the study, a Cox proportional hazard model was calculated on imputed, stacked data, variables were selected by AIC and BIC and a bound for the linear predictor was calculated by maximum logrank statistics. For validation, further Cox models were estimated by cross-validation and the Cindex and integrated Brier score were calculated. Furthermore, another cross-validation was performed to validate the bound for the linear predictor.

The structure of the data and the code is as follows:

- **Data**:  All required data
  - **cutoff.Rda**: Maximum cutoff of the linear predictor for the AIC (*max_cutoff_aic*) and BIC (*max_cutoff_bic*) model. Furthermore, the corresponding graph of the survival time curves for the two groups found (*cutoff_plot*; plot is included in Plots/Cutoffs.pdf)
  - **imputed dataset.Rda**: (*df_imp*) 20-fold stacked, imputed data set from the pseudo data
  - **models.Rda**: Original Coxph models, selected using AIC (*cox_aic*) and BIC (*cox_bic*)
  - **prep.RData**: prepared dataset from pseudo data (*df*)
  - **pseudo_data.RData**: generated pseudo data set from the original data (*df*)
  - **session_info.Rda**: Version information of R and packages
  - **validation.Rda**: From original data for cross-validation: optimal cutpoint for the AIC (*opt_cut_aic*) and BIC (*opt_cut_bic*) model, survival plots for the cutpoints for the AIC (*surv_val_aic*; plot is included in Plots/Cutoffs_validated_AICmodel.pdf) and BIC (*surv_val_bic*; plot is included in Plots/Cutoffs_validated_BICmodel.pdf) model, Brier and Cindex with respect to the AIC and BIC model (*performance_plot*; plot is included in Plots/Performance_validated_cox.pdf) and the corresponding values (res).
  
- **Plots**
  - **Cutoffs_validated_AICmodel.pdf**: Survival curves for the two strata found from the five-fold cross-validation for the AIC model (Graphic is saved in *Data/validation.Rda* and calculated in *R Code/validation.R*)
  - **Cutoffs_validated_BICmodel.pdf**: Survival curves for the two strata found from the five-fold cross-validation for the BIC model (Graphic is saved in *Data/validation.Rda* and calculated in *R Code/validation.R*)
  - **Cutoffs.pdf**: Survival curves for the two strata found from the original data (Graphic is saved in *cutoffs.Rda* and calculated in *cutoff.R*)
  - **Performance_validated_cox.pdf**: Brier and Cindex for five-fold cross-validation (Graphic is saved in *Data/validation.Rda* and calculated in *R Code/validation.R*)
  - **schoenfeldres_AICmodel.pdf**: Schoenfeld residuals for the original AIC model (Graphic is calculated in *R Code/ diagnostic.R*)
  - **schoenfeldres_BICmodel.pdf**: Schoenfeld residuals for the original BIC model (Graphic is calculated in *R Code/diagnostic.R*)
  
- **R Code**
  - **cutoffs.R**: Required data: *imputed dataset.Rda*, *models.Rda*. Calculation of the cut-off values based on the maximum logrank statistics for the models. 
  - **diagnostic.R**: Required data: *imputed dataset.Rda*, *models.Rda*. Calculation of the schoenfeld residuals and testing the proportional hazard assumption.
  - **functions.R**: Required data: none. Functions for imputing data and index formation for cross-validation.
  - **generation_of_pseudo_data.R**: Required data: original data (not included here). Generation of a pseudo data set from the original data
  - **graphic_parameters.R**: Required data: none. Colors, heights and widths for generating graphics.
  - **libraries.R**: Required data: none. All packages used for fast loading.
  - **model.R**: Required data: *pseudo_data.RData*. Calculation of the model with stepfunction according to AIC and BIC.
  - **preparation.R**: Required data: none. Preparation of the original data (not included here).
  - **validation.R**: Required data: *pseudo_data.RData*. Validation of the model and the cut-off values.
  
  

