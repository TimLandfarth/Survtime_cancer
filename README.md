# Survival time for new cancer therapy
The code and pseudo data are provided here for the reproducibility of the survival time analysis of a new cancer therapy.

For the study, a Cox proportional hazard model was calculated on imputed, stacked data, variables were selected by AIC and BIC and a bound for the linear predictor was calculated by maximum logrank statistics. For validation, further Cox models were estimated by cross-validation and the Cindex and integrated Brier score were calculated. Furthermore, another cross-validation was performed to validate the bound for the linear predictor.

The structure of the data and the code is as follows:

