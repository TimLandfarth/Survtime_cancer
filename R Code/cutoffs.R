# 0. Data, libraries and functions----
## 0.1 Loading data----
#### imputed Data
load(file = "./Data/imputed dataset.Rda")

#### Models
load(file = "./Data/models.Rda")

## 0.2 Functions----
source("./R Code/functions.R")

## 0.3 Libraries----
source("./R Code/libraries.R")

## 0.4 Graphical values----
source("./R Code/graphic_parameters.R")

# 1. Determine the optimal cutpoint of linear predictor----
## 1.1 AIC model----
#### prediction of the cox_aic model and implementation into the imputated dataframe
df_imp$predict_lp_aic <- predict(cox_aic, type = "lp")

#### Ordering the imputated dataframe on basis of the predicted values
df_imp <- df_imp[order(df_imp$predict_lp_aic),]

#### Generation of the parameter for determining the two groups
df_imp$grp_aic <- NA

#### Generation of the parameter for determining the maximum logrank statistics
logrank_aic <- NA

#### Display of how far the loop has been calculated
pb = txtProgressBar(min = 1, max = length(unique(df_imp$predict_lp_aic)), initial = 1)

#### Setting seed
set.seed(12345)

#### For all observations in the imputed data set
for(i in 1:length(unique(df_imp$predict_lp_aic))){
  
  #### Display of how far the loop has been calculated
  setTxtProgressBar(pb,i)
  
  #### All observations are first sorted into one group and then the first i 
  #### observations are sorted according to the linear predictor
  df_imp$grp_aic <- 1
  df_imp$grp_aic[which(df_imp$predict_lp_aic <= unique(df_imp$predict_lp_aic)[i])] <- 0
  df_imp$grp_aic <- factor(df_imp$grp_aic)
  
  #### Calculation of the logrank (chi squared) for the two groups.
  #### The covariates are not considered further here, as they have no influence
  #### on the logrank statistics.
  logrank_aic[i] <- survdiff(Surv(time, as.numeric(status)-1) ~ grp_aic, data = df_imp)$chisq
}

close(pb)
rm(i)

#### Maximum logrank is the following linear predictor
unique(df_imp$predict_lp_aic)[which(logrank_aic == max(logrank_aic))]

#### With the maximum logrank
max(logrank_aic)

#### Selecting the groups that had the maximum logrank
df_imp$grp_aic <- 1
df_imp$grp_aic[df_imp$predict_lp_aic <= unique(df_imp$predict_lp_aic)[which(logrank_aic == max(logrank_aic))]] <- 0

p1 <- survminer::ggsurvplot(survfit(Surv(time, as.numeric(status)-1) ~ grp_aic, data = df_imp),
                      xlim = c(0,2000))$plot+
  theme_tufte()+
  theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
        legend.position = "bottom")+
  scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
  xlab("Time [days]")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))
  
## 1.2 for the bic model----
#### prediction of the cox_aic model and implementation into the imputated dataframe
df_imp$predict_lp_bic <- predict(cox_bic, type = "lp")

#### Ordering the imputated dataframe on basis of the predicted values
df_imp <- df_imp[order(df_imp$predict_lp_bic),]

#### Generation of the parameter for determining the two groups
df_imp$grp_bic <- NA

#### Generation of the parameter for determining the maximum logrank statistics
logrank_bic <- NA

#### Display of how far the loop has been calculated
pb = txtProgressBar(min = 1, max = length(unique(df_imp$predict_lp_bic)), initial = 1)

#### Setting seed
set.seed(12345)

#### For all observations in the imputed data set
for(i in 1:length(unique(df_imp$predict_lp_bic))){
  
  setTxtProgressBar(pb,i)
  
  #### All observations are first sorted into one group and then the first i 
  #### observations sorted according to the linear predictor
  df_imp$grp_bic <- 1
  df_imp$grp_bic[which(df_imp$predict_lp_bic <= unique(df_imp$predict_lp_bic)[i])] <- 0
  df_imp$grp_bic <- factor(df_imp$grp_bic)
  
  #### Calculation of the logrank (chi squared) for the two groups.
  #### The covariates are not considered further here, as they have no influence
  #### on the logrank statistics.
  logrank_bic[i] <- survdiff(Surv(time, as.numeric(status)-1) ~ grp_bic, data = df_imp)$chisq
}

close(pb)
rm(i)

#### Maximum logrank is the following linear predictor
unique(df_imp$predict_lp_bic)[which(logrank_bic == max(logrank_bic))]

#### With the maximum logrank
max(logrank_bic)

#### Selecting the groups that had the maximum logrank
df_imp$grp_bic <- 1
df_imp$grp_bic[df_imp$predict_lp_bic <= unique(df_imp$predict_lp_bic)[which(logrank_bic == max(logrank_bic))]] <- 0

p2 <- survminer::ggsurvplot(survfit(Surv(time, as.numeric(status)-1) ~ grp_bic, data = df_imp),
                      xlim = c(0,2000))$plot+
  theme_tufte()+
  theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
        legend.position = "bottom",
        axis.title.y = element_blank())+
  scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
  xlab("Time [days]")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))

cutoff_plot <- cowplot::plot_grid(p1,p2, align = "hv", axis = "tblr", labels = c("AIC", "BIC"), hjust = -3, vjust = 2)
max_cutoff_bic <- unique(df_imp$predict_lp_bic)[which(logrank_bic == max(logrank_bic))]
max_cutoff_aic <- unique(df_imp$predict_lp_aic)[which(logrank_aic == max(logrank_aic))]

#pdf(width = pdf_w_h[1]*2, height = pdf_w_h[2] ,file = "./../Plots/Cutoffs.pdf")
cutoff_plot
#dev.off()

save(max_cutoff_aic, max_cutoff_bic, cutoff_plot, file = "./../Data/cutoff.Rda")
