# 0. Data, libraries and functions----
## 0.1 Loading data----
#load(file = "./Data/prep.RData")
load(file = "./Data/pseudo_data.RData")

## 0.2 Loading functions----
source("./R Code/functions.R")

## 0.3 Loading libraries----
source("./R Code/libraries.R")

## 0.4 Loading graphical values----
source("./R Code/graphic_parameters.R")

# 1. Validation of Model----
#### Formula for the saturated model----
formula <- "Surv(time, as.numeric(status)-1) ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein +
ldh + mg"

#### Crossvalidationfolds
K = 5

#### Number of multiple imputations
m = 20

#### Variable generation for models from the BIC and AIC 
#### and the performance measures Cindex and Integrated Brier Score
cindex_model.aic <- cindex_model.bic <- brier_model.aic <- brier_model.bic <- NA
set.seed(12345)

#### Loop for the number of repetitions for the entire calculation
for(j in 1:10) {
  #### Variable generation for temporary variables for storing the Cindex,
  #### Ibrier for the models from AIC and BIC in each loop 
  cindex_model.aic._temp <- cindex_model.bic._temp <- NA
  brier_model.aic._temp <- brier_model.bic._temp <- NA
  
  #### While loop, so that there are always values for the individual scores.
  while (any(is.na(brier_model.aic._temp)) |
         any(is.na(brier_model.bic._temp)) |
         any(is.na(cindex_model.aic._temp)) |
         any(is.na(cindex_model.bic._temp))) {
    
    #### Generation of indexes for K-fold cross-validation
    foldind <- makeCVdiv(data = df,
                         yname = "status",
                         ncv = K)
    
    #### Small display to know which iteration you are in
    print(paste("outer loop = ", j*100/10, "%"))
    
    ## Display of how far the inner loop has been calculated
    pb = txtProgressBar(min = 1, max = K, initial = 1)
    
    #### Cross-validation loop
    for (k in 1:K) {
      #### m multiple imputation of the data set for the test data and the training data
      dataimputedtest <- imputedata(data = df[foldind == k, ], m = m)
      dataimputedtrain <-
        do.call("rbind", imputedata(data = df[foldind != k, ], m = m))
      
      #### Calculation of the model without covariates from the training data set
      modcoxph_without <-
        coxph(Surv(time, as.numeric(status) - 1) ~ 1,
              data = dataimputedtrain,
              x = TRUE)
      
      #### Calculation of the saturated model with all covariates from the training data set
      modcoxph_saturated <-
        coxph(as.formula(formula), data = dataimputedtrain, x = TRUE)
      
      #### Calculation of the forward AIC with corrected loglik (k)
      modcoxph_aic <-
        step(
          modcoxph_without,
          direction = "forward",
          scope = list(lower = modcoxph_without, upper = modcoxph_saturated),
          k = 2 * length(dataimputedtrain) / mean(!is.na(dataimputedtrain[, !c(names(df) %in% c("status", "time"))])) ,
          trace = FALSE
        )
      
      #### Set the performance measures to zero 
      cindexs <- briers <- 0
      
      #### Calculation of the Cindexes and the Brierscores for the calculated AIC model,
      #### for an imputed data set individually, across all imputed data sets.
      for (i in 1:length(dataimputedtest)) {
        cindexs[i] <-
          SurvMetrics::Cindex(modcoxph_aic, predicted = dataimputedtest[[i]])
        briers[i] <- pec::crps(pec::pec(modcoxph_aic, formula = modcoxph_aic$formula, data = dataimputedtest[[i]]))[2]
      }

      #### Arithmetic mean of the calculated Cindexes and Briers for the AIC model
      cindex_model.aic._temp[k] <- mean(cindexs)
      brier_model.aic._temp[k] <- mean(briers)
      
      #### Calculation of the forward BIC with corrected loglik (k)
      modcoxph_bic <-
        step(
          modcoxph_without,
          direction = "forward",
          scope = list(lower = modcoxph_without, upper = modcoxph_saturated),
          k = log(nrow(dataimputedtrain)) * length(dataimputedtrain) / mean(!is.na(dataimputedtrain[, !c(names(dataimputedtrain) %in% c("status", "time"))])),
          trace = FALSE
        )
      
      #### Calculation of the Cindexes and the Brierscores for the calculated AIC model,
      #### for an imputed data set individually, across all imputed data sets.
      for (i in 1:length(dataimputedtest)) {
        cindexs[i] <-
          SurvMetrics::Cindex(modcoxph_bic, predicted = dataimputedtest[[i]])
        briers[i] <- pec::crps(pec::pec(modcoxph_bic, formula = modcoxph_bic$formula, data = dataimputedtest[[i]]))[2]
      }
      
      #### Arithmetic mean of the calculated Cindexes and Briers for the BIC model
      cindex_model.bic._temp[k] <- mean(cindexs)
      brier_model.bic._temp[k] <- mean(briers)
      setTxtProgressBar(pb,k)
      
    }
    close(pb)
  }
  
  #### Combining the various individual performances for the two different models
  cindex_model.aic <- c(cindex_model.aic, cindex_model.aic._temp)
  cindex_model.bic <- c(cindex_model.bic, cindex_model.bic._temp)
  brier_model.bic <-  c(brier_model.bic, brier_model.bic._temp)
  brier_model.aic <- c(brier_model.aic, brier_model.aic._temp)
}

#### Combining the performances into a single dataframe
res <- data.frame(cindex_model.aic = cindex_model.aic[-1], cindex_model.bic = cindex_model.bic[-1], brier_model.bic = brier_model.bic[-1],
            brier_model.aic = brier_model.aic[-1])

#### Recombination of the dataframe for a plot  
res_plot <- gather(res, cond, measurement, cindex_model.aic:brier_model.aic, factor_key=TRUE)
res_plot$cond <- as.character(res_plot$cond)
res_plot_c <- res_plot[res_plot$cond %in% c("cindex_model.aic", "cindex_model.bic"),]
res_plot_b <- res_plot[res_plot$cond %in% c("brier_model.aic", "brier_model.bic"),]

#### Boxplot for the different Models and Performances
p1 <- ggplot(res_plot_c, aes(x = cond, y = measurement, fill = cond))+
  geom_boxplot(show.legend = FALSE)+
  theme_tufte()+
  theme(axis.line = element_line(color = gray2[1], linewidth = axissize),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_x_discrete(labels = c("AIC model", "BIC model"))+
  scale_fill_manual(values = c("cindex_model.aic" = col_s[2], "cindex_model.bic" = col_s[4]))+
  xlab("")+
  ylab("cindex")
p2 <- ggplot(res_plot_b, aes(x = cond, y = measurement, fill = cond))+
  geom_boxplot(show.legend = FALSE)+
  theme_tufte()+
  theme(axis.line = element_line(color = gray2[1], linewidth = axissize),
        axis.ticks.x = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_x_discrete(labels = c("AIC model", "BIC model"))+
  scale_fill_manual(values = c("brier_model.aic" = col_s[2], "brier_model.bic" = col_s[4]))+
  xlab("")+
  ylab("integrated brier score")
performance_plot <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")  

#pdf(width = pdf_w_h[1], height = pdf_w_h[2]*1.5, file = "./Plots/Performance_validated_cox.pdf")
performance_plot
#dev.off()
 
#### summary of the ibrier for the bic model 
summary(res$brier_model.bic)

#### summary of the ibrier for the aic model
summary(res$brier_model.aic)

#### summary of the cindex for the bic model
summary(res$cindex_model.bic)

#### summary of the cindex for the aic model
summary(res$cindex_model.aic)

# 2. Validation of the cutoff values----
#### Formula for the saturated model----
formula <- "Surv(time, as.numeric(status)-1) ~ leuc + lymph + n_gran + mono + e_gran + crp + albumin + protein +
ldh + mg"

#### Crossvalidationfolds
K = 5

#### Number of multiple imputations
m = 20
cutpoints_BIC <- cutpoints_AIC <- list()
plot_aic <- list()
plot_bic <- list()
opt_cut_aic <- NA
opt_cut_bic <- NA
logrank_aic <- NA
dat_imp_test_aic <- dat_imp_test_bic <- list()

#### Setting seed
set.seed(12345)

#### Generation of indexes for K-fold cross-validation
foldind <- makeCVdiv(data = df,
                     yname = "status",
                     ncv = K)

#### Cross-validation loop
for (k in 1:K) {
  #### Small display to know which iteration you are in
  print(paste("loop = ", k*100/K, "%"))
  
  ## Generation of the parameter for determining the two groups
  dataimputedtest <- df %>% mutate(grp_aic = NA)
  dataimputedtrain <- df %>% mutate(grp_aic = NA)
  dataimputedtest <- df %>% mutate(grp_bic = NA)
  dataimputedtrain <- df %>% mutate(grp_bic = NA)
  
  while (length(unique(dataimputedtest$grp_aic)) == 1 | length(unique(dataimputedtest$grp_bic)) == 1) {
    #### m multiple imputation of the data set for the test data and the training data
    dataimputedtest <-
      do.call("rbind", imputedata(data = df[foldind == k,], m = 1))
    dataimputedtrain <-
      do.call("rbind", imputedata(data = df[foldind != k,], m = m))
    
  #### Calculation of the model without covariates from the training data set
  modcoxph_without <-
    coxph(Surv(time, as.numeric(status) - 1) ~ 1,
          data = dataimputedtrain,
          x = TRUE)
  
  #### Calculation of the saturated model with all covariates from the training data set
  modcoxph_saturated <-
    coxph(as.formula(formula), data = dataimputedtrain, x = TRUE)
  
  #### Calculation of the forward AIC with corrected loglik (k)
  modcoxph_aic <-
    step(
      modcoxph_without,
      direction = "forward",
      scope = list(lower = modcoxph_without, upper = modcoxph_saturated),
      k = 2 * length(dataimputedtrain) / mean(!is.na(dataimputedtrain[,!c(names(dataimputedtrain) %in% c("status", "time"))])) ,
      trace = FALSE
    )
  
  #### Predicting the model according to the AIC model
  dataimputedtrain$predict_lp_aic <-
    predict(modcoxph_aic, type = "lp")
  
  #### Ordering the imputated dataframe on basis of the predicted values
  dataimputedtrain <-
    dataimputedtrain[order(dataimputedtrain$predict_lp_aic), ]
  
  #### Generation of the parameter for determining the maximum logrank statistics
  logrank_aic <- NA
  pb = txtProgressBar(min = 1, max = length(unique(dataimputedtrain$predict_lp_aic))-1, initial = 1) 
  
  for (i in 1:(length(unique(dataimputedtrain$predict_lp_aic)) - 1)) {
    setTxtProgressBar(pb,i)
    
    #### All observations are first sorted into one group and then the first i
    #### observations are sorted according to the linear predictor
    dataimputedtrain$grp_aic <- 1
    dataimputedtrain$grp_aic[which(dataimputedtrain$predict_lp_aic <= unique(dataimputedtrain$predict_lp_aic)[i])] <-
      0
    dataimputedtrain$grp_aic <- factor(dataimputedtrain$grp_aic)
    
    #### Calculation of the logrank (chi squared) for the two groups.
    #### The covariates are not considered further here, as they have no influence
    #### on the logrank statistics.
    logrank_aic[i] <-
      survdiff(Surv(time, as.numeric(status) - 1) ~ grp_aic, data = dataimputedtrain)$chisq
  }
  close(pb)
  
  #### Selecting the groups that had the maximum logrank
  dataimputedtest$predict_lp_aic <-
    predict(modcoxph_aic, newdata = dataimputedtest, type = "lp")
  
  #### Ordering the imputated dataframe on basis of the predicted values
  dataimputedtest <-
    dataimputedtest[order(dataimputedtest$predict_lp_aic), ]
  dataimputedtest$grp_aic <- 1
  dataimputedtest$grp_aic[dataimputedtest$predict_lp_aic <= unique(dataimputedtrain$predict_lp_aic)[which(logrank_aic == max(logrank_aic))]] <-
    0
    plot_aic[[k]] <-
    survminer::ggsurvplot(survfit(Surv(time, as.numeric(status) - 1) ~ grp_aic, data = dataimputedtest),
                          xlim = c(0,2000), break.x.by = 250)
    opt_cut_aic[k] <-
    unique(dataimputedtrain$predict_lp_aic)[which(logrank_aic == max(logrank_aic))]
    dat_imp_test_aic[[k]] <- dataimputedtest
  
  #### Calculation of the forward BIC with corrected loglik (k)
  modcoxph_bic <-
    step(
      modcoxph_without,
      direction = "forward",
      scope = list(lower = modcoxph_without, upper = modcoxph_saturated),
      k = log(nrow(dataimputedtrain)) * length(dataimputedtrain) / mean(!is.na(dataimputedtrain[,!c(names(dataimputedtrain) %in% c("status", "time"))])),
      trace = FALSE
    )
  dataimputedtrain$predict_lp_bic <-
    predict(modcoxph_bic, type = "lp")
  
  #### Ordering the imputated dataframe on basis of the predicted values
  dataimputedtrain <-
    dataimputedtrain[order(dataimputedtrain$predict_lp_bic), ]
  
  #### Generation of the parameter for determining the two groups
  dataimputedtrain$grp_bic <- NA
  
  #### Generation of the parameter for determining the maximum logrank statistics
  logrank_bic <- NA
  pb = txtProgressBar(min = 1, max = length(unique(dataimputedtrain$predict_lp_bic))-1, initial = 1) 

  for (i in 1:(length(unique(dataimputedtrain$predict_lp_bic)) - 1)) {
    setTxtProgressBar(pb,i)
    
    #### All observations are first sorted into one group and then the first i
    #### observations sorted according to the linear predictor are sorted into the
    #### other group
    dataimputedtrain$grp_bic <- 1
    dataimputedtrain$grp_bic[which(dataimputedtrain$predict_lp_bic <= unique(dataimputedtrain$predict_lp_bic)[i])] <-
      0
    dataimputedtrain$grp_bic <- factor(dataimputedtrain$grp_bic)
    
    #### Calculation of the logrank (chi squared) for the two groups.
    #### The covariates are not considered further here, as they have no influence
    #### on the logrank statistics.
    logrank_bic[i] <-
      survdiff(Surv(time, as.numeric(status) - 1) ~ grp_bic, data = dataimputedtrain)$chisq
  }
  close(pb)
  #### Selecting the groups that had the maximum logrank
  dataimputedtest$predict_lp_bic <-
    predict(modcoxph_bic, newdata = dataimputedtest, type = "lp")
  
  #### Ordering the imputated dataframe on basis of the predicted values
  dataimputedtest <-
    dataimputedtest[order(dataimputedtest$predict_lp_bic), ]
  dataimputedtest$predict_lp_bic
  dataimputedtest$grp_bic <- 1
  dataimputedtest$grp_bic[dataimputedtest$predict_lp_bic <= unique(dataimputedtrain$predict_lp_bic)[which(logrank_bic == max(logrank_bic))]] <-
    0
  
  #### Saving plots, cutpoints and the used dataset
  plot_bic[[k]] <-
    survminer::ggsurvplot(survfit(Surv(time, as.numeric(status) - 1) ~ grp_bic, data = dataimputedtest),
                          xlim = c(0,2000), break.x.by = 250)
  opt_cut_bic[k] <-
    unique(dataimputedtrain$predict_lp_bic)[which(logrank_bic == max(logrank_bic))]
  
  dat_imp_test_bic[[k]] <- dataimputedtest
  }
}


#### Survival plots for the various calculated cutpoints in the AIC case
surv_val_aic <- cowplot::plot_grid(plot_aic[[1]]$plot+
                     theme_tufte()+
                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                           legend.position = "none",
                           axis.title.x = element_blank())+
                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                     xlab("Time [days]")+
                     scale_y_continuous(expand = c(0,0))+
                     scale_x_continuous(expand = c(0,0)),
                   plot_aic[[2]]$plot+
                     theme_tufte()+
                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                           legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.x = element_blank())+
                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                     xlab("Time [days]")+
                     scale_y_continuous(expand = c(0,0))+
                     scale_x_continuous(expand = c(0,0)),
                   plot_aic[[3]]$plot+
                     theme_tufte()+
                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                           legend.position = "none",
                           axis.title.y = element_blank())+
                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                     xlab("Time [days]")+
                     scale_y_continuous(expand = c(0,0))+
                     scale_x_continuous(expand = c(0,0)),
                   plot_aic[[4]]$plot+
                     theme_tufte()+
                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                           legend.position = "none")+
                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                     xlab("Time [days]")+
                     scale_y_continuous(expand = c(0,0))+
                     scale_x_continuous(expand = c(0,0)),
                   plot_aic[[5]]$plot+
                     theme_tufte()+
                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                           legend.position = "bottom",
                           axis.title.y = element_blank())+
                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                     xlab("Time [days]")+
                     scale_y_continuous(expand = c(0,0))+
                     scale_x_continuous(expand = c(0,0)), ncol = 3, nrow = 2, align = "hv", axis = "tblr")

#### Survival plots for the various calculated cutpoints in the BIC case
surv_val_bic <- cowplot::plot_grid(plot_bic[[1]]$plot+
                                     theme_tufte()+
                                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                                           legend.position = "none",
                                           axis.title.x = element_blank())+
                                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                                     xlab("Time [days]")+
                                     scale_y_continuous(expand = c(0,0))+
                                     scale_x_continuous(expand = c(0,0)),
                                   plot_bic[[2]]$plot+
                                     theme_tufte()+
                                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                                           legend.position = "none",
                                           axis.title.y = element_blank(),
                                           axis.title.x = element_blank())+
                                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                                     xlab("Time [days]")+
                                     scale_y_continuous(expand = c(0,0))+
                                     scale_x_continuous(expand = c(0,0)),
                                   plot_bic[[3]]$plot+
                                     theme_tufte()+
                                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                                           legend.position = "none",
                                           axis.title.y = element_blank())+
                                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                                     xlab("Time [days]")+
                                     scale_y_continuous(expand = c(0,0))+
                                     scale_x_continuous(expand = c(0,0)),
                                   plot_bic[[4]]$plot+
                                     theme_tufte()+
                                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                                           legend.position = "none")+
                                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                                     xlab("Time [days]")+
                                     scale_y_continuous(expand = c(0,0))+
                                     scale_x_continuous(expand = c(0,0)),
                                   plot_bic[[5]]$plot+
                                     theme_tufte()+
                                     theme(axis.line = element_line(linewidth  = axissize, color = gray2[2]),
                                           legend.position = "bottom",
                                           axis.title.y = element_blank())+
                                     scale_color_manual(values = c(col_s[2:3]), labels = c("low risk", "high risk"))+
                                     xlab("Time [days]")+
                                     scale_y_continuous(expand = c(0,0))+
                                     scale_x_continuous(expand = c(0,0)), ncol = 3, nrow = 2, align = "hv", axis = "tblr")

#pdf(width = pdf_w_h[1]*1.5, height = pdf_w_h[2] ,file = "./../Plots/Cutoffs_validated_AICmodel.pdf")
surv_val_aic
#dev.off()

#pdf(width = pdf_w_h[1]*1.5, height = pdf_w_h[2] ,file = "./../Plots/Cutoffs_validated_BICmodel.pdf")
surv_val_bic
#dev.off()

save(opt_cut_aic, opt_cut_bic, surv_val_aic, surv_val_bic, performance_plot, res,
     file = "./Data/validation.Rda")
