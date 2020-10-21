library(lme4)
library(tidyverse)
library(foreach)
library(doParallel)
library(AUC)
source("001_data_import.R")

# Import data for LR
df <- import_data_model(model = "LR")

eval_LR <- function(df, step_ahead = 1) {
  
  # This function evaluates the Logistic Regression model of the Ahn et Al.(2017)
  # With a forward chaining method
  
  # Inputs  
  # df : dataframe computed by import_data_model(model = "LR") from "001_data_import.R"
  # step-ahead : number of prediction steps at each step of the chain
  
  # Ouputs : 
  # df containing for each size of training set,  scores averaged over all filtered patients :  
  #         the brier and AUROC scores
  #         for the log reg model and two benchmark models
  
  # Maximum size of the training set
  max_length = max(df$day) - step_ahead 
  # Make clusters for parallel computing
  n_cores = 25
  registerDoParallel(n_cores)
  # Loop over all the training set sizes n_x
  scores <- foreach (n_x = seq(10, max_length , 5), .combine = rbind) %dopar% {
    # train test split
    training_set  <-  filter(df, day < n_x)
    test_set <- filter(df, (n_x <= day) & (day <= (n_x + step_ahead - 1)))[, c(1:6, 13:23)]
    y = test_set$symptom
    # Fit the glmm
    glmm <- glmer(symptom~ temperature  + PM10 +  NO2 + DTR + RH + O3 + SCORAD  + rainfall +  factor(DOW) + factor(fever) + factor(sex) + (1|ID), data = training_set, binomial)
    # Predict with glmm and two benchmark modesl
    y_pred_glm = predict(glmm, newdata = test_set, type = "response")
    y_pred_unif = rep(0.5, length(y_pred_glm))   # Uniform forecast
    y_pred_hist <- unlist(lapply(1:length(y_pred_glm),
                                 function(obs) {
                                   ID_ <- test_set$ID[obs]
                                   return(sum(filter(training_set, ID == ID_)$symptom) / length(filter(training_set, ID == ID_)$symptom))
                                 }))
    # Compute metrics
    auroc_glm = auc(roc(y_pred_glm, factor(y)))
    auroc_unif = auc(roc(y_pred_unif, factor(y)))
    auroc_hist = auc(roc(y_pred_hist, factor(y)))
    metrics = data.frame(training_size = rep(n_x, length(y)),
                         brier = (as.integer(y) - y_pred_glm)^2,
                         brier_unif = (as.integer(y) - y_pred_unif)^2,
                         brier_hist = (as.integer(y) - y_pred_hist)^2,
                         auroc = rep(auroc_glm, length(y)),
                         auroc_unif = rep(auroc_unif, length(y)),
                         auroc_hist = rep(auroc_hist, length(y))
    )
    return(metrics)
  }
  stopImplicitCluster()
  return(as_tibble(scores))
}

# Run LR forward chaining evaluation
eval_LR_res <- eval_LR(df, 1)  