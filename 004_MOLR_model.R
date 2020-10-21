library(lme4)
library(tidyverse)
library(foreach)
library(doParallel)
library(AUC)
source("001_data_import.R")


df <- import_data_model(model = "MOLR")

eval_MOLR <- function(df, step_ahead) {
  
  # Input : - df : dataframe computed by df_ordinal_models
  #         - step_ahead : number of steap ahead in the forward chaining setting 
  
  
  # Output : RPS of the MOLR, the chance and historical forecast for the 6 symptoms and for all days of the forward chaining
  
  max_length = max(df$day) - step_ahead 
  #Make clusters for parallel computing
  n_cores = 25
  registerDoParallel(n_cores)
  #Loop over all the training sizes
  metrics <- foreach (n_x = seq(10, max_length, 15), .combine = rbind) %dopar% {
    
    #Set the training and test set
    #Ancienne façon 
    training_set  <-  filter(df, day <= n_x)
    test_set <- arrange(filter(df, (n_x < day) & (day <= (n_x + step_ahead))), ID, j)
    #Fit the MOLR
    model_itching <- glmer(itching_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(itching_t0), data = training_set, family = binomial(link = 'logit'))
    model_sleep <- glmer(sleep_t1_cumu_y ~ factor(j) + (1 | ID) +   factor(sleep_t0) , data = training_set, family = binomial(link = 'logit'))
    model_redness <- glmer(redness_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(redness_t0) , data = training_set, family = binomial(link = 'logit'))
    model_dry <- glmer(dry_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(dry_t0) , data = training_set, family = binomial(link = 'logit'))
    model_oozing <- glmer(oozing_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(oozing_t0), data = training_set, family = binomial(link = 'logit'))
    model_edema <- glmer(edema_t1_cumu_y ~ factor(j) + (1 | ID) + factor(edema_t0) , data = training_set, family = binomial(link = 'logit'))
    
    #Add the predictions to the dataset
    df_metric <- mutate(test_set, ypred_itching = predict(model_itching, newdata = test_set, type = "response"),
                        ypred_sleep = predict(model_sleep, newdata = test_set, type = "response"),
                        ypred_redness = predict(model_redness, newdata = test_set, type = "response"),
                        ypred_dry = predict(model_dry, newdata = test_set, type = "response"),
                        ypred_oozing = predict(model_oozing, newdata = test_set, type = "response"),
                        ypred_edema = predict(model_edema, newdata = test_set, type = "response"))
    
    IDs <- unique(df_metric$ID)
    RPS <- do.call(rbind, 
                   lapply(1:length(IDs),
                          function(N_ID) {
                            df_metric_ID <- arrange(filter(df_metric, ID == IDs[N_ID]), j)
                            df_metric_ID$ypred_chance <- c(0.2, 0.4, 0.6, 0.8)
                            
                            #Compute the predicted probabilities of the historical forecast for this ID and this day
                            df_metric_ID$ypred_hist_itc <- unlist(lapply(0:3, function(x) length(filter(training_set, ID == IDs[N_ID] & itching_t1 <= x)$itching_t1) / length(filter(training_set, ID == IDs[N_ID])$itching_t1)))
                            df_metric_ID$ypred_hist_sle <- unlist(lapply(0:3, function(x) length(filter(training_set, ID == IDs[N_ID] & sleep_t1 <= x)$sleep_t1) / length(filter(training_set, ID == IDs[N_ID])$sleep_t1)))
                            df_metric_ID$ypred_hist_red <- unlist(lapply(0:3, function(x) length(filter(training_set, ID == IDs[N_ID] & redness_t1 <= x)$redness_t1) / length(filter(training_set, ID == IDs[N_ID])$redness_t1)))
                            df_metric_ID$ypred_hist_dry <- unlist(lapply(0:3, function(x) length(filter(training_set, ID == IDs[N_ID] & dry_t1 <= x)$dry_t1) / length(filter(training_set, ID == IDs[N_ID])$dry_t1)))
                            df_metric_ID$ypred_hist_ooz <- unlist(lapply(0:3, function(x) length(filter(training_set, ID == IDs[N_ID] & oozing_t1 <= x)$oozing_t1) / length(filter(training_set, ID == IDs[N_ID])$oozing_t1)))
                            df_metric_ID$ypred_hist_ede <- unlist(lapply(0:3, function(x) length(filter(training_set, ID == IDs[N_ID] & edema_t1 <= x)$edema_t1) / length(filter(training_set, ID == IDs[N_ID])$edema_t1)))
                            
                            return(tibble(n_x = n_x,
                                          RPS_itc =  1/4 * sum((df_metric_ID$ypred_itching - df_metric_ID$itching_t1_cumu_y)^2),
                                          RPS_sle  = 1/4 * sum((df_metric_ID$ypred_sleep - df_metric_ID$sleep_t1_cumu_y)^2),
                                          RPS_red  = 1/4 * sum((df_metric_ID$ypred_redness - df_metric_ID$redness_t1_cumu_y)^2),
                                          RPS_dry  = 1/4 * sum((df_metric_ID$ypred_dry - df_metric_ID$dry_t1_cumu_y)^2),
                                          RPS_ooz  = 1/4 * sum((df_metric_ID$ypred_oozing - df_metric_ID$oozing_t1_cumu_y)^2),
                                          RPS_ede  = 1/4 * sum((df_metric_ID$ypred_edema - df_metric_ID$edema_t1_cumu_y)^2),
                                          
                                          #Compute the RPS of the chance forecast for this ID and this day
                                          RPS_chance_itc  =  1/4 * sum((df_metric_ID$ypred_chance - df_metric_ID$itching_t1_cumu_y)^2),
                                          RPS_chance_sle  =  1/4 * sum((df_metric_ID$ypred_chance - df_metric_ID$sleep_t1_cumu_y)^2),
                                          RPS_chance_red  =  1/4 * sum((df_metric_ID$ypred_chance - df_metric_ID$redness_t1_cumu_y)^2),
                                          RPS_chance_dry  =  1/4 * sum((df_metric_ID$ypred_chance - df_metric_ID$dry_t1_cumu_y)^2),
                                          RPS_chance_ooz  =  1/4 * sum((df_metric_ID$ypred_chance - df_metric_ID$oozing_t1_cumu_y)^2),
                                          RPS_chance_ede  =  1/4 * sum((df_metric_ID$ypred_chance - df_metric_ID$edema_t1_cumu_y)^2),
                                          
                                          #Compute the RPS of the historical forecast for this ID and this day
                                          RPS_hist_itc  =  1/4 * sum((df_metric_ID$ypred_hist_itc - df_metric_ID$itching_t1_cumu_y)^2),
                                          RPS_hist_sle  =  1/4 * sum((df_metric_ID$ypred_hist_sle - df_metric_ID$sleep_t1_cumu_y)^2),
                                          RPS_hist_red  =  1/4 * sum((df_metric_ID$ypred_hist_red - df_metric_ID$redness_t1_cumu_y)^2),
                                          RPS_hist_dry  =  1/4 * sum((df_metric_ID$ypred_hist_dry - df_metric_ID$dry_t1_cumu_y)^2),
                                          RPS_hist_ooz  =  1/4 * sum((df_metric_ID$ypred_hist_ooz - df_metric_ID$oozing_t1_cumu_y)^2),
                                          RPS_hist_ede  =  1/4 * sum((df_metric_ID$ypred_hist_ede - df_metric_ID$edema_t1_cumu_y)^2)))
                          }
                   ))
    # Compute the brier score for the prediction of S_t1 
    brier <- do.call(rbind, 
                     lapply(1:length(IDs),
                            function(N_ID){
                              df_metric_ID <- arrange(filter(df_metric, ID == IDs[N_ID]), j)
                              
                              itching0 <- filter(df_metric_ID, j == 0)$ypred_itching
                              itching1 <- filter(df_metric_ID, j == 1)$ypred_itching
                              sleep0 <- filter(df_metric_ID, j == 0)$ypred_sleep
                              sleep1 <- filter(df_metric_ID, j == 1)$ypred_sleep
                              dry0 <- filter(df_metric_ID, j == 0)$ypred_dry
                              edema0 <- filter(df_metric_ID, j == 0)$ypred_edema
                              oozing0 <- filter(df_metric_ID, j == 0)$ypred_oozing
                              redness0 <- filter(df_metric_ID, j == 0)$ypred_redness
                              
                              p_S = (1 - itching0 * sleep0
                                     - itching0 * (sleep1 - sleep0)
                                     - sleep0 * (itching1 - itching0)) *
                                (1 -  dry0 * edema0 * oozing0 * redness0
                                 -  dry0 * edema0 * redness0 * (1 - oozing0)
                                 -  edema0 * oozing0 * redness0 * (1 - dry0)
                                 -  dry0 * oozing0 * redness0 * (1 - edema0)
                                 -  dry0 * oozing0 * edema0 * (1 - redness0))
                              return(brier = (as.integer(df_metric_ID$S_t1[1]) - p_S)^2)
                            }
                     )
    )
    return(mutate(RPS, brier = brier))
  }
  stopImplicitCluster()
  return(metrics)
}

eval_MOLR_res <- eval_MOLR(df, step_ahead = 1)