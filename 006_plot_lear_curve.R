library(tidyverse)

# This script allows to plot learning curves of RPS and Brier score, for the sign level s model resp. symptom state S models.
# Learning curves need to be previously produced with the following scripts : 
#  - 003_MOLR_model.R for the MOLR model at the sign level s model
#  - 002_LR_model.R for the LR model at the symptom state S model

# We assumed learning curves results (outputs of "eval_MOLR()" and "eval_LR()") functions
# Are stored in a file "Results_storage"

# ------------ PLOT RPS (sign s level ) -------------------

# Load forward chaining evaluation for MOLR, uniform and chance forecast
MOLR <- readRDS("Results_storage/RPS_none.RDS")    # NAME OF THE FILE TO COMPLETE
N_MOLR <- length(MOLR$n_x)
signs <- c("itching", "sleep", "redness", "dryness", "oozing", "edema")
# Color blind palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Stack RPS of different models and signs
plot_MOLR <- tibble(day = c(rep(MOLR$n_x, 6 * 3)),
                    symptoms = rep(unlist(lapply(signs,
                                                 function(sign){
                                                   rep(sign, N_MOLR)
                                                 })), 3)
                    ,
                    Model = factor(c(rep("Our model", 6 * l_MOLR),
                                     rep("Uniform forecast", 6 * l_MOLR),
                                     rep("Historical forecast", 6 * l_MOLR)),
                                   levels = c("Our model", "Uniform forecast", "Historical forecast")),
                    RPS = c(MOLR$RPS_itc, MOLR$RPS_sle, MOLR$RPS_red, MOLR$RPS_dry, MOLR$RPS_ooz,MOLR$RPS_ede, 
                            MOLR$RPS_chance_itc, MOLR$RPS_chance_sle, MOLR$RPS_chance_red, MOLR$RPS_chance_dry, MOLR$RPS_chance_ooz, MOLR$RPS_chance_ede,
                            MOLR$RPS_hist_itc, MOLR$RPS_hist_sle, MOLR$RPS_hist_red, MOLR$RPS_hist_dry, MOLR$RPS_hist_ooz, MOLR$RPS_hist_ede)
)

#Plot RPS learning curves
ggplot(data = plot_MOLR) +
  geom_smooth(aes(x = day, y = RPS, color = Model, fill = Model),
              alpha = 0.2, linetype = 'solid', method = 'loess', span = 0.5) +
  facet_wrap(~ symptoms) +
  labs( x = 'day', y = "RPS")+
  theme_bw() + 
  theme(axis.title.y = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11), 
        axis.title.x = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.position = "top") +
  labs(fill = "", color = "") +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette)

# ------------ PLOT Brier score (AD symptom state S level ) -------------------

MOLR <- readRDS("Results_storage/MOLR.RDS")
LR <- readRDS("Results_storage/LR.RDS")
N_MOLR <- length(MOLR$n_x)
N_LR <- length(LR$n_x)

BS = tibble(day = c(rep(LR$n_x, 3), MOLR$n_x),
            Model = factor(c(rep("Uniform forecast", N_LR),
                             rep("Historical forecast", l_GLM),
                             rep("Kim et al.(2017)", N_LR),
                             rep("Our model", N_MOLR)),
                           levels = c("Our model", "Kim et al.(2017)","Uniform forecast", "Historical forecast")),
            Value = c(rep(0.25, N_LR),
                      LR$Brier_hist,
                      LR$Brier, 
                      MOLR$Brier))

ggplot(BS, aes(x = day, y = Value, color = Model)) +
  geom_smooth(aes(fill = Model), alpha = 0.2, linetype = 'solid', method = 'loess', se = TRUE,span = 0.5) +
  labs(color = "",
       fill = "",
       x = 'day',
       y = "Brier Score") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11), 
        axis.title.x = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.position = "top")+
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette)
