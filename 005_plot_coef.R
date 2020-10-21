library(lme4)
library(tidyverse)
library(foreach)
library(doParallel)
library(AUC)
library(cowplot)
library(ggpubr)
library(latex2exp)
source("001_data_import.R")
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

# ------------ MOLR model ---------------

df <- import_data_model(model = "MOLR")

#Normalize variables we are studying to get standardized coefficients
df$temperature <- (df$temperature - mean(df$temperature)) / sd(df$temperature)
df$rainfall <- (df$rainfall - mean(df$rainfall)) / sd(df$rainfall)
df$RH <- (df$RH - mean(df$RH)) / sd(df$RH)
df$DTR <- (df$DTR - mean(df$DTR)) / sd(df$DTR)
df$PM10 <- (df$PM10 - mean(df$PM10)) / sd(df$PM10)
df$NO2 <- (df$NO2 - mean(df$NO2)) / sd(df$NO2)
df$O3 <- (df$O3 - mean(df$O3)) / sd(df$O3)

itc_cov <- glmer(itching_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(itching_t0) + temperature + rainfall + DTR + RH  + PM10 + NO2 + O3 + TCS + fever + SCORAD , data = df, family = binomial(link = 'logit'))
sle_cov <- glmer(sleep_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(sleep_t0) + temperature + rainfall + DTR + RH  + PM10 + NO2 + O3 + TCS + fever + SCORAD, data = df, family = binomial(link = 'logit'))
red_cov <- glmer(redness_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(redness_t0) +  temperature + rainfall + DTR + RH  + PM10 + NO2 + O3 + TCS + fever + SCORAD, data = df, family = binomial(link = 'logit'))
dry_cov <- glmer(dry_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(dry_t0) +  temperature + rainfall + DTR + RH  + PM10 + NO2 + O3 + TCS + fever + SCORAD, data = df, family = binomial(link = 'logit'))
ooz_cov <- glmer(oozing_t1_cumu_y ~ factor(j) + (1 | ID) +  factor(oozing_t0) + temperature + rainfall + DTR + RH  + PM10 + NO2 + O3 + TCS + fever + SCORAD, data = df, family = binomial(link = 'logit'))
ede_cov <- glmer(edema_t1_cumu_y ~ factor(j) + (1 | ID) + factor(edema_t0) +  temperature + rainfall + DTR + RH  + PM10 + NO2 + O3 + TCS + fever + SCORAD , data = df, family = binomial(link = 'logit'))

# Extract coefficients
MOLR_cov <- tibble(Symptom = c(rep("Itching", 8), rep("Sleep", 8), rep("Redness", 8), rep("Dry", 8), rep("Oozing", 8), rep("Edema", 8)),
                   Variable = factor(rep(c("Temp", "RF", "DTR", "RH", "PM10", "NO2", "O3", "TCS"), 6), levels = c("TCS", "Temp", "RF", "RH", "DTR", "PM10", "O3", "NO2")),
                   Coef = - c(as.vector(summary(itc_cov)$coefficients[9:16, 1]),
                              as.vector(summary(sle_cov)$coefficients[9:16, 1]),
                              as.vector(summary(red_cov)$coefficients[9:16, 1]),
                              as.vector(summary(dry_cov)$coefficients[9:16, 1]),
                              as.vector(summary(ooz_cov)$coefficients[9:16, 1]),
                              as.vector(summary(ede_cov)$coefficients[9:16, 1])),
                   SE = c(as.vector(summary(itc_cov)$coefficients[9:16, 2]),
                          as.vector(summary(sle_cov)$coefficients[9:16, 2]),
                          as.vector(summary(red_cov)$coefficients[9:16, 2]),
                          as.vector(summary(dry_cov)$coefficients[9:16, 2]),
                          as.vector(summary(ooz_cov)$coefficients[9:16, 2]),
                          as.vector(summary(ede_cov)$coefficients[9:16, 2])))

# Builds paths and labels to load validation metrics for the MOLR models with covariates
paths <- paste("Results_storage/RPS_", c("all", "TCS", "T", "RF", "RH", "DTR", "PM10",  "O3", "NO2"), ".RDS", sep = "")
lbls <- c("all", "TCS", "Temp", "RF", "RH", "DTR", "PM10",  "O3", "NO2")
RPS_lvls <- paste0("RPS_", c("itc", "sle", "red", "dry", "ooz", "ede"))
# MOLR model without covariate. We evaluated MM models with covaraites only for days n_x that are multiple of 10 
base <- select(filter(readRDS("Results_storage/RPS_none.RDS"), n_x %% 10 ==0), RPS_lvls)

# Compute pairwise RPS score difference between base model (MOLR without covariates) and augmented models 
# We evaluated MOLR models with covaraites only for days n_x that are multiple of 10 

MOLR_comparaison <- do.call(rbind,
                            lapply(1:length(paths),
                                   function(i) {
                                     df <- select(filter(readRDS(paths[i]), n_x %% 10 == 0 ),  RPS_lvls)
                                     diff <- base - df
                                     mean <- as.double(colMeans(diff))
                                     SE <- apply(diff, 2, sd) / sqrt(length(df$RPS_sle))
                                     return(tibble(sign = c("itching", "sleep", "redness", "dryness", "oozing", "edema"), 
                                                   covariate = factor(rep(lbls[i], 6), levels = lbls), mean = mean, SE = SE))
                                   }))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot() +
  geom_pointrange(data = MOLR_comparaison,
                  mapping = aes(x = covariate, y = mean, ymin = mean - SE, ymax = mean + SE, color = sign), 
                  position = position_dodge(width = 0.5))+
  theme_bw() +
  coord_flip() +
  geom_hline(yintercept = c(0), linetype = "dotted") +
  geom_line(data = tibble(x = factor(c("all", "TCS", "Temp", "RF", "RH","PM10",  "O3", "NO2"),
                                     levels = c("all", "TCS", "Temp", "RF", "RH","PM10",  "O3", "NO2")),
                          y = rep(0, 8)),
            mapping = aes(x = x, y = y)) +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  theme(axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 18),
        legend.text = element_text(size= 12),
        axis.title.x = element_text(size = 16),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.position = "top") +
  labs(y = TeX("$RPS$ - $RPS_{cov}$"), 
       x = "",
       color = "",
       shape = "") +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette)

#Combine with coefficients plot 


p2 <- ggplot() +
  geom_pointrange(data = MOLR_cov,
                  mapping = aes(x = Variable, y = Coef, ymin = Coef - SE, ymax = Coef + SE, color = Symptom), 
                  position = position_dodge(width = 0.5)) +
  theme_bw() +
  coord_flip() +
  geom_hline(yintercept = c(0), linetype = "dotted") +
  geom_line(data = tibble(x = factor(c("TCS", "Temp", "RF", "RH", "DTR", "PM10", "O3", "NO2"),
                                     levels = c("TCS", "Temp", "RF", "RH", "DTR", "PM10", "O3", "NO2")),
                          y = rep(0, 8)),
            mapping = aes(x = x, y = y)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 18),
        legend.text = element_text(size= 12), 
        axis.title.x = element_text(size = 16),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  labs(y = TeX("Standardized coefficient"), 
       x = "",
       color = "",
       shape = "") +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette)

ggarrange(p1, p2, ncol = 2,
          common.legend = TRUE, legend = "top", labels = c("a", "b"))
# ggsave("Plots/MOLR_coef_RPS.jpg", width = 15, height = 7, units = "cm", dpi = 300, scale = 1.5)
