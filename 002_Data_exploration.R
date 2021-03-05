library(tidyverse)
library(markovchain)
library(corrplot)
source("001_data_import.R")
df <- import_data()

# ----------- histogram of signs --------------------- 

len <- length(df$ID)
concat_symptoms = tibble(type = c(rep('dryness', len), rep('itching', len), rep('sleep', len), rep('redness', len), rep('oozing', len), rep('edema', len)),
                         value = c(df$dry, df$itching, df$sleep, df$redness, df$oozing, df$edema))

stats_symptoms <- concat_symptoms %>%
  group_by(type, value) %>%
  summarize(count = n()) %>%
  group_by(type) %>%
  mutate(density = count / sum(count))

ggplot(data = stats_symptoms,
       aes(x = value, y = density)) +
  geom_col(width = 0.9) +
  facet_wrap(~type) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Sign score ", y = "Frequency") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 20)

if (FALSE) {
  ggsave(file.path("Plots", "hist_signs.jpg"),
         width = 13, height = 8, units = "cm", dpi = 300, scale = 2)
}

# -------------- N(t) : number of patients observations at day t ----------

nt <- df %>%
  group_by(day = as.integer(day)) %>%
  summarize(n_t = n())

ggplot(data = nt)+
  geom_point(aes(day, n_t))+
  labs(x = latex2exp::TeX("$\\textit{t}$ (day)"),
       y = latex2exp::TeX("$\\textit{N(t)}$")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 20, face = "bold"),
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.position = "top")

if (FALSE) {
  ggsave(file.path("Plots", "n(t).jpg"),
         width = 13, height = 8, units = "cm", dpi = 300, scale = 2)
}

# -------------- Patient symptoms trajectories --------------------------

plot_trajectories <- function(df, patient_id) {
  # Plot symptoms trajectories
  #
  # Args:
  # df: Dataframe
  # patient_id: Patient ID
  #
  # Returns:
  # Ggplot
  
  old_lbl <- c("itching", "sleep", "redness", "dry", "oozing", "edema", "symptom")
  lbl <- c("itching", "sleep", "redness", "dryness", "oozing", "edema", "AD state")
  
  stopifnot(is.data.frame(df),
            all(colnames(c("ID", "day", old_lbl) %in% colnames(df))),
            patient_id %in% unique(df[["ID"]]))
  
  plot <- df %>%
    filter(ID == patient_id)
  
  # Add missing values for "blanks" to appear
  day_mis <- setdiff(1:max(plot$day), plot$day)
  if (length(day_mis) > 0) {
    plot <- bind_rows(plot,
                      data.frame(ID = patient_id, day = day_mis)) %>%
      arrange(day)
  }
  
  plot <- plot %>%
    rename(dryness = dry, `AD state` = symptom) %>%
    select(all_of(c("day", lbl))) %>%
    pivot_longer(cols = all_of(lbl), names_to = "symptom", values_to = "score")
  
  p1 <- ggplot(data = filter(plot, symptom != "AD state"),
               aes(day, score)) +
    geom_path(size = 1) +
    geom_point() +
    facet_wrap(~symptom, ncol = 1, strip.position = "right", scales = "free") +
    scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
    labs(y = "") +
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 17),
          strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 9))
  
  p2 <- ggplot(data = filter(plot, symptom == "AD state"),
               aes(day, score)) +
    geom_path(size = 1) +
    geom_point() +
    facet_wrap(~symptom, ncol = 1, strip.position = "right", scales = "free")+
    scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
    labs(x = latex2exp::TeX("$\\textit{t}$ (day)"), y = "") +
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 10),
          axis.title.x =  element_text(size = 17),
          axis.title.y = element_text(size = 17),
          strip.text = element_text(size = 11),
          strip.text.x = element_text(size = 9))
  
  ggpubr::ggarrange(p1, p2, ncol = 1,
                    common.legend = TRUE, legend = "top", heights = c(3, 1))
}

if (FALSE) {
  # 179 for Figure 1
  lapply(c("068", "084", "108", "120", "126", "134", "179", "180", "188", "193", "195"),
         function(x) {
           plot_trajectories(df, paste0("AD-", x))
           ggsave(file.path("Plots", paste0("trajectory_", x, ".jpg")),
                  width = 13, height = 7, units = "cm", dpi = 300, scale = 2.2)
         })
}

# ---------------   Sign correlations --------------

col2 <- colorRampPalette(c('red', 'white', 'blue'))

M <- cor(df[, c('itching', 'sleep', 'redness', 'dry', 'oozing', 'edema') ], method = 'spearman')
corrplot(M, method = "number", type = 'upper', order = 'AOE', col =  col2(100), cl.lim = c(0, 1))