
# author : Valentin Delorieux


# import_data() : import and clean the dataset
# import_data_model(model) : import the dataset and transform it to be fitted by the specified model

library(tidyverse)

import_data <- function() {
  
  # This script import the dataset from the CSV file 2017Plos_breakdown.csv
  # It filters the patient AD - 057
  # It retrieves all the missing features using the mean of the features over the other patients for the same day
  # It removes duplicates
  # It creates a regularised time scale "day"
  # It a flag "first_obs"
  
  
  # Ouput :
  # Tibble of the cleaned data with 2 new calculated colums 
  # day : defined as the number of days relatively to the first observation for each patient
  # first_obs : flag for dates following periods of non recording 
  
  
  df <- read_csv2(
    "2017Plos_breakdown.csv",
    col_types = cols(
      date = col_date(format = "%m/%d/%Y"),
      season = col_character(),
      DOW = col_character(),
      fever = col_character(),
      sex = col_character(),
      symptom = col_integer()
    )
  )
  
  df <- filter(df, ID != "AD-057")   # Contains only 1 data point
  
  # Temperature --------------
  vars <- c("temperature", "DTR", "rainfall", "O3", "NO2", "PM10")
  N_obs <- length(df$ID)
  for (var in vars){
    for (n_obs in seq(1, N_obs)) {
      if (is.na(df[[n_obs, var]])) {
        df[n_obs, var] <- df %>% filter(date == date[[n_obs]]) %>% select(var) %>% as_vector() %>% mean(na.rm = TRUE)
      }
    }
  }
  
  # Fever and TCS ; replace missing values with 0
  # We checked for each missing fever if putting 0 was relevant with surrounding
  df <- df %>% replace_na(list(fever = 0, TCS = 0))
  
  # Replace aberrant values for fever
  df[which(df$fever > 1), "fever"] = "1"
  
  # Deletion of the duplicate values
  df <-  unique(df)
  
  # Creation of a patient-specific day scale taking into account the missing dates but adjusting for start of the trial
  df <- df %>%
    group_by(ID) %>%
    mutate(day = as.integer(difftime(date, min(date), units = "days") + 1))
  
  # Create new flag column : first_obs(t,k) = 1 if we don't have a data point at t-1 for patient k 
  df <- df %>% 
    arrange(ID, date) %>%
    mutate(first_obs = 0) 
  
  df[[1, "first_obs"]] = 1
  for (i in seq(2, length(df$date))) {
    # Put 1 for dates after a gap
    if (df$date[[i-1]] + 1 != df$date[[i]]) {df$first_obs[[i]] = 1}
    # Put 1 for first observation of each patient
    if (df$day[[i]] == 1) {df$first_obs[[1]] = 1}
  }
  
  # Remove observations with missing symptoms
  df <- drop_na(df)
  
  return(df)
}

import_data_model <- function(model = c("MOLR", "LR")) {
  
  model <- match.arg(model)
  
  # If model is logistic regression (model of Ahn et. Al(2017))
  if (model == "LR") {
    return(import_data())
  }
  
  # If model is mixed effect ordinal logistic regression
  if (model == "MOLR") {
    df <- import_data()
    
    # Gather all pairs of consecutive observations in a tibble, with symptom state S, sign scores s_i, and environmental factors
    df_model <- do.call(rbind, lapply(1:(length(df$ID)-1),
                                      function(n_obs) {
                                        if (df$ID[n_obs] == df$ID[n_obs + 1] & df$first_obs[n_obs + 1] == 0) {
                                          return(tibble(ID = df$ID[n_obs],
                                                        day = as.integer(df$day[n_obs]),
                                                        TCS = df$TCS[n_obs],
                                                        temperature = round(df$temperature[n_obs], 1),
                                                        rainfall = round(df$rainfall[n_obs], 1),
                                                        DTR = round(df$DTR[n_obs], 1),
                                                        RH = round(df$RH[n_obs], 1),
                                                        PM10 = round(df$PM10[n_obs], 1),
                                                        NO2 = round(df$NO2[n_obs], 1),
                                                        O3 = round(df$O3[n_obs], 1),
                                                        SCORAD = as.integer(df$SCORAD[n_obs]),
                                                        fever = as.integer(df$fever[n_obs]),
                                                        itching_t0 = df$itching[n_obs],
                                                        itching_t1 = df$itching[n_obs + 1],
                                                        sleep_t0 = df$sleep[n_obs],
                                                        sleep_t1 = df$sleep[n_obs + 1],
                                                        redness_t0 = df$redness[n_obs],
                                                        redness_t1 =  df$redness[n_obs + 1],
                                                        dry_t0 =  df$dry[n_obs],
                                                        dry_t1 = df$dry[n_obs + 1],
                                                        oozing_t0 = df$oozing[n_obs],
                                                        oozing_t1 = df$oozing[n_obs + 1],
                                                        edema_t0 = df$edema[n_obs],
                                                        edema_t1 = df$edema[n_obs + 1],
                                                        S_t0 = as.integer(df$symptom[n_obs]),
                                                        S0_t0 = as.integer(df$symptom[n_obs] == "0"),
                                                        S1_t0 = as.integer(df$symptom[n_obs] == "1"),
                                                        S_t1 = as.integer(df$symptom[n_obs + 1])))
                                        }
                                      }))
    
    # 15 patients do no have two consecutive data points on first two days of the beginning of the trial
    # We artificially add one observation for each of these patients to avoid mixed effect artefacts on MOLR fitting 
    IDs <- unique(df_model$ID)
    ID_okay <- unique(filter(df_model, day == 1)$ID)
    for (i in seq_along(IDs)){
      if (!(IDs[i] %in% ID_okay)) {
        df_model <- add_row(df_model,
                            ID = IDs[i], 
                            TCS = 0,
                            temperature = round(mean(df$temperature), 1),
                            rainfall = round(mean(df$rainfall), 1),
                            DTR = round(mean(df$DTR), 1),
                            RH = round(mean(df$RH), 1),
                            fever = 0,
                            SCORAD = round(mean(df$SCORAD, 1)),
                            PM10 = round(mean(df$PM10), 1),
                            NO2 = round(mean(df$NO2), 1),
                            O3 = round(mean(df$O3), 1),
                            day = 1, 
                            S_t0 = 0, S0_t0 = 1, S1_t0 = 0, S_t1 = 0,
                            itching_t0 = 0, itching_t1 = 0,
                            sleep_t0 = 0, sleep_t1 = 0,
                            redness_t0 = 0, redness_t1 = 0,
                            dry_t0 = 0, dry_t1 = 0,
                            oozing_t0 = 0, oozing_t1 = 0,
                            edema_t0 = 0, edema_t1 = 0
        )
      }
    }
    
    # An ordinal logistic regression can be performed by jointly fitting the cumulative distributions P(y<=j) with (0<= j <=3)
    # We stack here the dataset 4 times to fit these 4 distributions
    # The outcomes of each logistic regression will be the cumulative distributions sign_t1_cumu_y
    len <- length(df_model$ID)
    stacked <- rbind(df_model, df_model, df_model, df_model)
    stacked$j <- c(rep(0, len), rep(1,len), rep(2, len), rep(3, len))
    stacked <- mutate(stacked, 
                      itching_t1_cumu_y = as.integer(itching_t1 <= j),
                      sleep_t1_cumu_y = as.integer(sleep_t1 <= j),
                      redness_t1_cumu_y = as.integer(redness_t1 <= j),
                      dry_t1_cumu_y = as.integer(dry_t1 <= j),
                      oozing_t1_cumu_y = as.integer(oozing_t1 <= j),
                      edema_t1_cumu_y = as.integer(edema_t1 <= j))
    
    return(arrange(stacked, ID, day, j))
  }
  
}
