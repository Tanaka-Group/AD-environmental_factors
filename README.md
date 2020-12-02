# Impact of environmental factors in predicting daily severity scores of atopic dermatitis

This repository presents the code written for the paper by **Hurault et al. (2020)**, [Impact of environmental factors in predicting daily severity scores of atopic dermatitis](https://www.medrxiv.org/content/10.1101/2020.10.27.20220947v1) (preprint).
The code is written in R language for statistical computing.

## File structure 

The dataset that was used in the study originates from the paper ["Short-term effects of weather and air pollution on atopic dermatitis symptoms in children: A panel study in Korea" by Kim et al. (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0175229) and can be requested from the [Zenodo database](https://zenodo.org/deposit/124411/).

With the dataset file named `2017Plos_breakdown.csv`, the data can be loaded and process with the functions in [`001_data_import.R`](001_data_import.R).

Data exploration is performed in [`002_Data_exploration.R`](002_Data_exploration.R).

Model validation is conducted in [`003_LR_model.R`](003_LR_model.R) for the mixed effect logistic regression model that was previously published [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0175229) by Kim et al.(2017).; and in [`004_MOLR_model.R`](004_MOLR_model.R) for the mixed effect ordinal logistic regression model.

The project also contains visualisation scripts:

- [`005_plot_coef.R`](005_plot_coef.R): the script used to produce the plots presenting the coefficients associated to environmental factors and the difference of performance between models without environmental covariates and models with environmental covariates.
- [`006_plot_lear_curve.R`](006_plot_lear_curve.R): the script used to plot the learning curves.
