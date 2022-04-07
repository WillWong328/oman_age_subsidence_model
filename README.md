# oman_age_subsidence_model
An age-depth and subsidence model based on the Nafun Group stratigraphy of the MIQRAT-1 well of Oman


This model has four components:
1) A bayesian age-depth updater that enforces the law of superposition.
2) A bayesian age-model generator using Modified BChron (Trayler et al., 2019)
3) A decompaction and backstripping component for making porosity-loss and sediment loading corrections.
4) A bayesian thermal subsidence model generator.


All model plots are generated from the main.R script. Age-depth data is inputted via the data_set.csv
file, and the depths of important lithological boundaries are inputted through the important_depths.csv
file. All files and scripts need to placed in the same folder for this program to run.

The sample data set here is for the Nafun Group stratigraphy based on the MIQRAT-1 well of Oman. Other 
sections and/or drill cores can also be used with this model, but data and specific settings in the 
main.R script will need to be adjusted accordingly. More details of this model can be found in the
Masters Thesis of Wong (2022).
