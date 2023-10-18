# _targets.R file

#library
library(targets)

#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis",
                            "rgdal","raster","rosm","terra","dplyr","sf",
                            "lubridate","lme4","rstan","pROC"),
               error = "continue") 

#Targets
list(
  # with beta fixed
  tar_target(
    db.clim.file,
    "output/db_EuForest.csv", # a changer --> mette chemin du tar_object
    format="file"
  ),
  # fit model with random effect
  tar_target(
    output.clim.file,
    "output/df.outputClim.csv",
    format="file"
  ),
  tar_target(
    fit.random,
    fit.allspecies(db.clim.file,
                   output.clim.file,
                   output="output/",
                   mod.folder="fit_allrandom/")
  ),
  NULL
  
)