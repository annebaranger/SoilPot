# _targets.R file

#library
library(targets)
# library(tarchetypes)
# lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis",
                            "rgdal","raster","rosm","terra","dplyr","sf",
                            "lubridate","lme4","rstan","caret","pROC"),
               error = "continue") 

#Targets
list(
  # Load necessary files
  tar_target(
    df.species.file,
    "output/df.species.csv",
    format="file"
  ),
  tar_target(
    df.species,
    read.csv(df.species.file)
  ),
  
  # with beta fixed
  tar_target(
    db.clim.file,
    "output/db_EuForest.csv", # a changer --> mette chemin du tar_object
    format="file"
  ),
  # run models
  # tar_target(
  #   fit.mod,
  #   fit.logistic(db.clim.file,
  #                df.species,
  #                soil.depth="real",
  #                output="fit_mod6/")
  # ),
  tar_target(
    df.output,
    get.output(db.clim.file,
               df.species,
               output.path="fit_mod6/",
               "output/df.ouput2.csv")
  ),
  tar_target(
    df.output.auc,
    compute.auc(df.output,
                db.clim.file)
  ),
  tar_target(
    df.mod.select,
    select.model(df.output.auc)
  ),
  
  
  # with varying beta
  tar_target(
    db.clim.file.varb,
    "output/db_EuForest_varb.csv", # a changer --> mette chemin du tar_object
    format="file"
  ),
  tar_target(
    fit.mod.varb,
    fit.logistic(db.clim.file.varb,
                 df.species,
                 soil.depth="real",
                 output="fit_mod_varb/")
  ),
  tar_target(
    df.output.varb,
    get.output(db.clim.file.varb,
               df.species,
               output.path="fit_mod_varb/",
               "output/df.ouput2.csv")
  ),
  # tar_target(
  #   df.output.auc,
  #   compute.auc(df.output,
  #               db.clim.file)
  # ),
  tar_target(
    df.mod.select.varb,
    select.model(df.output.varb)
  ),

  # tar_render(
  #   report,
  #   "Report_jan.Rmd"
  # ),
  
  # fit with climatic variable
  tar_target(
    fit.mod.clim,
    fit.logistic.clim(db.clim.file,
                      df.species,
                      output="fit_modClim/",
                      file.path="output/df.outputClim.csv")
  ),
  NULL

)