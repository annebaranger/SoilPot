# _targets.R file

#library
library(targets)
library(future)
# library(tarchetypes)

#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis",
                            "rgdal","raster","rosm","terra","dplyr","sf",
                            "lubridate","lme4","rstan","pROC"),
               error = "continue") 
tar_option_set(memory = "transient")
future::plan(future::multisession, workers = 6)

#Targets
list(
  # Load necessary files
  tar_target(
    df.species.file,
    "target_occ/objects/df.species",
    format="file"
  ),
  tar_target(
    df.species,
    readRDS(df.species.file)
  ),
  
  # with beta fixed
  tar_target(
    occurence.file,
    "target_occ/objects/occurence", 
    format="file"
  ),  
  tar_target(
    occurence,
    readRDS(occurence.file)
  ),
  # run model with safety margins
  tar_target(
    fit.safmarg,
    fit.logistic(occurence,
                 var.hsm="psi_eraday_real",
                 var.fsm="tmin_era",
                 df.species,
                 file.path="output/output_safmarg_era.csv",
                 output="fit_safmarg_era/"
                 )
  ),
  # tar_target(
  #   df.output,
  #   get.output(db.clim.file,
  #              df.species,
  #              output.path="fit_mod6/",
  #              "output/df.ouput2.csv")
  # ),
  # tar_target(
  #   df.output.auc,
  #   compute.auc(df.output,
  #               db.clim.file)
  # ),
  tar_target(
    mod.select.safmarg,
    select.model(fit.safmarg)
  ),
  
  
  # with varying beta
  # tar_target(
  #   db.clim.file.varb,
  #   "output/db_EuForest_varb.csv", # a changer --> mette chemin du tar_object
  #   format="file"
  # ),
  # tar_target(
  #   fit.mod.varb,
  #   fit.logistic(db.clim.file.varb,
  #                df.species,
  #                soil.depth="real",
  #                output="fit_mod_varb/")
  # ),
  # tar_target(
  #   df.output.varb,
  #   get.output(db.clim.file.varb,
  #              df.species,
  #              output.path="fit_mod_varb/",
  #              "output/df.ouput2.csv")
  # ),
  # tar_target(
  #   df.output.auc,
  #   compute.auc(df.output,
  #               db.clim.file)
  # ),
  # tar_target(
  #   df.mod.select.varb,
  #   select.model(df.output.varb)
  # ),

  # tar_render(
  #   report,
  #   "Report_jan.Rmd"
  # ),
  
  # fit with climatic variable
  tar_target(
    fit.clim,
    fit.logistic.clim(occurence,
                 var.hsm="psi_eraday_real",
                 var.fsm="tmin_era",
                 df.species,
                 file.path="output/output_clim_chelsa.csv",
                 output="fit_clim_chelsa/"
    )
  ),
  tar_target(
    mod.select.clim,
    select.model(fit.clim)
  ),
  
  # tar_target(
  #   fit.mod.clim,
  #   fit.logistic.clim(db.clim.file,
  #                     df.species,
  #                     output="fit_modClim/",
  #                     file.path="output/df.outputClim.csv")
  # ),
  NULL

)