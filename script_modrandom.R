# _targets.R file

#libraries
library(targets)
library(future)

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
  # species
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
  # fit model with random effect
  tar_target(
    output.clim.file,
    "output/output_clim_chelsa.csv",
    format="file"
  ),
  tar_target(
    fit.random.safmarg,
    fit.allspecies.safmarg(occurence,
                   var.hsm="psi_eraday_real",
                   var.fsm="tmin_era",
                   df.species,
                   output.clim.file,
                   output="output/",
                   mod.folder="fit_allrandom/")
  ),
  tar_target(
    fit.random.clim,
    fit.allspecies.clim(occurence,
                           var.hsm="psi_eraday_real",
                           var.fsm="tmin_era",
                           df.species,
                           output.clim.file,
                           output="output/",
                           mod.folder="fit_allrandom/")
  ),
  tar_target(
    auc.modrandom,
    auc.modrandom(occurence,
                var.hsm="psi_eraday_real",
                var.fsm="tmin_era",
                df.species,
                fit.random.safmarg,
                fit.random.clim,
                output="output/")
  ),
  NULL
  
)