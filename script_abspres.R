# _targets.R file

#library
library(targets)

#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") #"KrigR","gdalUtils",

# df.traits <- tar_read(df.traits,store="target_safetymargin")


#Targets
list(
  tar_target(
    df.traits.file,
    "target_safetymargin/objects/df.traits",
    format="file"
  ),
  tar_target(
    df.traits,
    readRDS(df.traits.file)
  ),
  # regular psi_min with beta=0.97
  tar_target(
    psi_min,
    "output/psihorday_real.csv",
    format="file"
  ),
  # regular psi_min with varying beta
  tar_target(
    psi_varb,
    "output/psihorday_varb.csv",
    format="file"
  ),
  tar_target(
    psi_min_100,
    "output/psihorday_100.csv",
    format="file"
  ),
  tar_target(
    frost.index,
    "output/budburst_tquant.csv",
    format="file"
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 9 - Load Mauri Data ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    df.mauri.sfm,
    get.mauri(dir.occ="data/EUForestsMauri/EUForestspecies.csv",
              psi_min=psi_min,
              psi_min_100=psi_min_100,
              frost.index=frost.index,
              file.path="output/db_EuForest_raw.csv")
  ),
  
  tar_target(
    df.mauri,
    get.occurence(df.mauri.sfm,
                  df.traits,
                  file.path="output/db_EuForest.csv")
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 10 - Idem, with other psi_min computation ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    df.mauri.sfm.varb,
    get.mauri(dir.occ="data/EUForestsMauri/EUForestspecies.csv",
              psi_min=psi_varb,
              psi_min_100=psi_min_100,
              frost.index=frost.index,
              file.path="output/db_EuForest_raw_varb.csv")
  ),
  
  tar_target(
    df.mauri.varb,
    get.occurence(df.mauri.sfm.varb,
                  df.traits,
                  file.path="output/db_EuForest_varb.csv")
  ),

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 11 - Compute niche index ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    species.list,
    list.files("data/chorological_maps_dataset/")
  ),
  tar_target(
    df.preval,
    get.prevalence(species.list,
                   df.mauri)
  ),
  tar_target(
    df.niche,
    get.niche(species.list,
              df.mauri)
  ),
  tar_target(
    df.shadetol,
    get.shadetol(species.list,
                 df.mauri)
  ),
  tar_target(
    df.species,
    get.species(species.list,
                df.preval,
                df.shadetol,
                df.niche)
  ),
  NULL
)