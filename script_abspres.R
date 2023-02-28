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
    "output/df_trait_filtered.csv",
    format="file"
  ),
  tar_target(
    df.traits,
    read.csv(df.traits.file)
  ),
  tar_target(
    psi_min,
    "output/psihorday_real.csv",
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
              df.traits,
              psi_min,
              frost.index)
  ),

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 10 - Compute niche index ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    df.preval,
    get.prevalence(df.traits,
                   df.mauri.sfm)
  ),
  tar_target(
    df.niche,
    get.niche(df.traits,
              df.mauri.sfm)
  ),
  tar_target(
    df.shadetol,
    get.shadetol(df.traits)
  ),
  tar_target(
    df.species,
    get.species(df.traits,
                df.preval,
                df.shadetol,
                df.niche)
  ),
  NULL
)