# _targets.R file

#library
library(targets)

#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") #"KrigR","gdalUtils",

#Targets
list(
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
  NULL
)