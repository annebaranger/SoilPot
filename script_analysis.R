# _targets.R file

#library
library(targets)
library(tarchetypes)
# lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis",
                            "rgdal","raster","rosm","terra","dplyr","sf",
                            "lubridate","lme4","rstan"),
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
  tar_target(
    db.clim.file,
    "output/db_EuForest.csv",
    format="file"
  ),
  
  # run models
  # tar_target(
  #   fit.mod,
  #   fit.logistic(db.clim.file,
  #                df.species)
  # ),
  tar_target(
    df.output,
    get.output(db.clim.file,
               df.species,
               output.path="fit_mod6/")
  ),
  tar_target(
    df.mod.select,
    select.model(df.output)
  ),
  # tar_render(
  #   report,
  #   "Report_jan.Rmd"
  # ),
  NULL

)