# _targets.R file

#library
library(targets)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") #"KrigR","gdalUtils",

# Load objects from other projects
europe <- tar_read(europe, store = "target_data")
SWCtot <- tar_read(SWCtot,store = "target_data")
textureERA5 <- tar_read(textureERA5,store = "target_data")
textureESDAC <- tar_read(textureESDAC,store = "target_data")


#Targets
list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 2 - Downloading soil water content to targetted depth ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' @description SWC weighted until 289cm, extrema among months
  #'
  tar_target(
    SWC289,
    extr_swc(weight_swc(SWCtot,
                        depth=289),
             "month",
             "1949-12-01",
             europe)
  ),
  
  #' @description SWC weighted until 100cm, timeseries, monthly
  #'
  tar_target(
    SWC100t,
    weight_swc(SWCtot,
               depth=100)
  ),
  #' @description SWC weighted until 100cm, extrema among months
  #'
  tar_target(
    SWC100,
    extr_swc(SWC100t,
             "month",
             "1949-12-01",
             europe)
  ),
  
  #' @description SWC weighted until 50cm, extrema among months
  #'
  tar_target(
    SWC50,
    extr_swc(weight_swc(SWCtot,
                        depth=50),
             "month",
             "1949-12-01",
             europe)
  ),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 3 - Compute psi min with weighted swc####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  #' @description ERA5 x d50
  #'
  tar_target(
    psi_ERA_50,
    compute_psiweighted(textureERA5,SWC50,
                        file.output="output/psi_ERA_50.csv")
  ),
  #' @description ERA5 x d100
  #'
  tar_target(
    psi_ERA_100,
    compute_psiweighted(textureERA5,SWC100,
                        file.output="output/psi_ERA_100.csv")
  ),
  #' @description ERA5 x d289
  #'
  tar_target(
    psi_ERA_289,
    compute_psiweighted(textureERA5,SWC289,
                        file.output="output/psi_ERA_289.csv")
  ),
  #' @description ESDAC x topsoil x d50
  
  tar_target(
    psi_ESDACt_50,
    compute_psiweighted(textureESDAC$texture$topsoil,SWC50,
                        file.output="output/psi_ESDACt_50.csv")
  ),
  #' @description ESDAC x topsoil x d100
  #'
  tar_target(
    psi_ESDACt_100,
    compute_psiweighted(textureESDAC$texture$topsoil,SWC100,
                        file.output="output/psi_ESDACt_100.csv")
  ),
  #' @description ESDAC x subsoil x d50
  #'
  tar_target(
    psi_ESDACs_50,
    compute_psiweighted(textureESDAC$texture$subsoil,SWC50,
                        file.output="output/psi_ESDACs_50.csv")
  ),
  #' @description ESDAC x subsoil x d100
  #'
  tar_target(
    psi_ESDACs_100,
    compute_psiweighted(textureESDAC$texture$subsoil,SWC100,
                        file.output="output/psi_ESDACs_100.csv")
  ),
  #' @description ESDAC x subsoil x d289
  #'
  tar_target(
    psi_ESDACs_289,
    compute_psiweighted(textureESDAC$texture$subsoil,SWC289,
                        file.output="output/psi_ESDACs_289.csv")
  ),

  NULL
)