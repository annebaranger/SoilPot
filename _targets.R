# _targets.R file

#library
library(targets)
#lapply(c("KrigR", "ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
source("R/functions_plot.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("KrigR", "ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf"))

#Targets
list(
## Create europe mask according to countries of interest ##
###########################################################
  tar_target(
    europe,
    europe_extent()
  ),
## Terraclimate, useless ##
###########################
  tar_target(
    terraclimate,
    terraclimate_data()
  ),
## Get texture data from 2 databases ##
#######################################
  tar_target(
    textureESDAC,
    get_textureESDAC(dir.data="data",
                 dir.soil="STU_EU_Layers",
                 europe)
  ),
  tar_target(
    textureERA5,
    get_textureERA5(dir.data="data",
                    dir.file="texture-2020.nc",
                    europe)
  ),
## Get forest cover at the right resolution ##
##############################################
  tar_target(
    forestcover,
    forest_cover(dir.data="data",
                 dir.file="TCD_2018_100m_eu_03035_v020/DATA/TCD_2018_100m_eu_03035_V2_0.tif",
                 textureESDAC$topsoil)
  ),
## Compute SWC according to depth required ##
#############################################
  #All horizons
  tar_target(
    SWC289,
    volumetric_content(dir.data="data",
                       dir.file="swc-1950-2021.nc",
                       europe,
                       depth=289)
  ),
  # 3 first horizons
  tar_target(
    SWC100,
    volumetric_content(dir.data="data",
                       dir.file="swc-1950-2021.nc",
                       europe,
                       depth=100)
  ),
  # 3 first horizons until 50
  tar_target(
    SWC50,
    volumetric_content(dir.data="data",
                       dir.file="swc-1950-2021.nc",
                       europe,
                       depth=50)
  ),
## Compute psi_min according to texture and swc ##
##################################################
  # ERA5 x d50
  tar_target(
    psi_ERA_50,
    soil_potential(textureERA5,SWC50[[1]])
  ),
  # ERA5 x d100
  tar_target(
    psi_ERA_100,
    soil_potential(textureERA5,SWC100[[1]])
  ),
  # ESDAC x topsoil x d50
  tar_target(
    psi_ESDACt_50,
    soil_potential(textureESDAC$topsoil,SWC50[[1]])
  ),
  # ESDAC x topsoil x d100
  tar_target(
    psi_ESDACt_100,
    soil_potential(textureESDAC$topsoil,SWC100[[1]])
  ),
  # ESDAC x subsoil x d50
  tar_target(
    psi_ESDACs_50,
    soil_potential(textureESDAC$subsoil,SWC50[[1]])
  ),
  # ESDAC x subsoil x d100
  tar_target(
    psi_ESDACs_100,
    soil_potential(textureESDAC$subsoil,SWC100[[1]])
  ),
  # ESDAC x subsoil x d289
  tar_target(
    psi_ESDACs_289,
    soil_potential(textureESDAC$subsoil,SWC289[[1]])
  ),
## Mask only forest plots ##
############################
  tar_target(
    psiforest,
    psi_forest(psi_ESDACs_289,
               forestcover,
               40)
  ),
  tar_target(
    chronology_weighted,
    chronology_swc(SWC=SWC289[[2]])
  ),
  tar_target(
    HSMs,
    HSM_distribution(dir.data="data",
                      dir.p50="p50select.csv",
                      europe,
                      psi_ESDACs_100)
  )
)
