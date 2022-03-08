# _targets.R file

#library
library(targets)

#Options
source("R/functions_data.R")
source("R/functions_plot.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("KrigR", "ggplot2","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils"))

#Targets
list(
  tar_target(
    era_download,
    era_data(Dir.Data="data",
             Dir.Shapes="data/Shapes",
             countries = "France",
             variables=c("Volumetric_soil_water_layer_1","Volumetric_soil_water_layer_2","Volumetric_soil_water_layer_3","Volumetric_soil_water_layer_4"),
             date_start="2000-01-01",
             date_end="2015-01-01",
             time_step="month")
  ),
  tar_target(
    texture,
    texture_data(Dir.Data="data",
                 Dir.Soil="STU_EU_Layers")
  ),
  tar_target(
    SWC,
    volumetric_content(Dir.Data="data",
                       Variables=c("Volumetric_soil_water_layer_1","Volumetric_soil_water_layer_2","Volumetric_soil_water_layer_3","Volumetric_soil_water_layer_4"),
                       Abv=c("SWC1","SWC2","SWC3","SWC4"),
                       texture)
  ),
  tar_target(
    psi_min,
    soil_potential(texture,SWC[[1]])
  ),
  tar_target(
    chronology_weighted,
    chronology_swc(SWC=SWC[[2]])
  ),
  tar_target(
    chronology_mean,
    chronology_swc(SWC=SWC[[3]])
  )
)
