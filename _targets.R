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
  # Create europe mask according to countries of interest
  tar_target(
    europe,
    europe_extent()
  ),
  # tar_target(
  #   era_download,
  #   era_data(Dir.Data="data",
  #            Dir.Shapes="data/Shapes",
  #            extent = europe,
  #            variables=c("Volumetric_soil_water_layer_1","Volumetric_soil_water_layer_2","Volumetric_soil_water_layer_3","Volumetric_soil_water_layer_4"),
  #            date_start="1950-01-01",
  #            date_end="2021-12-31",
  #            time_step="month")
  # ),
  tar_target(
    terraclimate,
    terraclimate_data()
  ),
  tar_target(
    texture,
    texture_data(Dir.Data="data",
                 Dir.Soil="STU_EU_Layers",
                 europe)
  ),
  tar_target(
    forestcover,
    forest_cover(Dir.Data="data",
                 Dir.file="TCD_2018_100m_eu_03035_v020/DATA/TCD_2018_100m_eu_03035_V2_0.tif",
                 texture)
  ),
  tar_target(
    SWC,
    volumetric_content(Dir.Data="data")
  ),
  tar_target(
    psi_min,
    soil_potential(texture,SWC[[1]])
  ),
  tar_target(
    psiforest,
    psi_forest(psi_min,
               forestcover,
               40)
  ),
  tar_target(
    chronology_weighted,
    chronology_swc(SWC=SWC[[2]])
  ),
  tar_target(
    HSMs,
    HSM_distribution(Dir.Data="data",
                      Dir.p50="p50select.csv",
                      europe,
                      psi_min)
  )
)
