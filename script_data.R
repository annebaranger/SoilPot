# _targets.R file

#library
library(targets)
# lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") 

#Targets
list(
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 1 - Downloading climatic and pedologic data ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  #' @description Create europe mask according to countries of interest
  #'
  tar_target(
    europe,
    europe_extent()
  ),

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 2 - Downloading soil water content to targetted depth ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' @description SWC min from ERA5, monthly timestep, all horzions, perc05
  #'
  tar_target(
    swc_era_month,
    load_swc(dir.data="data/ERA5-land/monthly/",
             dir.file="era5_svwcperc05_",
             vars=c("h1","h2","h3","h4"),
             extension="_1950-2021_inv",
             format="grib",
             europe)
  ),
  
  #' @description SWC min from ERA5, daily timestep, all horzions, mean of 5 min
  #'
  tar_target(
    swc_era_day,
    load_swc(dir.data="data/ERA5-land/daily/",
             dir.file="swcd-1950-2021-",
             vars=c("layer1","layer2","layer3","layer4"),
             extension="",
             format="csv",
             europe)
  ),
  #' @description SWC min from cerra, daily timestep, 11 horzions, perc05
  #'
  tar_target(
    swc_cerra_day,
    load_swc(dir.data="data/cerra-land/liquid_vol_content/min/",
             dir.file="cerra_lvslperc05_",
             vars=paste0("h",1:11),
             extension="_1984-2021_inv",
             format="grib",
             europe)
  ),
  tar_target(
    tmin_cerra,
    get_frostindex_cerra(europe,
                         dir.file="data/cerra-land/skin_temp/cerra_perc05_1984-2021_inv.grib",
                         output.tmin="output/tmin_cerra.csv")
  ),
  tar_target(
    tmin_era,
    get_frostindex_chelsa(europe,
                          dir.file="data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc",
                          output.tmin="output/tmin_chelsa.csv")
  ),
  NULL
)
