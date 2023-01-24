# _targets.R file

#library
library(targets)
# lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") #"KrigR","gdalUtils",

# Load objects from other projects
europe <- tar_read(europe, store = "target_data")

#Targets
list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 5 - Frost index####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Frost index (CHELSA, FDG)
    tar_target(
      frost.index.quant,
      get.frostindex(europe,
                     dir.clim="data/CHELSA/CHELSA_EUR11_tasmin_month_quant05_19802005.nc",
                     dir.fdg="data/eea_2000-2016/SOS_2000_2016_DOY.BIL",
                     output.index="output/budburst_tquant.csv",
                     output.budburst="output/budburst_1.csv",
                     year_start=2000)
    ),
    # Frost index (CHELSA, FDG)
    tar_target(
      frost.index.min,
      get.frostindex(europe,
                       dir.clim="data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc",
                       dir.fdg="data/eea_2000-2016/SOS_2000_2016_DOY.BIL",
                       output.index="output/budburst_tmin.csv",
                       output.budburst="output/budburst_2.csv",
                      year_start=2000)
    ),
  NULL
)