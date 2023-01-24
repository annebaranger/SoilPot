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
psihorday_real <- tar_read(psihorday_real,store="target_psi2")

#Targets
list(
# Load necessary files
  tar_target(
    rast.temp,
    "output/budburst_tquant.csv",
    format="file"
  ),
  tar_target(
    rast.fdg.mean,
    "output/budburst_2.csv",
    format="file"
  ),
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 6 - Safety margins ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Frost safety margins
  tar_target(
    safety.margins,
    compute.sfm(df.traits,
                rast.temp,
                rast.fdg.mean,
                psihorday_real,
                europe)
  ),

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 7 - Load traits ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    df.LT50,
    get.LT50()
  ),
  tar_target(
    df.P50,
    get.P50()
  ),
  tar_target(
    df.traits,
    get.traits(df.P50=df.P50,
               df.LT50=df.LT50$df.LT50sp.cor)
  ),
 NULL
)