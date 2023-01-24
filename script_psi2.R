# _targets.R file

#library
library(targets)
# lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") #"KrigR","gdalUtils",

# Load objects from other projects
europe <- tar_read(europe, store = "target_data")
SWC_day <- tar_read(SWC_day,store="target_data")
SWCtot <- tar_read(SWCtot,store = "target_data")


#Targets
list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 4 - Compute psi_min swc and param per horizon####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # With SWC min on monthly timestep (R code), using all horizons and SUREAU
    tar_target(
      psihormonth_100,
      compute_psihor(SWCtot,
                     3,
                     "month",
                     "1949-12-01",
                     dir.data="data",
                     dir.file="EU_SoilHydroGrids_1km",
                     europe,
                     file.output="output/psihormonth_100.csv")
    ),
    # With SWC min on daily timestep (Python code), using real depth and SUREAU
    tar_target(
      psihorday_real,
      compute_psihorday(SWC_day,
                        europe,
                        dir.data="data",
                        dir.file="EU_SoilHydroGrids_1km",
                        depth="real",
                        file.output="output/psihorday_real.csv")
    ),
    # With SWC min on daily timestep (Python code), until 100cm and SUREAU
    tar_target(
      psihorday_100,
      compute_psihorday(SWC_day,
                        europe,
                        dir.data="data",
                        dir.file="EU_SoilHydroGrids_1km",
                        depth=100,
                        file.output="output/psihorday_100.csv")
    ),
  NULL
)