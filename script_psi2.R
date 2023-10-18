# _targets.R file

#library
library(targets)
library(tarchetypes)
# lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") 

# Load objects from other projects
europe <- tar_read(europe, store = "target_data")

values <- data.frame(LAImax=c(2,3,4,5,5,5,5,6,7,8),
                     betaRootProfile=c(0.966,0.966,0.966,0.914,0.942,0.966,0.976,0.966,0.966,0.966))
          # expand.grid(LAImax=3:6,
          #             betaRootProfile=c(0.914,0.943,0.966,0.976))

#Targets
list(
  ## LOAD NECESSARY TARGETS ##
  tar_target(
    swc_era_month_file,
    "target_data/objects/swc_era_month",
    format="file"
  ),
  tar_target(
    swc_era_month,
    readRDS(swc_era_month_file)
  ),
  tar_target(
    swc_month_289,
    weight_swc(swc_era_month,
               depth=289)
  ),
  tar_target(
    swc_month_100,
    weight_swc(swc_era_month,
               depth=100)
  ),
  tar_target(
    swc_era_day_file,
    "target_data/objects/swc_era_day",
    format="file"
  ),
  tar_target(
    swc_era_day,
    readRDS(swc_era_day_file)
  ),
  tar_target(
    swc_cerra_day_file,
    "target_data/objects/swc_cerra_day",
    format="file"
  ),
  tar_target(
    swc_cerra_day,
    readRDS(swc_cerra_day_file)
  ),
  ## COMPUTE PSI_MIN ##
  tar_target(
    psi_eramonth_real,
    compute_psi_sureau(swc_era_month,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.07,0.28,1,2.89),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_era_month_real_fixed.csv"
    )
  ),
  tar_target(
    psi_eraday_real,
    compute_psi_sureau(swc_era_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.07,0.28,1,2.89),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_era_day_real_fixed.csv"
    )
  ),
  tar_target(
    psi_eraday_100,
    compute_psi_sureau(swc_era_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=100,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.07,0.28,1,2.89),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_era_day_100_fixed.csv"
    )
  ),
  tar_target(
    psi_cerraday_real,
    compute_psi_sureau(swc_cerra_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_cerra_day_real_fixed.csv"
    )
  ),
  tar_target(
    psi_cerraday_100,
    compute_psi_sureau(swc_cerra_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=100,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_cerra_day_100_fixed.csv"
    )
  ),
    #  # With SWC min on daily timestep (Python code), using real depth and SUREAU
    # tar_target(
    #   psihorday_real,
    #   compute_psihorday(SWC_day,
    #                     europe,
    #                     dir.data="data",
    #                     dir.file="EU_SoilHydroGrids_1km",
    #                     depth="real",
    #                     file.output="output/psihorday_real.csv")
    # ),
    # With SWC min on daily timestep (Python code), until 100cm and SUREAU
    # tar_target(
    #   psihorday_100,
    #   compute_psihorday(SWC_day,
    #                     europe,
    #                     dir.data="data",
    #                     dir.file="EU_SoilHydroGrids_1km",
    #                     depth=100,
    #                     file.output="output/psihorday_100.csv")
    # ),
  # sensitivity analysis
  # tar_map(
  #   values=values,
  #   tar_target(sensitivity,
  #              compute_sensitivity(SWC_day,
  #                       europe,
  #                       dir.data="data",
  #                       dir.file="EU_SoilHydroGrids_1km",
  #                       depth="real",
  #                       LAImax,
  #                       betaRootProfile))
  #   ),
  # tar_target(
  #   psihorday_rbeta,
  #   compute_psihorday_beta(SWC_day,
  #                     europe,
  #                     dir.data="data",
  #                     dir.file="EU_SoilHydroGrids_1km",
  #                     LAImax=5,
  #                     file.output="output/psihorday_varb.csv")
  # ),
  NULL
)