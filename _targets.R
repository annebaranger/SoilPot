# _targets.R file

#library
library(targets)
#lapply(c("ggplot2","stringr","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils","sf","lubridate","lme4"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
source("R/functions_plot.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4")) #"KrigR","gdalUtils",

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
#' @description Terraclimate, useless
#'
  tar_target(
    terraclimate,
    terraclimate_data()
  ),
#' @description Get Chelsa data
#'
  tar_target(
    chelsabio6,
    chelsabio6_data("data/CHELSA/","CHELSA_bio6_1981-2010_V.2.1.tif",europe)
  ),

#' @description Get IFN data
#'
  tar_target(
    db.fni,
    get.fni()
  ),

#' @description Get texture data from 2 databases
#'
#' ESDAC
  tar_target(
    textureESDAC,
    get_textureESDAC(dir.data="data",
                 dir.soil="STU_EU_Layers",
                 europe)
  ),

#' ERA5, pars from Wosten 1999, subsoil
  tar_target(
    textureERA5,
    get_textureERA5(dir.data="data",
                    dir.file="texture-2020.nc",
                    europe,
                    # pars from Wosten 1999, for SUBSOIL
                    pars=data.frame(texture=c(6,5,4,3,2,1),
                                    psi_e=c(-0.790569415,-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
                                    b=c(2.6411388301,2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933),
                                    teta_r=c(0.010,0.0250,0.0100,0.0100,0.0100,0.0100),
                                    teta_s=c(0.766,0.3660,0.3920,0.4120,0.4810,0.5380),
                                    alpha=c(0.0130,0.0430,0.0249,0.0082,0.0198,0.0168),
                                    n=c(1.2039,1.5206,1.1689,1.2179,1.0861,1.0730),
                                    m=c(0.1694,0.3424,0.1445,0.1789,0.0793,0.0680)))
  ),

#' ERA5, pars from Toth 2015, subsoil
  tar_target(
    textureERA5Toth,
    get_textureERA5(dir.data="data",
                    dir.file="texture-2020.nc",
                    europe,
                    # pars from Toth 2015, for SUBSOIL
                    pars=data.frame(texture=c(6,5,4,3,2,1),
                                    psi_e=c(-0.790569415,-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
                                    b=c(2.6411388301,2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933),
                                    teta_r=c(0.111,0.045,0.000,0.000,0.000,0.000),
                                    teta_s=c(0.697,0.438,0.459,0.432,0.478,0.522),
                                    alpha=c(0.0069,0.0478,0.0309,0.0094,0.0403,0.0112),
                                    n=c(1.4688,1.3447,1.1920,1.2119,1.1176,1.1433),
                                    m=c(0.3192,0.2563,0.1611,0.1749,0.1053,0.1253)))
  ),

#' @description Get forest cover at the right resolution
#'
  tar_target(
    forestcover,
    forest_cover(dir.data="data",
                 dir.file="TCD_2018_100m_eu_03035_v020/DATA/TCD_2018_100m_eu_03035_V2_0.tif",
                 textureESDAC$topsoil)
  ),

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Downloading soil water content to targetted depth ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @description SWC time series, all horizons, from 1950 to 2022, monthly
#'
tar_target(
    SWCtot,
    get_SWC(dir.data="data",
            dir.file="swc-1950-2021.nc"
            )
  ),

#' @description SWC weighted until 289cm, extrema among months
#'
  tar_target(
    SWC289,
    extr_swc(weight_swc(dir.data="data",
                        dir.file="swc-1950-2021.nc",
                        SWCtot,
                        depth=289),
             "month",
             "1949-12-01",
             europe)
  ),

#' @description SWC weighted until 100cm, timeseries, monthly
#'
  tar_target(
    SWC100t,
    weight_swc(dir.data="data",
               dir.file="swc-1950-2021.nc",
               SWCtot,
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
    extr_swc(weight_swc(dir.data="data",
                        dir.file="swc-1950-2021.nc",
                        SWCtot,
                        depth=50),
             "month",
             "1949-12-01",
             europe)
  ),

#' @description SWC min with daily timestep, all horizons, computed with python
#'
  tar_target(
    SWC_day,
    load_swc(dir.data="data",
             dir.file="ERA5-land/swcd-1950-2021-",
             vars=c("layer1","layer2","layer3","layer4"))
  ),

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Compute psi min with weighted swc####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @description ERA5 x d50
#'
  tar_target(
    psi_ERA_50,
    compute_psiweighted(textureERA5,SWC50)
  ),
#' @description ERA5 x d100
#'
  tar_target(
    psi_ERA_100,
    compute_psiweighted(textureERA5,SWC100)
  ),
#' @description ERA5 x d289
#'
  tar_target(
    psi_ERA_289,
    compute_psiweighted(textureERA5,SWC289)
  ),
#' @description ESDAC x topsoil x d50
#'
  tar_target(
    psi_ESDACt_50,
    compute_psiweighted(textureESDAC$topsoil,SWC50)
  ),
#' @description ESDAC x topsoil x d100
#'
  tar_target(
    psi_ESDACt_100,
    compute_psiweighted(textureESDAC$topsoil,SWC100)
  ),
#' @description ESDAC x subsoil x d50
#'
  tar_target(
    psi_ESDACs_50,
    compute_psiweighted(textureESDAC$subsoil,SWC50)
  ),
#' @description ESDAC x subsoil x d100
#'
  tar_target(
    psi_ESDACs_100,
    compute_psiweighted(textureESDAC$subsoil,SWC100)
  ),
#' @description ESDAC x subsoil x d289
#'
  tar_target(
    psi_ESDACs_289,
    compute_psiweighted(textureESDAC$subsoil,SWC289)
  ),


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
                   europe)
  ),
  # With SWC min on daily timestep (Python code), using real depth and SUREAU
  tar_target(
    psihorday_real,
    compute_psihorday(SWC_day,
                      europe,
                      dir.data="data",
                      dir.file="EU_SoilHydroGrids_1km",
                      depth="real")
  ),
  # With SWC min on daily timestep (Python code), until 100cm and SUREAU
  tar_target(
    psihorday_100,
    compute_psihorday(SWC_day,
                      europe,
                      dir.data="data",
                      dir.file="EU_SoilHydroGrids_1km",
                      depth=100)
  ),
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 5 - Compute safety margins####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # HSM (SWC month, param 3Dhydro, SUREAU) x FSM (CHELSA BOI6)
  tar_target(
    safety.margins,
    compute.sm(psihormonth_100,
               chelsabio6)
  ),
  # HSM (SWC day, param 3Dhydro, SUREAU, depth 100) x FSM (CHELSA BOI6)
  tar_target(
    safety.margins.day.100,
    compute.sm(psihorday_100,
               chelsabio6)
  ),
  # HSM (SWC day, param 3Dhydro, SUREAU, depth real) x FSM (CHELSA BOI6)
  tar_target(
    safety.margins.day.real,
    compute.sm(psihorday_real,
               chelsabio6)
  ),
  # Frost index (CHELSA, FDG)
  tar_target(
    frost.index,
    get.frostindex(europe,
                   year_start=2000)
  ),
  # Frost index (CHELSA, FDG)
  # tar_target(
  #   frost.index.bis,
  #   get.frostindex.bis(europe,
  #                  year_start=1995)
  # ),
  # Frost safety margins
  tar_target(
    frost.safety.margins,
    compute.fsm(dir.species="data/Species traits/trait.select.csv",
                frost.index$rast.temp,
                frost.index$rast.fdg.mean)
  ),

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 6 - Load traits ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    df.traits,
    get.traits()
  ),
  tar_target(
    list.LT50,
    get.LT50()
  ),

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 7 - Plot safety margins####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # tar_target(
  #   plots.sftymargins,
  #   plot.sftym(dir.data="data",
  #              dir.p50="Species traits/trait.select.csv",
  #              europe,
  #              safety.margins.day)
  # ),
#%%%%%%%%%%%%%%%%%%
#### Section Annexe
#%%%%%%%%%%%%%%%%%%

#' @description  Mask only forest plots
#'
  tar_target(
    psiforest,
    psi_forest(psi_ESDACs_289,
               forestcover,
               40)
  ),
#' @description timeseries
#'
  tar_target(
    chronology_weighted,
    chronology_swc(SWC=SWC100t)
  ),
  # tar_target(
  #   HSMs,
  #   HSM_distribution(dir.data="data",
  #                     dir.p50="p50select.csv",
  #                     europe,
  #                     psi_ESDACs_100)
  # )
  NULL
)