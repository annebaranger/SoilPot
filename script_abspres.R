# _targets.R file

#library
library(targets)

#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","sf","lubridate","lme4"),
               error = "continue") #"KrigR","gdalUtils",

# df.traits <- tar_read(df.traits,store="target_safetymargin")


#Targets
list(
  tar_target(
    traits.file,
    "target_traits/objects/traits_max",
    format="file"
  ),
  tar_target(
    traits,
    readRDS(traits.file)
  ),
  # regular psi_min with beta=0.97
  tar_target(
    psi_cerraday_real_file,
    "target_psi2/objects/psi_cerraday_real",
    format="file"
  ),
  tar_target(
    psi_cerraday_real,
    readRDS(psi_cerraday_real_file)
  ),
  # tar_target(
  #   psi_cerraday_beta_file,
  #   "target_psi2/objects/",
  #   format="file"
  # ),
  # tar_target(
  #   psi_cerraday_real,
  #   readRDS(psi_cerraday_real_file)
  # ),
  tar_target(
    psi_eraday_real_file,
    "target_psi2/objects/psi_eraday_real",
    format="file"
  ),
  tar_target(
    psi_eraday_real,
    readRDS(psi_eraday_real_file)
  ),
  tar_target(
    tmin_cerra_file,
    "target_data/objects/tmin_cerra",
    format="file"
  ),
  tar_target(
    tmin_cerra,
    readRDS(tmin_cerra_file)
  ),
  tar_target(
    tmin_chelsa_file,
    "target_data/objects/tmin_chelsa",
    format="file"
  ),
  tar_target(
    tmin_chelsa,
    readRDS(tmin_chelsa_file)
  ),
  tar_target(
    clim_chelsa_file,
    "target_data/objects/clim_chelsa",
    format="file"
  ),
  tar_target(
    clim_chelsa,
    readRDS(clim_chelsa_file)
  ),

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 9 - Load Mauri Data ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    db.mauri,
    get.mauri(dir.occ="data/EUForestsMauri/EUForestspecies.csv")
  ),
  tar_target(
    occurence,
    get.occ.clim(db.mauri,
                 clim.list=list(psi_cerraday_real=psi_cerraday_real[,c("x","y","psi")],
                                psi_eraday_real=psi_eraday_real[,c("x","y","psi")],
                                tmin_cerra=tmin_cerra,
                                tmin_chelsa=tmin_chelsa,
                                mat=clim_chelsa[,c("x","y","mat")],
                                map=clim_chelsa[,c("x","y","map")],
                                pet=clim_chelsa[,c("x","y","pet")],
                                sgdd=clim_chelsa[,c("x","y","sgdd")],
                                wai=clim_chelsa[,c("x","y","wai")]),
                 file.path="output/df.occurence.clim.csv")
  ),
  # tar_target(
  #   df.mauri,
  #   get.occurence(df.mauri.sfm,
  #                 df.traits,
  #                 file.path="output/db_EuForest.csv")
  # ),
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # #### Section 10 - Idem, with other psi_min computation ####
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # tar_target(
  #   df.mauri.sfm.varb,
  #   get.mauri(dir.occ="data/EUForestsMauri/EUForestspecies.csv",
  #             psi_min=psi_varb,
  #             psi_min_100=psi_min_100,
  #             frost.index=frost.index,
  #             file.path="output/db_EuForest_raw_varb.csv")
  # ),
  # 
  # tar_target(
  #   df.mauri.varb,
  #   get.occurence(df.mauri.sfm.varb,
  #                 df.traits,
  #                 file.path="output/db_EuForest_varb.csv")
  # ),
  # 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #### Section 11 - Compute niche index ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    species.list,
    list.files("data/chorological_maps_dataset/")
  ),
  tar_target(
    df.preval,
    get.prevalence(species.list,
                   occurence)
  ),
  tar_target(
    df.niche,
    get.niche(species.list,
              occurence,
              psi="psi_cerraday_real",
              tmin="tmin_cerra")
  ),
  tar_target(
    df.shadetol,
    get.shadetol(species.list,
                 occurence)
  ),
  tar_target(
    df.species,
    get.species(species.list,
                df.preval,
                df.shadetol,
                df.niche,
                file.output="output/df.species2.csv")
  ),
  NULL
)
