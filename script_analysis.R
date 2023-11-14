# _targets.R file

#library
library(targets)
library(future)
# library(tarchetypes)

#Options
source("R/functions_analyses.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","data.table","tidyr","viridis",
                            "rgdal","raster","rosm","terra","dplyr","sf",
                            "lubridate","lme4","rstan","pROC"),
               error = "continue") 
tar_option_set(memory = "transient")
future::plan(future::multisession, workers = 6)

#Targets
list(
  # Load necessary files
  tar_target(
    df.species.file,
    "target_occ/objects/df.species",
    format="file"
  ),
  tar_target(
    df.species,
    readRDS(df.species.file)
  ),
  
  # with beta fixed
  tar_target(
    occurence.file,
    "target_occ/objects/occurence", 
    format="file"
  ),  
  tar_target(
    occurence,
    readRDS(occurence.file)
  ),
  # with beta varying
  tar_target(
    occurence_beta.file,
    "target_occ/objects/occurence", 
    format="file"
  ),  
  tar_target(
    occurence_beta,
    readRDS(occurence_beta.file)
  ),
  # run model with safety margins
  
  # tar_target(
  #   fit.safmarg,
  #   fit.logistic(occurence,
  #                var.hsm="psi_eraday_real",
  #                var.fsm="tmin_era",
  #                df.species,
  #                file.path="output/output_safmarg_era.csv",
  #                output="fit_safmarg_era/"
  #                )
  # ),
  tar_target(
    sp.list,
    df.species |> 
      filter(species %in% unique(occurence$species)) |> 
      filter(!is.na(px)) |> 
      filter(!is.na(lt50)) |> 
      pull(species)
  ),
  tar_target(
    fit.safmarg.sp,
    fit.logistic.dynamic(occurence=occurence,
                         var.hsm="psi_eraday_real",
                         var.fsm="tmin_era",
                         df.species=df.species,
                         sp=sp.list,
                         output="fit_safmarg_era/"),
    pattern=map(sp.list)
  ),
  tar_target(
    saf.marg.output,
    fwrite(fit.safmarg.sp,file="output/output_safmarg_era.csv")
  ),
  # tar_target(
  #   df.output,
  #   get.output(db.clim.file,
  #              df.species,
  #              output.path="fit_mod6/",
  #              "output/df.ouput2.csv")
  # ),
  # tar_target(
  #   df.output.auc,
  #   compute.auc(df.output,
  #               db.clim.file)
  # ),
  tar_target(
    mod.select.safmarg,
    select.model(fit.safmarg.sp)
  ),
  
  
  # run model with safety margins
  tar_target(
    fit.safmarg.beta.sp,
    fit.logistic.dynamic(occurence=occurence_beta,
                         var.hsm="psi_eraday_real",
                         var.fsm="tmin_era",
                         df.species=df.species,
                         sp=sp.list,
                         output="fit_safmarg_era_beta/"),
    pattern=map(sp.list)
  ),
  tar_target(
    saf.marg_beta.output,
    fwrite(fit.safmarg.beta.sp,file="output/output_safmarg_era.csv")
  ),
  tar_target(
    mod.select.safmarg_beta,
    select.model(fit.safmarg.beta.sp)
  ),
  
  # with varying beta
  # tar_target(
  #   db.clim.file.varb,
  #   "output/db_EuForest_varb.csv", # a changer --> mette chemin du tar_object
  #   format="file"
  # ),
  # tar_target(
  #   fit.mod.varb,
  #   fit.logistic(db.clim.file.varb,
  #                df.species,
  #                soil.depth="real",
  #                output="fit_mod_varb/")
  # ),
  # tar_target(
  #   df.output.varb,
  #   get.output(db.clim.file.varb,
  #              df.species,
  #              output.path="fit_mod_varb/",
  #              "output/df.ouput2.csv")
  # ),
  # tar_target(
  #   df.output.auc,
  #   compute.auc(df.output,
  #               db.clim.file)
  # ),
  # tar_target(
  #   df.mod.select.varb,
  #   select.model(df.output.varb)
  # ),

  # tar_render(
  #   report,
  #   "Report_jan.Rmd"
  # ),
  
  # fit with climatic variable
  tar_target(
    fit.clim,
    fit.logistic.clim(occurence,
                 var.hsm="psi_eraday_real",
                 var.fsm="tmin_era",
                 df.species,
                 file.path="output/output_clim_chelsa.csv",
                 output="fit_clim_chelsa/"
    )
  ),
  tar_target(
    mod.select.clim,
    select.model(fit.clim)
  ),
  NULL

)