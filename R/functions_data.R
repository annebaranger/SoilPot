#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation, subsetting and computation of 
#               required variables
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Defining spatial extent of the study ####
#' @description Function used to create a polygon of study spatial extent
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Europe extent
#' 
#' @description This is a function to extract countries borders from an R 
#' package 
#' @param countries default selection encompass geographical definition of 
#' europe, except Russia
#' @return SpatialPolygonsDataFrame containing all countries of selection
#'
europe_extent <- function(countries=c("Albania","Austria","Belgium","Bulgaria",
                                   "Bosnia and Herzegovina","Switzerland",
                                   "Czech Republic","Germany","Denmark","Spain",
                                   "Estonia","Finland","France"
                                   ,"United Kingdom","Greece","Croatia",
                                   "Hungary","Ireland","Iceland","Italy",
                                   "Kosovo","Lithuania","Luxembourg","Latvia",
                                   "Macedonia","Montenegro","Netherlands",
                                   "Norway","Poland","Portugal","Romania",
                                   "Republic of Serbia","Slovakia","Slovenia",
                                   "Sweden","Andorra","Liechtenstein","Monaco",
                                   "Malta","San Marino","Vatican")
                          ){
  library(rworldmap) # package to get countries polygon
  europe <- getMap(resolution="high")
  europe <- europe[europe@data$ADMIN.1%in%countries,] # select countries of interest
  europe <-  crop(vect(europe),ext(-10,46,32,72))
  europe <- st_as_sf(europe)
  return(europe)
}

#' Forest cover
#' 
#' @description This is a function to extract forest cover map over europe
#' @note forest cover map .tif file was downloaded at 1x1km resolution from 
#' https://www.eea.europa.eu/data-and-maps/figures/forest-map-of-europe-1
#' @param dir.data directory of project data
#' @param dir.file file path of forest cover
#' @return sf object, polygons of forest cover in europe
#'dir.file="TCD_2018_100m_eu_03035_v020/DATA/TCD_2018_100m_eu_03035_V2_0.tif"
forest_cover <- function(dir.data="data",
                         dir.file,
                         texture){
  forestcover=rast(file.path(dir.data,dir.file))
  forestcover=project(forestcover,subset(rast(texture,crs="epsg:4326"),1),mask=TRUE)
  forestcover=as.data.frame(forestcover,xy=TRUE) %>% 
    mutate(Cover=as.numeric(sub("%.*", "", Class_Name))) %>% 
    select(-Class_Name)
  return(forestcover)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Downloading climatic and pedologic data ####
#' @description Functions used to download/load required climatic and pedologic
#'  data
#' @authors Anne Baranger, Matthieu Combaud (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to compute the water aridity index based on chelsa mean temperature data
#'
#' @param dir.chelsa
#' @param df.loc
#' 
get_waisgdd <- function(dir.chelsa="data/CHELSA/",
                           df.loc){
  pr_chelsa=extract(rast("data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif"),
                    df.loc)[2]
  pet_chelsa=extract(rast("data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif"),
                     df.loc)[2]
  sgdd_chelsa=extract(rast("data/CHELSA/CHELSA_gdd5_1981-2010_V.2.1.tif"),
                      df.loc)[2]
  df.loc=cbind(df.loc,
               pr=pr_chelsa,
               pet=pet_chelsa,
               sgdd=sgdd_chelsa)
  colnames(df.loc)=c("x","y","pr","pet","sgdd")
  df.loc <- df.loc %>% 
    mutate(wai=(pr-12*pet)/pet) 
  return(df.loc)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Computing essential variable over spatial area ####
#' @description Functions used to compute required variables to apply Campbell 
#' equations. Texture from ESDAC and ERA5. SWC daily or monthly timepath.
#' Parameters from Toth.
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load min computed externally with python, from ERA5land
#' 
#' @description Load and reshape SWC data computed with python program, per horizon
#' @note Python code need to be run prior executing this function
#' @param dir.data
#' @param dir.var
#' @param vars list of vars corresponding to each horizons computed
#' @param format either csv, grib or nc
#' @return dataframe with min swc per horizon

load_swc <- function(dir.data="data",
                     dir.file="ERA5-land/daily/swcd-1950-2021-",
                     vars=c("layer1","layer2","layer3","layer4"),
                     extension="",
                     format,
                     europe){
  if (format=="csv"){ # swc daily min computed u
    rast.model=rast(paste0(dir.data,"/",dir.file,vars,extension,".nc"))[[1]]
    rast.swc=rast(nrows=dim(rast.model)[1],
                  ncol=dim(rast.model)[2],
                  xmin=ext(rast.model)$xmin,
                  xmax=ext(rast.model)$xmax,
                  ymin=ext(rast.model)$ymin,
                  ymax=ext(rast.model)$ymax,
                  nlyrs=4)
    for (i in 1:length(vars)){
      rast.min=as.matrix(read.table(paste0(dir.data,"/",dir.file,vars[i],"_min.csv"),
                                    header=FALSE,
                                    sep=",",
                                    dec = "."))
      rast.min[rast.min==-32767] <- NA
      rast.swc[[i]]=rast(nrows=dim(rast.min)[1],
                         ncol=dim(rast.min)[2],
                         xmin=ext(rast.model)$xmin,
                         xmax=ext(rast.model)$xmax,
                         ymin=ext(rast.model)$ymin,
                         ymax=ext(rast.model)$ymax,
                         crs=crs(rast.model),
                         vals=c(t(rast.min)))
      
    }
    names(rast.swc)=paste0("h",1:length(vars))
  }
  
  if (format=="grib"){
    rast.swc=rast(paste0(dir.data,"/",dir.file,vars,extension,".grib"))
    names(rast.swc)=vars
    rast.swc=crop(project(rast.swc,
                          "EPSG:4326"),
                  vect(europe))
  }

  return(as.data.frame(rast.swc,xy=TRUE))
}

#' Compute weighted soil volumetric water content
#' 
#' @description Using ERA5-land SWC and a given depth, the function weight the
#' different horizons according to their thickness
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param SWCtot dataframe of SWCtot cropped to the accurate extent
#' @param depth numeric indicating to which depth swc is to be considered
#' @return dataframe with weighted swc over different horizons
#'  
weight_swc <- function(swc,
                       depth){
  swc<-rast(swc,crs="epsg:4326")
  
  if(depth>100){
    x1=7
    x2=21
    x3=72
    x4=depth-100
  }
  if(depth<=100&depth>28){
    x1=7
    x2=21
    x3=depth-28
    x4=0
  }
  if(depth<=28&depth>7){
    x1=7
    x2=depth-7
    x3=0
    x4=0
  }
  if(depth<=7){
    x1=depth
    x2=0
    x3=0
    x4=0
  }
  swc_weighted <- (x1*swc[[grepl("h1", names(swc))]]+
                     x2*swc[[grepl("h2", names(swc))]]+
                     x3*swc[[grepl("h3", names(swc))]]+
                     x4*swc[[grepl("h4", names(swc))]])/(x1+x2+x3+x4)
  names(swc_weighted)="swc"
  return(as.data.frame(swc_weighted,xy=TRUE))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Computing Psi_min over the spatial area ####
#' @description Functions used to compute psi_min with different
#' methods, using weighted swc or wieghted psi over horizons
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute psi with sureau module
#' 
#' @description function calls swc over different horizons, and use sureau model
#' to integrate their difference in conductivity in the overall psi. It take into 
#' account different possible root depth profile, maximum root depth and parameters
#' of sureau
#' @param swc 5th percentile of swc time serie, for "obs" horizons
#' @param europe Europe spatial extent
#' @param dir.hydro directory of hydraulic parameters
#' @param depth NULL if varying according to ESDAC data, numeric if set to a max
#' @param dir.depth directory of depth raster
#' @param dir.ecoregions directory of WWF ecoregions
#' @param LAImax SUREAU param: max LAI, varying from 2 to 6
#' @param fRootToLeaf SUREAU param: Root to leaf ratio
#' @param rootRadius SUREAU param: root radius 
#' @param beta SUREAU param: root profile, NULL if varying with ecoregions, or
#'  numeric if set constant 
#' @param obs vector of horizons depth of swc file
#' @param ref vector of horizons depth of hydraulic param
##' @param max_depth maximum depth considered, to extend last swc horizon
#' @param file.output file path where to write result
#' @return Psi_min dataframe, that contains values of Psi_min computed by weighting
#' the different horizons


compute_psi_sureau <- function(swc,
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
                               file.output
){
  # Rast model for resolution issues
  rast.mod<-crop(rast(paste0(dir.hydro,
                             "MRC_alp_sl2.tif")),
                 vect(europe),
                 mask=TRUE)
  europe<-vect(europe)
  
  
  # depth used
  print("Load depth")
  if (is.null(depth_max)){ # varying depth
    rast.depth=rast(dir.depth)
    crs(rast.depth) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
    rast.depth <- terra::project(rast.depth,"epsg:4326",method="near")
    rast.depth <- resample(rast.depth,
                           rast.mod,
                           method="near")
    names(rast.depth)="depth"
  } else{
    rast.depth=rast(rast.mod,nlyr=1,vals=depth_max)
    names(rast.depth)="depth"
  }
  
  
  ## Look for biome beta
  print("Compute betarootprofile")
  if (is.null(beta)){
    ecoregions=read_sf(dsn=dir.ecoregions, #read ecoregions
                       layer="wwf_terr_ecos") |> 
      select(BIOME,geometry)
    sf_use_s2(FALSE)
    ecoregions=st_crop(ecoregions,sf::st_bbox(europe))
    points=st_as_sf(swc[,c("x","y")],coords=c("x","y"),crs=st_crs(ecoregions))
    points_biome=st_join(points,ecoregions) |> # associate each ecoR to beta
      mutate(beta=case_when(BIOME%in%c(4,8,9,98)~0.966,
                            BIOME%in%c(5)~0.976,
                            BIOME%in%c(6)~0.943,
                            BIOME==11~0.914,
                            BIOME==12~0.964))
    beta=resample(rast(cbind(swc[,c("x","y")],
                             betaRootProfile=points_biome$beta)),
                  rast.mod,
                  method="near")
    rm(ecoregions,points,points_biome)
  } else{
    beta=resample(rast(cbind(swc[,c("x","y")],
                             betaRootProfile=beta)),
                  rast.mod,
                  method="near")
  }
  
  ## Compute psi and ksoil and ksoiltostem ##
  print("compute SUREAU module")
  # "downscale" SWC 
  swc=resample(rast(swc,crs="epsg:4326"),
               rast.mod,
               method="near")
  
  # Compute root density characteristics per horizons according to depth
  compute_horizon_density <- function(d, betaRootProfile, horizon_start, horizon_end) {
    # Use ifelse for a vectorized condition
    result <- ifelse(d <= horizon_start, 0,
                     switch(
                       as.character(horizon_end),
                       "5" = ifelse(d < 5, 1, (1 - betaRootProfile^5) / (1 - betaRootProfile^d)),
                       ifelse(d<horizon_end,
                              (betaRootProfile^horizon_start - betaRootProfile^d) / (1 - betaRootProfile^d),
                              (betaRootProfile^horizon_start - betaRootProfile^horizon_end) / (1 - betaRootProfile^d)
                       )
                     )
    )
    return(result)
  }
  
  ## Compute SUREAU ##
  RAI = LAImax*fRootToLeaf
  compute.B_GC <- function(hor_value,
                           depth,
                           RAI, 
                           rootRadius,
                           horizon_start,
                           horizon_end
  ) {
    d=ifelse(depth<=horizon_start,
             NA,
             ifelse(depth>horizon_end,
                    horizon_end,
                    depth))
    La <- RAI * hor_value / (2 * pi * rootRadius)
    Lv <- RAI * hor_value / (2 * pi * rootRadius * d)
    b <- 1 / sqrt(pi * Lv)
    B_GC <- La * 2 * pi / log(b / rootRadius)
    return(B_GC)
  }
  
  rast.depth.samp <- c(rast.depth, beta) %>%
    as.data.frame(xy = TRUE, na.rm = FALSE) %>%
    mutate(depth = as.numeric(na_if(depth, 0))) %>%
    mutate(
      hor.1 = compute_horizon_density(d=depth, betaRootProfile, 0, 5),
      hor.2 = compute_horizon_density(d=depth, betaRootProfile, 5, 15),
      hor.3 = compute_horizon_density(d=depth, betaRootProfile, 15, 30),
      hor.4 = compute_horizon_density(d=depth, betaRootProfile, 30, 60),
      hor.5 = compute_horizon_density(d=depth, betaRootProfile, 60, 100),
      hor.6 = compute_horizon_density(d=depth, betaRootProfile, 100, 200),
      hor.7 = compute_horizon_density(d=depth, betaRootProfile, 200, Inf)
    ) |> 
    # compute bcg per horizon
    mutate(hor.1=compute.B_GC(hor.1,depth,RAI,rootRadius,0,5),
           hor.2=compute.B_GC(hor.2,depth,RAI,rootRadius,5,15),
           hor.3=compute.B_GC(hor.3,depth,RAI,rootRadius,15,30),
           hor.4=compute.B_GC(hor.4,depth,RAI,rootRadius,30,60),
           hor.5=compute.B_GC(hor.5,depth,RAI,rootRadius,60,100),
           hor.6=compute.B_GC(hor.6,depth,RAI,rootRadius,100,200),
           hor.7=compute.B_GC(hor.7,depth,RAI,rootRadius,200,Inf)) |> 
    select(-betaRootProfile)
  rast.depth.samp=rast(rast.depth.samp,crs="epsg:4326")
  rm(rast.depth,beta)
  #rast(SWCmin,crs="epsg:4326")
  
  
  # correspondance between horizons fril ERA5 and 3D hydrosoil
  calculate_weights <- function(ref, obs,max_depth=3) {
    
    if(ref[length(ref)]<max_depth){
      ref_ext  <- c(ref,max_depth)
    } else(ref_ext  <- ref)
    if(obs[length(obs)]<ref_ext[length(ref_ext)]){
      obs_ext  <- c(obs,ref_ext[length(ref_ext)])
    } else(obs_ext  <- obs)    
    
    # Create a matrix filled with zeros
    weight_matrix <- matrix(0, nrow =length(obs_ext)-1, ncol = length(ref_ext)-1)
    
    for(i in 1:(length(obs_ext)-1)) {
      for(j in 1:(length(ref_ext)-1)) {
        if(obs_ext[i]<ref_ext[j+1] & ref_ext[j]<obs_ext[i+1]){
          weight_matrix[i,j]=min(1,
                                 (obs_ext[i+1]-ref_ext[j])/(ref_ext[j+1]-ref_ext[j]),
                                 (ref_ext[j+1]-obs_ext[i])/(ref_ext[j+1]-ref_ext[j]),
                                 (obs_ext[i+1]-obs_ext[i])/(ref_ext[j+1]-ref_ext[j]),
                                 na.rm=TRUE)
        }
      }
    }
    if (sum(colSums(weight_matrix))!=length(ref)){print("error")}
    return(weight_matrix)
  }
  
  weight_matrix<-calculate_weights(ref,
                                   obs)
  
  # add a fictive layer if deepest horizon doesn't reach max_depth
  if(obs[length(obs)]<max(ref[length(ref)],max_depth)){
    swc=c(swc,swc[[nlyr(swc)]])
  }
  
  swc_weighted=rast(swc,nlyr=length(ref),vals=NA)
  
  
  print("Computing weighted swc")
  for(i in 1:nlyr(swc_weighted)){
    print(i)
    swc_weighted[[i]]=sum(weight_matrix[,i]*swc,na.rm = TRUE)
  }
  
  
  
  # Gather spatialized PDF parameters with SWC min/max for each of the 7 horizons
  
  print("Croping hydraulic pars")
  pars.files=list.files(dir.hydro)
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  mrc.rast=crop(rast(file.path(dir.hydro,mrc.files)),
                europe)
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  hcc.rast=crop(rast(file.path(dir.hydro,hcc.files)),
                europe)  
  names(swc_weighted)=rep("swc_min",length(ref))
  names(rast.depth.samp)=c("depth",rep("BGC",length(ref)))
  
  
  print("Compute psi per horizon")
  psi_hor=as.data.frame(swc_weighted,xy=TRUE,na.rm=FALSE) |> select(x,y)
  for (i in 1:7){ # loop on the 7 horizons of terraclimate
    print(i)
    psi=c(swc_weighted[[i]],
          rast.depth.samp[[1]], # depth 
          rast.depth.samp[[i+1]], # BCG of horizons
          mrc.rast[[grep(paste0("sl",i),names(mrc.rast))]],
          hcc.rast[[grep(paste0("sl",i),names(hcc.rast))]]
    )
    names(psi)=str_replace(names(psi), paste0("_sl",i), "") #make generic names
    
    psi=as.data.frame(psi,xy=TRUE,na.rm=FALSE)
    setDT(psi)
    psi[, ':=' (
      REW_mrc = ifelse((swc_min - MRC_thr*10^(-4)) / ((MRC_ths - MRC_thr) * 10^(-4)) < 0, 0.01,
                       ifelse((swc_min - MRC_thr*10^(-4)) / ((MRC_ths - MRC_thr) * 10^(-4)) > 1, 1,
                              (swc_min - MRC_thr*10^(-4)) / ((MRC_ths - MRC_thr) * 10^(-4)))),
      REW_hcc = ifelse((swc_min - HCC_thr*10^(-4)) / ((HCC_ths - HCC_thr) * 10^(-4)) < 0, 0.01,
                       ifelse((swc_min - HCC_thr*10^(-4)) / ((HCC_ths - HCC_thr) * 10^(-4)) > 1, 1,
                              (swc_min - HCC_thr*10^(-4)) / ((HCC_ths - HCC_thr) * 10^(-4)))))]
    psi[, ':=' (
      psi_min = -(((1/REW_mrc)^(1/(MRC_m*10^(-4))) - 1)^(1/(MRC_n*10^(-4))) * (1/(MRC_alp*10^(-4))) * 9.78*10^(-2)),
      ksoil = HCC_K0 * (REW_hcc^(HCC_L*10^(-4))) * (1 - (1 - REW_hcc^(1/(HCC_m*10^(-4))))^(HCC_m*10^(-4)))^2,
      ksoilGC = (HCC_K0 * (REW_hcc^(HCC_L*10^(-4))) * (1 - (1 - REW_hcc^(1/(HCC_m*10^(-4))))^(HCC_m*10^(-4)))^2) * 1000 * BGC)]
    psi[, ':=' (
      psi_w = ksoilGC * psi_min
    )]
    colnames(psi)[!colnames(psi) %in% c("x","y","depth")]=paste0(colnames(psi)[!colnames(psi) %in% c("x","y","depth")],"_",i)
    psi_hor=cbind(psi_hor,
                  swc_min=psi$swc_min,
                  ksoilGC=psi$ksoilGC,
                  psi_w=psi$psi_w) |> 
      rename(!!paste0("swc_min_",i):="swc_min",
             !!paste0("ksoilGC_",i):="ksoilGC",
             !!paste0("psi_w_",i):="psi_w") 
    rm(psi)
  }
  rm(hcc.rast,mrc.rast)
  
  print("Compute psi total")
  all_psi_summary = psi_hor %>%
    select(x, y, matches("ksoilGC_|psi_w_")) %>%
    filter(!is.na(psi_w_1)) 
  
  all_psi_summary$sum_psi=rowSums(all_psi_summary[grepl("psi_w_",colnames(all_psi_summary))],na.rm = TRUE)
  all_psi_summary$sum_k=rowSums(all_psi_summary[grepl("ksoilGC_",colnames(all_psi_summary))],na.rm = TRUE)
  all_psi_summary$psi=all_psi_summary$sum_psi/all_psi_summary$sum_k 
  
  fwrite(all_psi_summary,file.output)
  
  return(all_psi_summary)
}


#' Compute psi min taking into account the different horizons, monthly timestep
#' 
#' @description function applies VG equation to compute minimum soil 
#' potential using swc per horizon and associated parameters extracted from 
#' 3D soil database. SWC computed with monthly timestep
#' @param SWCtot dataframe of SWCtot cropped to the accurate extent
#' @param hor number of horizons to account for
#' @param timepath month or day
#' @param date_begin "YYYY-MM-DD" one unit timepath before the beginning of the
#' time serie
#' @return Psi_min dataframe, that contains values of Psi_min computed by weighting
#' the different horizons

compute_psihor <- function(SWCtot,
                           hor,
                           timepath,
                           begin_date,
                           dir.data="data",
                           dir.file="EU_SoilHydroGrids_1km",
                           europe,
                           file.output){
  
  ## Compute SWC extr per horizons ##

  SWCtot <-rast(SWCtot,crs="epsg:4326")
  for (i in 1:hor){
    assign(paste0("SWC_h",i),
           SWCtot[[grepl(paste0("swvl",i), names(SWCtot))]])
  }
  SWC=mget(ls()[grepl(paste0("SWC_h"),ls())])
  rm(list=ls()[grepl("SWC_h",ls())])
  
  # Compute minimum/max per horizon
  SWC_extr=lapply(SWC,function(x)extr_swc(x,timepath,begin_date,europe))
  
  ## Compute SUREAU ##
  # necessary functions 
  kseriesum <- function(k1,k2) {return(1/(1/k1+1/k2))}
  compute.B_GC <- function (La, Lv, rootRadius=0.0004) {
    b <- 1 / sqrt(pi * Lv)
    return(La * 2 * pi / (log(b/rootRadius)))
  }
  
  # compute root characteristics
  betaRootProfile =0.97 
  ## check this
  hor_corr=matrix(data=c(1,2,2,3,3,4,4,1,2,3,4,5,6,7,0,0.05,0.15,0.3,0.6,1,2),ncol=3)
  hor_corr_inv=matrix(data=c(1,2,3,4,1,3,5,7),ncol=2)
  depth <- hor_corr[,3][2:(hor_corr_inv[hor,2]+1)]
  SoilVolume <- depth 
  rootDistribution <-numeric(length(depth)) #Three soil layers
  for (i in 1:(hor_corr_inv[hor,2]-1)){
    rootDistribution[i]=1-betaRootProfile^(depth[i]*100)-sum(rootDistribution[1:(i-1)])
  }
  rootDistribution[hor_corr_inv[hor,2]] = 1-sum(rootDistribution[1:(hor_corr_inv[hor,2]-1)])
  k_RootToStem   = 1 * rootDistribution
  
  fRootToLeaf=1
  LAImax = 5
  RAI = LAImax*fRootToLeaf
  rootRadius=0.0004
  La = RAI*rootDistribution / (2*pi*rootRadius)
  Lv = La/(SoilVolume)
  #compute.B_GC(La,Lv,rootRadius)
  
  ## Compute psi and ksoil and ksoiltostem
  
  # Gather spatialized PDF parameters with SWC min/max for each of the 7 horizons
  pars.files=list.files(file.path(dir.data,dir.file))
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  for (i in 1:hor_corr_inv[hor,2]){
    assign(paste0("psi_h",i),
           c(rast(SWC_extr[[hor_corr[i,1]]]),
             resample(rast(file.path(dir.data,dir.file,mrc.files[endsWith(mrc.files,paste0(i,".tif"))])),
                      rast(SWC_extr[[hor_corr[i,1]]]),
                      method="near"),
             resample(rast(file.path(dir.data,dir.file,hcc.files[endsWith(hcc.files,paste0(i,".tif"))])),
                      rast(SWC_extr[[hor_corr[i,1]]]),
                      method="near")
           )
    )
  }
  list=mget(ls()[grepl("psi_h",ls())])
  
  
  # Compute psi and other relevant variables for each horizons
  psi_min=lapply(seq_along(list),function(i){
    psi_df=list[[i]]
    names(psi_df)=str_replace(names(psi_df), paste0("_sl",i), "")
    psi_df=as.data.frame(psi_df,xy=TRUE) %>% 
      mutate(
        REW_mrc=case_when((SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))<0~0.01,
                          (SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))),
        REW_hcc=case_when((SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))<0~0.01,
                          (SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))),
        psi_min=-(((1/REW_mrc)
                   ^(1/(MRC_m*10^(-4)))
                   -1)
                  ^(1/(MRC_n*10^(-4)))
                  *(1/(MRC_alp*10^(-4)))
                  *9.78*10^(-2)),
        ksoil=(HCC_K0
               *(REW_hcc^(HCC_L*10^(-4)))
               *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                 ^(HCC_m*10^(-4)))^2),
        ksoilGC=(HCC_K0
               *(REW_hcc^(HCC_L*10^(-4)))
               *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                 ^(HCC_m*10^(-4)))^2)*1000 *compute.B_GC(La,Lv,rootRadius)[i],
        k_soiltostem=kseriesum(ksoilGC,k_RootToStem[i]),
        psi_w=ksoilGC*psi_min
      )
    
    return(psi_df)
  }
  )
  
  plot_psi=psi_min[[1]][,1:2]
  sum_psi=numeric(dim(psi_min[[1]])[1])
  sum_k=numeric(dim(psi_min[[1]])[1])
  for(i in 1:length(psi_min)){
    plot_psi=cbind(plot_psi,psi_min[[i]][,c(3,4,17:19)])
    names(plot_psi)=c(names(plot_psi)[1:(2+5*(i-1))],paste0(names(psi_min[[i]][,c(3,4,17:19)]),i))
    sum_psi=sum_psi+psi_min[[i]]$psi_w
    sum_k=sum_k+psi_min[[i]]$ksoilGC
  }
  plot_psi$psi=sum_psi/sum_k 
  write.csv(plot_psi,file.output)
  return(plot_psi)
}


#' Compute psi min weighted over horizons with daily timepath
#' 
#' @description function applies VG equation to compute minimum soil 
#' potential using swc per horizon and associated parameters extracted from 
#' 3D soil database. SWC computed with daily timestep
#' @param SWCmin dataframe of SWC min of each horizons
#' @param hor number of horizons to account for
#' @param europe Europe polygon
#' @param dir.data
#' @param dir.file directory of 3D soil database
#' @return Psi_min dataframe, that contains values of Psi_min computed with SUREAU 
#' method
compute_psihorday <- function(SWCmin,
                              europe,
                              dir.data="data",
                              dir.file="EU_SoilHydroGrids_1km",
                              depth,
                              file.output){
  
  ## Compute SUREAU ##
  # necessary functions 
  kseriesum <- function(k1,k2) {return(1/(1/k1+1/k2))}
  compute.B_GC <- function (La, Lv, rootRadius=0.0004) {
    b <- 1 / sqrt(pi * Lv)
    return(La * 2 * pi / (log(b/rootRadius)))
  }
  
  # compute root characteristics
  betaRootProfile =0.97 
  fRootToLeaf=1
  LAImax = 5
  RAI = LAImax*fRootToLeaf
  rootRadius=0.0004
  
  # depth used
  if (typeof(depth)=="character"){
    rast.depth=rast("data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst")
    crs(rast.depth) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
    rast.depth <- terra::project(rast.depth,"epsg:4326",method="near")
    rast.depth <- resample(rast.depth,
                           crop(rast(file.path(dir.data,
                                               dir.file,
                                               "MRC_alp_sl2.tif")),
                                vect(europe),
                                mask=TRUE),
                           method="near")
    names(rast.depth)="depth"
  } else {
    # create a raster model
    rast.depth=crop(rast(file.path(dir.data,
                                   dir.file,
                                   "MRC_alp_sl2.tif")),
                    vect(europe),
                    mask=TRUE)
    rast.depth[[1]][rast.depth[[1]]>=0] <- depth
    names(rast.depth)="depth"
  }
 
 
  ## Compute psi and ksoil and ksoiltostem ##
  
  # correspondance between horizons fril ERA5 and 3D hydrosoil
  hor_corr=matrix(data=c(1,2,2,3,3,4,4,1,2,3,4,5,6,7,0,0.05,0.15,0.3,0.6,1,2),ncol=3)
  hor_corr_inv=matrix(data=c(1,2,3,4,1,3,5,7),ncol=2)
  
  # "downscale" SWC 
  SWCmin=resample(rast(SWCmin,crs="epsg:4326"),
                  crop(rast(file.path(dir.data,
                                 dir.file,
                                 "MRC_alp_sl2.tif")),
                            vect(europe),
                            mask=TRUE),
                  method="near")
  
  # Compute root density characteristics per horizons according to depth
  rast.depth.samp <- rast.depth %>% 
    as.data.frame(xy=TRUE,na.rm=FALSE) %>% 
    #filter(depth!=0)%>% 
    ## Compute rescaled density of roots per horizons
    mutate(depth=as.numeric(depth),
           depth=na_if(depth,0),
           hor.1=case_when(depth<5~1,
                           TRUE~(1-0.97^5)/(1-0.97^depth)),
           hor.2=case_when(depth<=5~0,
                           depth>5 & depth < 15 ~ (0.97^5-0.97^depth)/(1-0.97^depth),
                           TRUE~(0.97^5-0.97^15)/(1-0.97^depth)),
           hor.3=case_when(depth<=15~0,
                           depth>15 & depth < 30 ~ (0.97^15-0.97^depth)/(1-0.97^depth),
                           TRUE~(0.97^15-0.97^30)/(1-0.97^depth)),
           hor.4=case_when(depth<=30~0,
                           depth>30 & depth < 60 ~ (0.97^30-0.97^depth)/(1-0.97^depth),
                           TRUE~(0.97^30-0.97^60)/(1-0.97^depth)),
           hor.5=case_when(depth<=60~0,
                           depth>60 & depth < 100 ~ (0.97^60-0.97^depth)/(1-0.97^depth),
                           TRUE~(0.97^60-0.97^100)/(1-0.97^depth)),
           hor.6=case_when(depth<=100~0,
                           depth>100 & depth < 200 ~ (0.97^100-0.97^depth)/(1-0.97^depth),
                           TRUE~(0.97^100-0.97^200)/(1-0.97^depth)),
           hor.7=case_when(depth<=200~0,
                           TRUE~(0.97^200-0.97^depth)/(1-0.97^depth))
    ) %>%
    # mutate(depth=as.numeric(depth),
    #        hor.1=case_when(depth<5~depth,
    #                        TRUE~5),
    #        hor.2=case_when(depth<=5~0,
    #                        depth>5 & depth < 15 ~ depth-5,
    #                        TRUE~10),
    #        hor.3=case_when(depth<=15~0,
    #                        depth>15 & depth < 30 ~ depth-15,
    #                        TRUE~15),
    #        hor.4=case_when(depth<=30~0,
    #                        depth>30 & depth < 60 ~ depth-30,
    #                        TRUE~30),
    #        hor.5=case_when(depth<=60~0,
    #                        depth>60 & depth < 100 ~ depth-60,
    #                        TRUE~40),
    #        hor.6=case_when(depth<=100~0,
    #                        depth>100 & depth < 200 ~ depth-100,
    #                        TRUE~100),
    #        hor.7=case_when(depth<=200~0,
    #                        TRUE~depth-200)
    # ) %>% 
    
    ## Compute  BCG per horizons
    mutate(hor.1=case_when(depth<5~compute.B_GC(La = RAI*hor.1/(2*pi*rootRadius),
                                                Lv = RAI*hor.1/(2*pi*rootRadius*depth)),
                           TRUE~compute.B_GC(La = RAI*hor.1/(2*pi*rootRadius),
                                             Lv = RAI*hor.1/(2*pi*rootRadius*5))),
           hor.2=case_when(depth<=5~0,
                           depth>5 & depth < 15 ~compute.B_GC(La = RAI*hor.2/(2*pi*rootRadius),
                                                              Lv = RAI*hor.2/(2*pi*rootRadius*(depth-5))),
                           TRUE~compute.B_GC(La = RAI*hor.2/(2*pi*rootRadius),
                                             Lv = RAI*hor.2/(2*pi*rootRadius*10))),
           hor.3=case_when(depth<=15~0,
                           depth>15 & depth < 30 ~ compute.B_GC(La = RAI*hor.3/(2*pi*rootRadius),
                                                                Lv = RAI*hor.3/(2*pi*rootRadius*(depth-15))),
                           TRUE~compute.B_GC(La = RAI*hor.3/(2*pi*rootRadius),
                                             Lv = RAI*hor.3/(2*pi*rootRadius*15))),
           hor.4=case_when(depth<=30~0,
                           depth>30 & depth < 60 ~ compute.B_GC(La = RAI*hor.4/(2*pi*rootRadius),
                                                                Lv = RAI*hor.4/(2*pi*rootRadius*(depth-30))),
                           TRUE~compute.B_GC(La = RAI*hor.4/(2*pi*rootRadius),
                                             Lv = RAI*hor.4/(2*pi*rootRadius*30))),
           hor.5=case_when(depth<=60~0,
                           depth>60 & depth < 100 ~compute.B_GC(La = RAI*hor.5/(2*pi*rootRadius),
                                                                Lv = RAI*hor.5/(2*pi*rootRadius*(depth-60))),
                           TRUE~compute.B_GC(La = RAI*hor.5/(2*pi*rootRadius),
                                             Lv = RAI*hor.5/(2*pi*rootRadius*40))),
           hor.6=case_when(depth<=100~0,
                           depth>100 & depth < 200 ~compute.B_GC(La = RAI*hor.6/(2*pi*rootRadius),
                                                                 Lv = RAI*hor.6/(2*pi*rootRadius*(depth-100))),
                           TRUE~compute.B_GC(La = RAI*hor.6/(2*pi*rootRadius),
                                             Lv = RAI*hor.6/(2*pi*rootRadius*100))),
           hor.7=case_when(depth<=200~0,
                           TRUE~compute.B_GC(La = RAI*hor.7/(2*pi*rootRadius),
                                             Lv = RAI*hor.7/(2*pi*rootRadius*(depth-200))))
    )
  rast.depth.samp=rast(rast.depth.samp,crs="epsg:4326")
  #rast(SWCmin,crs="epsg:4326")
  
  # Gather spatialized PDF parameters with SWC min/max for each of the 7 horizons
  names(SWCmin)=rep("SWC_min",4)
  names(rast.depth.samp)=c("depth",rep("BGC",7))
  pars.files=list.files(file.path(dir.data,dir.file))
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  for (i in 1:7){ # loop on the 7 horizons
    psi=c(SWCmin[[hor_corr[i,1]]],
          rast.depth.samp[[1]], # depth 
          rast.depth.samp[[i+1]], # BCG of horizons
          crop(rast(file.path(dir.data,dir.file,mrc.files[endsWith(mrc.files,paste0(i,".tif"))])),
               vect(europe),
               mask=TRUE),
          crop(rast(file.path(dir.data,dir.file,hcc.files[endsWith(hcc.files,paste0(i,".tif"))])),
               vect(europe),
               mask=TRUE)
             )
    names(psi)=str_replace(names(psi), paste0("_sl",i), "") #make generic names
    psi=as.data.frame(psi,xy=TRUE) %>% 
      mutate(
        REW_mrc=case_when((SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))<0~0.01,
                          (SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))),
        REW_hcc=case_when((SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))<0~0.01,
                          (SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))),
        psi_min=-(((1/REW_mrc)
                   ^(1/(MRC_m*10^(-4)))
                   -1)
                  ^(1/(MRC_n*10^(-4)))
                  *(1/(MRC_alp*10^(-4)))
                  *9.78*10^(-2)),
        ksoil=(HCC_K0
               *(REW_hcc^(HCC_L*10^(-4)))
               *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                 ^(HCC_m*10^(-4)))^2),
        ksoilGC=(HCC_K0
                 *(REW_hcc^(HCC_L*10^(-4)))
                 *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                   ^(HCC_m*10^(-4)))^2)*1000 *BGC,
        #k_soiltostem=kseriesum(ksoilGC,k_RootToStem[i]),
        psi_w=ksoilGC*psi_min
      )
    assign(paste0("psi_h",i),
           psi)
    rm(psi)
  }
  
  
  # Apply SUREAU wieghting method across horizons
  psi_min=psi_h1[c("x","y")]
  for (i in 1:7){
    psi_min=psi_min %>% 
      left_join(eval(parse(text=paste0("psi_h",i)))[c("x","y","ksoilGC","psi_w")],
                by=c("x","y"))
    colnames(psi_min)=c(colnames(psi_min)[1:(2+2*(i-1))],
                        paste0("ksoilGC_",i),
                        paste0("psi_w_",i))
  }
  sum_psi=rowSums(psi_min[grepl("psi_w_",colnames(psi_min))],na.rm = TRUE)
  sum_k=rowSums(psi_min[grepl("ksoilGC_",colnames(psi_min))],na.rm = TRUE)
  psi_min$psi=sum_psi/sum_k 
  
  write.csv(psi_min,file.output)

  return(psi_min)
}



#' Compute psi min weighted over horizons with daily timepath, and with varying beta
#' 
#' @description function applies VG equation to compute minimum soil 
#' potential using swc per horizon and associated parameters extracted from 
#' 3D soil database. SWC computed with daily timestep
#' @param SWCmin dataframe of SWC min of each horizons
#' @param hor number of horizons to account for
#' @param europe Europe polygon
#' @param dir.data
#' @param dir.file directory of 3D soil database
#' @return Psi_min dataframe, that contains values of Psi_min computed with SUREAU 
#' method
compute_psihorday_beta <- function(SWCmin,
                              europe,
                              dir.data="data",
                              dir.file="EU_SoilHydroGrids_1km",
                              LAImax=5,
                              file.output
                              ){
  
  ## Compute SUREAU ##
  # necessary functions 
  kseriesum <- function(k1,k2) {return(1/(1/k1+1/k2))}
  compute.B_GC <- function (La, Lv, rootRadius=0.0004) {
    b <- 1 / sqrt(pi * Lv)
    return(La * 2 * pi / (log(b/rootRadius)))
  }
  
  # compute root characteristics
  # betaRootProfile =0.97 
  fRootToLeaf=1
  # LAImax = 5
  RAI = LAImax*fRootToLeaf
  rootRadius=0.0004
  
  # depth used
  rast.depth=rast("data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst")
  crs(rast.depth) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  rast.depth <- terra::project(rast.depth,"epsg:4326",method="near")
  rast.depth <- resample(rast.depth,
                         crop(rast(file.path(dir.data,
                                             dir.file,
                                             "MRC_alp_sl2.tif")),
                              vect(europe),
                              mask=TRUE),
                         method="near")
  names(rast.depth)="depth"
  
  ## Look for biome beta
  ecoregions=read_sf(dsn="data/WWF/official",
                     layer="wwf_terr_ecos") |> 
    select(BIOME,geometry)
  sf_use_s2(FALSE)
  ecoregions=st_crop(ecoregions,sf::st_bbox(europe))
  points=st_as_sf(SWCmin[,c("x","y")],coords=c("x","y"),crs=st_crs(ecoregions))
  points_biome=st_join(points,ecoregions) |> 
    mutate(beta=case_when(BIOME%in%c(4,8,9,98)~0.966,
                          BIOME%in%c(5)~0.976,
                          BIOME%in%c(6)~0.943,
                          BIOME==11~0.914,
                          BIOME==12~0.964))
  beta=resample(rast(cbind(SWCmin[,c("x","y")],
                           betaRootProfile=points_biome$beta)),
                crop(rast(file.path(dir.data,
                                    dir.file,
                                    "MRC_alp_sl2.tif")),
                     vect(europe),
                     mask=TRUE),
                method="near")
  ## Compute psi and ksoil and ksoiltostem ##
  
  # correspondance between horizons fril ERA5 and 3D hydrosoil
  hor_corr=matrix(data=c(1,2,2,3,3,4,4,1,2,3,4,5,6,7,0,0.05,0.15,0.3,0.6,1,2),ncol=3)
  hor_corr_inv=matrix(data=c(1,2,3,4,1,3,5,7),ncol=2)
  
  # "downscale" SWC 
  SWCmin=resample(rast(SWCmin,crs="epsg:4326"),
                  crop(rast(file.path(dir.data,
                                      dir.file,
                                      "MRC_alp_sl2.tif")),
                       vect(europe),
                       mask=TRUE),
                  method="near")
  
  # Compute root density characteristics per horizons according to depth
  rast.depth.samp <- c(rast.depth,beta) %>% 
    as.data.frame(xy=TRUE,na.rm=FALSE) %>% 
    #filter(depth!=0)%>% 
    ## Compute rescaled density of roots per horizons
    mutate(depth=as.numeric(depth),
           depth=na_if(depth,0),
           hor.1=case_when(depth<5~1,
                           TRUE~(1-betaRootProfile^5)/(1-betaRootProfile^depth)),
           hor.2=case_when(depth<=5~0,
                           depth>5 & depth < 15 ~ (betaRootProfile^5-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^5-betaRootProfile^15)/(1-betaRootProfile^depth)),
           hor.3=case_when(depth<=15~0,
                           depth>15 & depth < 30 ~ (betaRootProfile^15-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^15-betaRootProfile^30)/(1-betaRootProfile^depth)),
           hor.4=case_when(depth<=30~0,
                           depth>30 & depth < 60 ~ (betaRootProfile^30-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^30-betaRootProfile^60)/(1-betaRootProfile^depth)),
           hor.5=case_when(depth<=60~0,
                           depth>60 & depth < 100 ~ (betaRootProfile^60-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^60-betaRootProfile^100)/(1-betaRootProfile^depth)),
           hor.6=case_when(depth<=100~0,
                           depth>100 & depth < 200 ~ (betaRootProfile^100-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^100-betaRootProfile^200)/(1-betaRootProfile^depth)),
           hor.7=case_when(depth<=200~0,
                           TRUE~(betaRootProfile^200-betaRootProfile^depth)/(1-betaRootProfile^depth))
    ) %>%
  ## Compute  BCG per horizons
  mutate(hor.1=case_when(depth<5~compute.B_GC(La = RAI*hor.1/(2*pi*rootRadius),
                                              Lv = RAI*hor.1/(2*pi*rootRadius*depth)),
                         TRUE~compute.B_GC(La = RAI*hor.1/(2*pi*rootRadius),
                                           Lv = RAI*hor.1/(2*pi*rootRadius*5))),
         hor.2=case_when(depth<=5~0,
                         depth>5 & depth < 15 ~compute.B_GC(La = RAI*hor.2/(2*pi*rootRadius),
                                                            Lv = RAI*hor.2/(2*pi*rootRadius*(depth-5))),
                         TRUE~compute.B_GC(La = RAI*hor.2/(2*pi*rootRadius),
                                           Lv = RAI*hor.2/(2*pi*rootRadius*10))),
         hor.3=case_when(depth<=15~0,
                         depth>15 & depth < 30 ~ compute.B_GC(La = RAI*hor.3/(2*pi*rootRadius),
                                                              Lv = RAI*hor.3/(2*pi*rootRadius*(depth-15))),
                         TRUE~compute.B_GC(La = RAI*hor.3/(2*pi*rootRadius),
                                           Lv = RAI*hor.3/(2*pi*rootRadius*15))),
         hor.4=case_when(depth<=30~0,
                         depth>30 & depth < 60 ~ compute.B_GC(La = RAI*hor.4/(2*pi*rootRadius),
                                                              Lv = RAI*hor.4/(2*pi*rootRadius*(depth-30))),
                         TRUE~compute.B_GC(La = RAI*hor.4/(2*pi*rootRadius),
                                           Lv = RAI*hor.4/(2*pi*rootRadius*30))),
         hor.5=case_when(depth<=60~0,
                         depth>60 & depth < 100 ~compute.B_GC(La = RAI*hor.5/(2*pi*rootRadius),
                                                              Lv = RAI*hor.5/(2*pi*rootRadius*(depth-60))),
                         TRUE~compute.B_GC(La = RAI*hor.5/(2*pi*rootRadius),
                                           Lv = RAI*hor.5/(2*pi*rootRadius*40))),
         hor.6=case_when(depth<=100~0,
                         depth>100 & depth < 200 ~compute.B_GC(La = RAI*hor.6/(2*pi*rootRadius),
                                                               Lv = RAI*hor.6/(2*pi*rootRadius*(depth-100))),
                         TRUE~compute.B_GC(La = RAI*hor.6/(2*pi*rootRadius),
                                           Lv = RAI*hor.6/(2*pi*rootRadius*100))),
         hor.7=case_when(depth<=200~0,
                         TRUE~compute.B_GC(La = RAI*hor.7/(2*pi*rootRadius),
                                           Lv = RAI*hor.7/(2*pi*rootRadius*(depth-200))))
  ) |> 
    select(-betaRootProfile)
  rast.depth.samp=rast(rast.depth.samp,crs="epsg:4326")
  #rast(SWCmin,crs="epsg:4326")
  
  # Gather spatialized PDF parameters with SWC min/max for each of the 7 horizons
  names(SWCmin)=rep("SWC_min",4)
  names(rast.depth.samp)=c("depth",rep("BGC",7))
  pars.files=list.files(file.path(dir.data,dir.file))
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  for (i in 1:7){ # loop on the 7 horizons
    psi=c(SWCmin[[hor_corr[i,1]]],
          rast.depth.samp[[1]], # depth 
          rast.depth.samp[[i+1]], # BCG of horizons
          crop(rast(file.path(dir.data,dir.file,mrc.files[endsWith(mrc.files,paste0(i,".tif"))])),
               vect(europe),
               mask=TRUE),
          crop(rast(file.path(dir.data,dir.file,hcc.files[endsWith(hcc.files,paste0(i,".tif"))])),
               vect(europe),
               mask=TRUE)
    )
    names(psi)=str_replace(names(psi), paste0("_sl",i), "") #make generic names
    psi=as.data.frame(psi,xy=TRUE) %>% 
      mutate(
        REW_mrc=case_when((SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))<0~0.01,
                          (SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))),
        REW_hcc=case_when((SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))<0~0.01,
                          (SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))),
        psi_min=-(((1/REW_mrc)
                   ^(1/(MRC_m*10^(-4)))
                   -1)
                  ^(1/(MRC_n*10^(-4)))
                  *(1/(MRC_alp*10^(-4)))
                  *9.78*10^(-2)),
        ksoil=(HCC_K0
               *(REW_hcc^(HCC_L*10^(-4)))
               *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                 ^(HCC_m*10^(-4)))^2),
        ksoilGC=(HCC_K0
                 *(REW_hcc^(HCC_L*10^(-4)))
                 *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                   ^(HCC_m*10^(-4)))^2)*1000 *BGC,
        #k_soiltostem=kseriesum(ksoilGC,k_RootToStem[i]),
        psi_w=ksoilGC*psi_min
      )
    assign(paste0("psi_h",i),
           psi)
    rm(psi)
  }
  
  
  # Apply SUREAU wieghting method across horizons
  psi_min=psi_h1[c("x","y")]
  for (i in 1:7){
    psi_min=psi_min %>% 
      left_join(eval(parse(text=paste0("psi_h",i)))[c("x","y","ksoilGC","psi_w")],
                by=c("x","y"))
    colnames(psi_min)=c(colnames(psi_min)[1:(2+2*(i-1))],
                        paste0("ksoilGC_",i),
                        paste0("psi_w_",i))
  }
  sum_psi=rowSums(psi_min[grepl("psi_w_",colnames(psi_min))],na.rm = TRUE)
  sum_k=rowSums(psi_min[grepl("ksoilGC_",colnames(psi_min))],na.rm = TRUE)
  psi_min$psi=sum_psi/sum_k 
  
  fwrite(psi_min,file.output)
  
  return(psi_min)
}



#' Compute psi min weighted over horizons with daily timepath
#' 
#' @description function applies VG equation to compute minimum soil 
#' potential using swc per horizon and associated parameters extracted from 
#' 3D soil database. SWC computed with daily timestep
#' @param SWCmin dataframe of SWC min of each horizons
#' @param hor number of horizons to account for
#' @param europe Europe polygon
#' @param dir.data
#' @param dir.file directory of 3D soil database
#' @return Psi_min dataframe, that contains values of Psi_min computed with SUREAU 
#' method
compute_sensitivity <- function(SWCmin,
                                europe,
                                dir.data="data",
                                dir.file="EU_SoilHydroGrids_1km",
                                depth,
                                LAImax,
                                betaRootProfile){
  
  ## Compute SUREAU ##
  # necessary functions 
  kseriesum <- function(k1,k2) {return(1/(1/k1+1/k2))}
  compute.B_GC <- function (La, Lv, rootRadius=0.0004) {
    b <- 1 / sqrt(pi * Lv)
    return(La * 2 * pi / (log(b/rootRadius)))
  }
  
  # compute root characteristics
  # betaRootProfile =0.97 
  fRootToLeaf=1
  # LAImax = 5
  RAI = LAImax*fRootToLeaf
  rootRadius=0.0004
  
  # depth used
  if (typeof(depth)=="character"){
    rast.depth=rast("data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst")
    crs(rast.depth) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
    rast.depth <- terra::project(rast.depth,"epsg:4326",method="near")
    rast.depth <- resample(rast.depth,
                           crop(rast(file.path(dir.data,
                                               dir.file,
                                               "MRC_alp_sl2.tif")),
                                vect(europe),
                                mask=TRUE),
                           method="near")
    names(rast.depth)="depth"
  } else {
    # create a raster model
    rast.depth=crop(rast(file.path(dir.data,
                                   dir.file,
                                   "MRC_alp_sl2.tif")),
                    vect(europe),
                    mask=TRUE)
    rast.depth[[1]][rast.depth[[1]]>=0] <- depth
    names(rast.depth)="depth"
  }
  
  
  ## Compute psi and ksoil and ksoiltostem ##
  
  # correspondance between horizons fril ERA5 and 3D hydrosoil
  hor_corr=matrix(data=c(1,2,2,3,3,4,4,1,2,3,4,5,6,7,0,0.05,0.15,0.3,0.6,1,2),ncol=3)
  hor_corr_inv=matrix(data=c(1,2,3,4,1,3,5,7),ncol=2)
  
  # "downscale" SWC 
  SWCmin=resample(rast(SWCmin,crs="epsg:4326"),
                  crop(rast(file.path(dir.data,
                                      dir.file,
                                      "MRC_alp_sl2.tif")),
                       vect(europe),
                       mask=TRUE),
                  method="near")
  
  # Compute root density characteristics per horizons according to depth
  rast.depth.samp <- rast.depth %>% 
    as.data.frame(xy=TRUE,na.rm=FALSE) %>% 
    #filter(depth!=0)%>% 
    ## Compute rescaled density of roots per horizons
    mutate(depth=as.numeric(depth),
           depth=na_if(depth,0),
           hor.1=case_when(depth<5~1,
                           TRUE~(1-betaRootProfile^5)/(1-betaRootProfile^depth)),
           hor.2=case_when(depth<=5~0,
                           depth>5 & depth < 15 ~ (betaRootProfile^5-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^5-betaRootProfile^15)/(1-betaRootProfile^depth)),
           hor.3=case_when(depth<=15~0,
                           depth>15 & depth < 30 ~ (betaRootProfile^15-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^15-betaRootProfile^30)/(1-betaRootProfile^depth)),
           hor.4=case_when(depth<=30~0,
                           depth>30 & depth < 60 ~ (betaRootProfile^30-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^30-betaRootProfile^60)/(1-betaRootProfile^depth)),
           hor.5=case_when(depth<=60~0,
                           depth>60 & depth < 100 ~ (betaRootProfile^60-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^60-betaRootProfile^100)/(1-betaRootProfile^depth)),
           hor.6=case_when(depth<=100~0,
                           depth>100 & depth < 200 ~ (betaRootProfile^100-betaRootProfile^depth)/(1-betaRootProfile^depth),
                           TRUE~(betaRootProfile^100-betaRootProfile^200)/(1-betaRootProfile^depth)),
           hor.7=case_when(depth<=200~0,
                           TRUE~(betaRootProfile^200-betaRootProfile^depth)/(1-betaRootProfile^depth))
    ) %>%
    ## Compute  BCG per horizons
    mutate(hor.1=case_when(depth<5~compute.B_GC(La = RAI*hor.1/(2*pi*rootRadius),
                                                Lv = RAI*hor.1/(2*pi*rootRadius*depth)),
                           TRUE~compute.B_GC(La = RAI*hor.1/(2*pi*rootRadius),
                                             Lv = RAI*hor.1/(2*pi*rootRadius*5))),
           hor.2=case_when(depth<=5~0,
                           depth>5 & depth < 15 ~compute.B_GC(La = RAI*hor.2/(2*pi*rootRadius),
                                                              Lv = RAI*hor.2/(2*pi*rootRadius*(depth-5))),
                           TRUE~compute.B_GC(La = RAI*hor.2/(2*pi*rootRadius),
                                             Lv = RAI*hor.2/(2*pi*rootRadius*10))),
           hor.3=case_when(depth<=15~0,
                           depth>15 & depth < 30 ~ compute.B_GC(La = RAI*hor.3/(2*pi*rootRadius),
                                                                Lv = RAI*hor.3/(2*pi*rootRadius*(depth-15))),
                           TRUE~compute.B_GC(La = RAI*hor.3/(2*pi*rootRadius),
                                             Lv = RAI*hor.3/(2*pi*rootRadius*15))),
           hor.4=case_when(depth<=30~0,
                           depth>30 & depth < 60 ~ compute.B_GC(La = RAI*hor.4/(2*pi*rootRadius),
                                                                Lv = RAI*hor.4/(2*pi*rootRadius*(depth-30))),
                           TRUE~compute.B_GC(La = RAI*hor.4/(2*pi*rootRadius),
                                             Lv = RAI*hor.4/(2*pi*rootRadius*30))),
           hor.5=case_when(depth<=60~0,
                           depth>60 & depth < 100 ~compute.B_GC(La = RAI*hor.5/(2*pi*rootRadius),
                                                                Lv = RAI*hor.5/(2*pi*rootRadius*(depth-60))),
                           TRUE~compute.B_GC(La = RAI*hor.5/(2*pi*rootRadius),
                                             Lv = RAI*hor.5/(2*pi*rootRadius*40))),
           hor.6=case_when(depth<=100~0,
                           depth>100 & depth < 200 ~compute.B_GC(La = RAI*hor.6/(2*pi*rootRadius),
                                                                 Lv = RAI*hor.6/(2*pi*rootRadius*(depth-100))),
                           TRUE~compute.B_GC(La = RAI*hor.6/(2*pi*rootRadius),
                                             Lv = RAI*hor.6/(2*pi*rootRadius*100))),
           hor.7=case_when(depth<=200~0,
                           TRUE~compute.B_GC(La = RAI*hor.7/(2*pi*rootRadius),
                                             Lv = RAI*hor.7/(2*pi*rootRadius*(depth-200))))
    )
  rast.depth.samp=rast(rast.depth.samp,crs="epsg:4326")
  #rast(SWCmin,crs="epsg:4326")
  
  # Gather spatialized PDF parameters with SWC min/max for each of the 7 horizons
  names(SWCmin)=rep("SWC_min",4)
  names(rast.depth.samp)=c("depth",rep("BGC",7))
  pars.files=list.files(file.path(dir.data,dir.file))
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  for (i in 1:7){ # loop on the 7 horizons
    psi=c(SWCmin[[hor_corr[i,1]]],
          rast.depth.samp[[1]], # depth 
          rast.depth.samp[[i+1]], # BCG of horizons
          crop(rast(file.path(dir.data,dir.file,mrc.files[endsWith(mrc.files,paste0(i,".tif"))])),
               vect(europe),
               mask=TRUE),
          crop(rast(file.path(dir.data,dir.file,hcc.files[endsWith(hcc.files,paste0(i,".tif"))])),
               vect(europe),
               mask=TRUE)
    )
    names(psi)=str_replace(names(psi), paste0("_sl",i), "") #make generic names
    psi=as.data.frame(psi,xy=TRUE) %>% 
      mutate(
        REW_mrc=case_when((SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))<0~0.01,
                          (SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-MRC_thr*10^(-4))/((MRC_ths-MRC_thr)*10^(-4))),
        REW_hcc=case_when((SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))<0~0.01,
                          (SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))>1~1,
                          TRUE~(SWC_min-HCC_thr*10^(-4))/((HCC_ths-HCC_thr)*10^(-4))),
        psi_min=-(((1/REW_mrc)
                   ^(1/(MRC_m*10^(-4)))
                   -1)
                  ^(1/(MRC_n*10^(-4)))
                  *(1/(MRC_alp*10^(-4)))
                  *9.78*10^(-2)),
        ksoil=(HCC_K0
               *(REW_hcc^(HCC_L*10^(-4)))
               *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                 ^(HCC_m*10^(-4)))^2),
        ksoilGC=(HCC_K0
                 *(REW_hcc^(HCC_L*10^(-4)))
                 *(1-(1-REW_hcc^(1/(HCC_m*10^(-4))))
                   ^(HCC_m*10^(-4)))^2)*1000 *BGC,
        #k_soiltostem=kseriesum(ksoilGC,k_RootToStem[i]),
        psi_w=ksoilGC*psi_min
      )
    assign(paste0("psi_h",i),
           psi)
    rm(psi)
  }
  
  
  # Apply SUREAU wieghting method across horizons
  psi_min=psi_h1[c("x","y")]
  for (i in 1:7){
    psi_min=psi_min %>% 
      left_join(eval(parse(text=paste0("psi_h",i)))[c("x","y","ksoilGC","psi_w")],
                by=c("x","y"))
    colnames(psi_min)=c(colnames(psi_min)[1:(2+2*(i-1))],
                        paste0("ksoilGC_",i),
                        paste0("psi_w_",i))
  }
  sum_psi=rowSums(psi_min[grepl("psi_w_",colnames(psi_min))],na.rm = TRUE)
  sum_k=rowSums(psi_min[grepl("ksoilGC_",colnames(psi_min))],na.rm = TRUE)
  psi_min$psi=sum_psi/sum_k 
  
  # write.csv(psi_min,file.output)
  
  return(psi_min)
}






#' Mask psi_min with only forest areas
#' 
#' @description function that enables to mask psi_min data with only forest 
#' areas according to a forest cover threshold
#' @param psi_min psi values dataset
#' @param forestcover forestcover data.frame, cleaned 
#' @param ts threshold of forest cover to be applied
#' @return Psi_min_masked dataframe, that contains values of Psi_min over areas
#' with more then "ts" forestcover
#' 
psi_forest <- function(psi_min,
                       forestcover,
                       ts){
  forestcover=forestcover %>% 
    filter(Cover>ts)
  psi_min=rast(psi_min,crs="epsg:4326")
  forestcover=project(rast(forestcover,crs="epsg:4326"),subset(psi_min,1))
  psi_forest=as.data.frame(mask(psi_min,forestcover),xy=TRUE)
  return(psi_forest)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 5 - Pedo-transfer functions ####
#' @description Functions used to explore parametric pedotransfer functions
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Compute Campbell PDF parameters
#' 
#' @description Computing the paramter psi_e and b of campbell pedotransfer 
#' function based on soil composition in clay, silt, sand
#' @param clay_ct [0-1] mass fraction of clay 
#' @param silt_ct [0-1] mass fraction of silt
#' @param sand_ct [0-1] mass fraction of sand
#' @return psi_e and b parameters
#' 
campbell_par <- function(clay_ct,silt_ct,sand_ct){
  a=clay_ct*log(0.001)+silt_ct*log(0.026)+sand_ct*log(1.025)
  b=((clay_ct*log(0.001)^2+silt_ct*log(0.026)^2+sand_ct*log(1.025)^2)-a^2)^(1/2)
  GMD=exp(a)
  GSD=exp(b)
  psi_e=-0.5*GMD^(-0.5)
  b=-2*psi_e+0.2*GSD
  return(data.frame(psi_e,b))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 6 - Compute frost safety margins ####
#' @description From date of growing season, LT50 and Tmin time series, compute
#' the safety margin as the difference between LT50 estimated at budburst and 
#' averaged minimum temperature before budburst
#' @note budburst date are from 
#' https://www.eea.europa.eu/data-and-maps/data/annual-start-of-vegetation-growing
#' @note other link for budburst
#'  https://land.copernicus.eu/pan-european/biophysical-parameters/high-resolution-vegetation-phenology-and-productivity/vegetation-indices
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute frost index
#' 
#' @description Compute the averaged minimum temperature during the 30 days beofrdffCompute frost according to budburst date
#' @param europe df of psimin computed over europe
#' @param year_date beginning year of budburst timeseries
#' @return dataframe of frost index
#'

get.frostindex <- function(europe,
                           dir.clim, #="data/CHELSA/CHELSA_EUR11_tasmin_month_quant05_19802005.nc",
                           dir.fdg, #="data/eea_2000-2016/SOS_2000_2016_DOY.BIL",
                           output.index, #="output/budburst_tquant.csv",
                           output.budburst, #="output/budburst_1.csv",
                           year_start=2000){
  
  europe=vect(europe)
  

  ### Load 95th quantile of Tmin by month
  rast.tasmin <- rast(dir.clim)
  rast.tasmin=classify(rast.tasmin, cbind(6553500, NA))
  rast.tasmin <- crop(rast.tasmin,europe,mask=TRUE)
  
  ### Load date of budburst from 2000 to 2016
  rast.fdg=rast(dir.fdg)
  year=dim(rast.fdg)[3]
  names(rast.fdg) <- paste0("fdg_",year_start:(year_start+year-1))
  

  ### Compute mean date of budburst
  rast.fdg.mean <- mean(rast.fdg,na.rm=TRUE)
  rast.fdg.mean <- crop(terra::project(rast.fdg.mean,rast.tasmin),europe,mask=TRUE)

  
  ### Rescale each year in accurate date of budburst 
  #rast.fdg.mean[rast.fdg.mean<0] <- 365+rast.fdg.mean
  rast.fdg.mean[rast.fdg.mean<0] <- 30
  
  
  # /!\ decide what to do with negative values #
  
  # rast.fdg[[i]][rast.fdg[[i]]>180] <- 180
  # for (i in 1:year){
  #   print(i)
  #   rast.fdg[[i]][rast.fdg[[i]]<0] <- 365+rast.fdg[[i]]
  #   rast.fdg[[i]][rast.fdg[[i]]>180] <- 180
  # }
  
  rm(rast.fdg)
  
  ### create a function for each species LT50 
  year.month=c(31,28,31,30,31,30,31,31,30,31,30,31)
  year.cumul=cumsum(year.month)
  rast.cumul=rep(rast.fdg.mean,12)
  rast.cumul.bis=rep(rast.fdg.mean,12)
  rast.cumul[[1]]=ifel(rast.fdg.mean/year.cumul[1]>=1,
                       yes=year.month[1],
                       no= rast.fdg.mean)
  rast.cumul.bis[[1]]=ifel((rast.fdg.mean-30)/year.cumul[1]>=1,
                           yes=year.month[1],
                           no= (rast.fdg.mean-30))
  for (i in 2:12){
    print(i)
    rast.cumul[[i]]=ifel(rast.fdg.mean/year.cumul[i]>=1,
                         yes =year.month[i],
                         no= trunc(rast.fdg.mean/year.cumul[i-1])*(rast.fdg.mean-year.cumul[i-1]))
    rast.cumul.bis[[i]]=ifel((rast.fdg.mean-30)/year.cumul[i]>=1,
                             yes =year.month[i],
                             no= trunc((rast.fdg.mean-30)/year.cumul[i-1])*((rast.fdg.mean-30)-year.cumul[i-1]))
    gc()
  }
  
  rast.cumul=rast.cumul-rast.cumul.bis
  rast.temp=sum(rast.cumul*(rast.tasmin))/30
  
  rm(rast.cumul,rast.cumul.bis)
  
  names(rast.temp)="tmin"
  names(rast.fdg.mean)="fdg"
  
  write.csv(as.data.frame(rast.temp,xy=TRUE),output.index)
  write.csv(as.data.frame(rast.fdg.mean,xy=TRUE),output.budburst)
  return(list(rast.temp=as.data.frame(rast.temp,xy=TRUE),rast.fdg.mean=as.data.frame(rast.fdg.mean,xy=TRUE)))
}

#' Load frost index from cerra
#' 
#' @description Load files of minimum temperature computed with cdo
#' @param europe 
#' @param dir.file
#' @param output.tmin
#' @return dataframe of minimum temperature
#'

get_frostindex_cerra<-function(europe,
                               dir.file,
                               output.tmin="output/tmin_cerra.csv"){
  tmin<-rast(dir.file)
  tmin<-crop(project(tmin,
                     "epsg:4326"),
             vect(europe))
  names(tmin)="tmin"
  fwrite(as.data.frame(tmin,xy=TRUE),output.tmin)
  return(as.data.frame(tmin,xy=TRUE))

}


#' Load frost index from chelsa
#' 
#' @description Load files of minimum temperature computed with cdo
#' @param europe 
#' @param dir.file
#' @param output.tmin
#' @return dataframe of minimum temperature
#'
get_frostindex_chelsa<-function(europe,
                               dir.file="data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc",
                               output.tmin="output/tmin_chelsa.csv"){
  tmin=min(rast("data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc"),na.rm=FALSE)
  tmin=classify(tmin, cbind(6553500, NA)) #set as NA default value
  names(tmin)="tmin"
  fwrite(as.data.frame(tmin,xy=TRUE),output.tmin)
  return(as.data.frame(tmin,xy=TRUE))
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 7 - Compute safety margins ####
#' @description Combination of both sft m
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute hydraulic and frost safety margins
#' 
#' @description Compute safety margin for selected species
#' @param df.traits df of traits of each species
#' @param psimin df of psimin computed over europe
#' @param rast.temp df of frost index at leaf out date over europe
#' @param rast.fdg df of leaf out day
#' @param europe europe extent polygons
#' @return one dataframe, for each species one column for HSM, 2 types of FSM spring,
#'  1 FSM winter, 1 LT50 estimates for spring
compute.sfm <- function(df.traits,
                        rast.temp,
                        rast.fdg,
                        psimin,
                        europe){
  #load index
  rast.temp=rast(fread(rast.temp)[,-1],crs="epsg:4326")
  rast.fdg=rast(fread(rast.fdg)[,-1],crs="epsg:4326")
  frost.index.winter=min(rast("data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc"),na.rm=FALSE)
  frost.index.winter=classify(frost.index.winter, cbind(6553500, NA)) #set as NA default value
  names(frost.index.winter)="tmin.winter"
  frost.index.winter=resample(crop(frost.index.winter,
                                   vect(europe),
                                   mask=TRUE),
                              rast.temp,
                              method="near")
  rast.psimin=resample(rast(psimin[,c("x","y","psi")],crs="epsg:4326"),
                       rast.temp,
                       method="near")
  # Load rast.temp & rast.fdg and bind them
  df.sfm=as.data.frame(c(rast.temp,rast.fdg,frost.index.winter,rast.psimin),xy=TRUE)
  
  
  # COmmpute fsm
  for (i in 1:dim(df.traits)[1]){
    LT50=df.traits$LT50.mean[i]
    P50=df.traits$PX.mu[i]
    df.sfm$LT50spring=-5+30*(5+LT50)/(2*df.sfm$fdg)
    df.sfm$fsm_LT50spring=df.sfm$tmin-df.sfm$LT50spring
    df.sfm$fsm_LT50moy=df.sfm$tmin-(LT50*0.4)
    df.sfm$fsm_LT50winter=df.sfm$tmin.winter-LT50
    df.sfm$hsm_=df.sfm$psi-P50*1000
    colnames(df.sfm)[grepl("^LT50spring$",colnames(df.sfm))] <- paste0("LT50spring.",df.traits$sp.ind[i]) 
    colnames(df.sfm)[grepl("fsm_",colnames(df.sfm))] <- c(paste0("fsm.spring.",df.traits$sp.ind[i]),
                                                          paste0("fsm.moy.",df.traits$sp.ind[i]),
                                                          paste0("fsm.winter",df.traits$sp.ind[i]))
    colnames(df.sfm)[grepl("hsm_",colnames(df.sfm))] <- c(paste0("hsm.",df.traits$sp.ind[i]))
  }
  
  write.csv(df.sfm,"output/df.sfm.csv")
  return(df.sfm)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 8 - Get LT50 and P50 traits ####
#' @description Load and clean LT50 database
#' @authors Maximilian Larter (INRAE - Biogeco) & Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get species list 
#' 
#' @description get european species list from euforgen
#' @return  list of European species 
get.specieslist <- function(){
  species.list=list.files("data/chorological_maps_dataset/")
  return(species.list)
}


#'Get LT50 database
#'
#'@description From Constance database and measures performed in Clermont campaign
#'2022, we build a dataset with consistent measures
#'@return Dataframe with LT50 per European species and sd associated when available

get.LT50 <- function(){
  ### Load LT50 database from litterature ###
  ###########################################
  df.LT50.lit <- read.csv2(file = "data/Species traits/df_LT50_sp_checked.csv", header = T,na.strings = c("", "NA")) # made by code in "R/matchsp_get_tax.R"
  df.LT50.lit$New_continent <- "Continent"
  # remove NA lines (where temperature was reported as "below xx°")
  df.LT50.lit %>% filter(!is.na(Temperature_of_resistance)) -> df.LT50.lit 
  
  ### cleaning data coding ###
  ############################
  df.LT50.lit %>%
    mutate(across(.cols=c("Freezing_rate_.K.h.1.","Temperature_of_resistance","Duration_of_min_temp_exposure_.hours.","Longitude","Latitude"),
                  ~as.numeric(.))) %>% 
    mutate(Country = case_when(!is.na(Country) ~ Country,
                               # Morin et al Tree Physiol 2007
                               Reference == "Morin et al Tree Physiol 2007" & Species == "Quercus robur" ~ "Germany",
                               Reference == "Morin et al Tree Physiol 2007" & Species == "Quercus ilex" ~ "France",
                               Reference == "Morin et al Tree Physiol 2007" & Species == "Quercus pubescens" ~ "France",
                               # Sakai et al Ecology 1981 correct country below
                               Reference == "Sakai et al Ecology 1981" & Location %in% c("Canberra Botanical Garden", "Craigieburn", "Mount Ginini", "Snowy Mountains", "Tasmania") ~ "Australia",
                               Reference == "Sakai et al Ecology 1981" & Location %in% c("Christchurch Botanical Garden", "Arthurs Pass", "Jacksons", "Porters Pass", "Otira", "Waimakariri") ~ "New Zealand",
                               Reference == "Sakai et al Ecology 1981" & Location %in% c("Mount Fuji", "Hakodate", "Mount Tsukuba", "Sapporo, Hokkaido", "Shiojiri", "Tokyo", "Yamabe, Hokkaido") ~ "Japan",
                               Reference == "Sakai et al Ecology 1981" & Species %in% c("Araucaria heterophylla", "Callitris endlicheri", "Diselma archeri",  "Microcachrys tetragona") ~ "Australia", 
                               Reference == "Sakai et al Ecology 1981" & Species %in% c("Callitris oblonga", "Dacrycarpus dacrydioides", "Dacrydium cupressinum", "Halocarpus bidwillii", "Libocedrus bidwillii", "Phyllocladus alpinus", "Pinus pseudostrobus", "Podocarpus nivalis") ~ "New Zealand",
                               # Bannister New Zealand J. of Bot. 1990 al grown in NZ
                               Reference == "Bannister New Zealand J. of Bot. 1990" ~ "New Zealand",
                               # Sakai Can J. Bot. 1983 "himalaya" species were sent from Bot garden in PNG or from field in Nepal
                               Reference == "Sakai Can J. Bot. 1983" & Location == "Himalaya" & Species %in% c("Abies spectabilis", "Larix potaninii", "Tsuga dumosa", "Juniperus monosperma", "Juniperus squamata", "Picea smithiana") ~ "Nepal",
                               Reference == "Sakai Can J. Bot. 1983" & Species %in% c("Picea orientalis") ~ "Japan",
                               Reference == "Sakai Can J. Bot. 1983" & Species %in% c("Picea smithiana") ~ "Nepal",
                               Reference == "Darrow et al New Zealand J. of Bot. 2001" ~ "New Zealand",
                               Comment == "Reported from Oohata & Sakai 1982" ~ "Japan",
                               Reference == "Read & Hill Aust. J. Bot. 1988" ~ "Australia",
                               Reference == "Harrison et al Plant Physiol. 1978" ~ "USA",
                               Reference == "Alberdi et al Phytochemistry 1989" ~ "Chile",
                               Reference ==  "Bannister New Zealand J. of Bot. 2003" ~ "New Zealand",
                               Reference ==  "Goldstein et al Oecologia 1985" ~ "Venezuela",
                               TRUE ~ Country),
           Organ = case_when(!is.na(Organ) ~ Organ,
                             Comment == "Reported in Bannister New Zealand J. of Bot. 2007" ~ "Leaf",
                             Reference == "Read & Hill Aust. J. Bot. 1988" ~ "Leaf",
                             Reference == "Durham et al Physiol. Plant. 1991" ~ "Leaf"),
           Freezing_rate_.K.h.1. = case_when(!is.na(Freezing_rate_.K.h.1.) ~ Freezing_rate_.K.h.1.,
                                             Reference == "Darrow et al New Zealand J. of Bot. 2001" ~ 4,
                                             Reference ==  "Goldstein et al Oecologia 1985" ~ 10,
                                             TRUE ~ Freezing_rate_.K.h.1.),
           Method = case_when(!is.na(Method) ~ Method,
                              Reference == "Darrow et al New Zealand J. of Bot. 2001" ~ "Visual scoring",
                              Reference == "Alberdi et al Phytochemistry 1989" ~ "Visual scoring",
                              Reference ==  "Bannister New Zealand J. of Bot. 2003" ~ "Visual scoring",
                              Reference ==  "Goldstein et al Oecologia 1985" ~ "TTC",
                              TRUE ~ Method)
    ) -> df.LT50.lit
  
  ### Load LT50 database from 2022 campaign ###
  #############################################
  df.LT50.2022 <- read.csv2("data/Species traits/camp_LT50_DB.csv") |> # made by code in "Analysis new campaign data.Rmd" : output file
    bind_rows(read.csv2("data/Species traits/new_data_2023_prelim.csv"))
  
    
    #clean rm species with bad fits:
  df.LT50.2022 %>% filter(!Species %in% c("Abies balsamea", 
                                          "Acer saccharum",
                                          "Betula papyrifera",
                                          "Larix laricina",
                                          "Tsuga canadensis", 
                                          "Pinus banksiana",
                                          "Picea mariana")) -> df.LT50.2022
  df.LT50.2022$Thawing_rate_.K.h.1. <- as.character(df.LT50.2022$Thawing_rate_.K.h.1.)
  
  ### merge together ###
  ######################
  full_join(df.LT50.lit, df.LT50.2022) -> df.LT50
  df.LT50 <- df.LT50 %>% select(-Family, -Genus, -Phyllum)
  
  # country lists
  S_A <- c("Argentina", "Brazil", "Chile", "South America", "Venezuela")
  Asia <- c("China", "Iran", "Israel" , "Japan", "Korea", "Russia", "Taiwan", "Turkey", "Asia Minor", "Himalaya", "Nepal")
  E_U <- c("Austria", "Croatia", "Czech Rep.", "Denmark", "England", "Finland", "France", "Germany", "Iceland", "Italy", "Poland", "Romania", "Serbia", "Slovakia", "Spain", "Sweden", "Swiss", "Switzerland", "Ukraine", "Europe")
  N_A <- c("Canada", "Mexico", "Quebec", "USA")
  Oce <- c("Australia", "New Guinea", "Papua New Guinea", "New Zealand")
  Africa <- c("South Africa", "Canary Islands")
  
  df.LT50 %>% mutate(New_continent = case_when(Country %in% E_U  | Provenance %in% E_U ~ "E_U",
                                               Country %in% Asia | Provenance %in% Asia  ~ "Asia",
                                               Country %in% N_A | Provenance %in% N_A  ~ "N_A",
                                               Country %in% S_A | Provenance %in% S_A   ~ "S_A",
                                               Country %in% Africa | Provenance %in% Africa  ~ "Africa",
                                               Country %in% Oce | Provenance %in% Oce ~ "Oceania",
                                               is.na(Country) ~ Provenance,
                                               TRUE ~ Provenance)) -> df.LT50
  rm(S_A,Africa,Asia,E_U,N_A,Oce)
  
  bud <- c("Lateral bud",  "Bud", "Bud base vascular tissue", "Leaf primordia", "Primordial shoot", "Floral primordial", "Procambium")
  flowerbud <- c("Flower Bud", "Flower bud", "Female flower", "Male flower")
  branch <- c("Cane", "Basal stem","Wood", "Twig", "Shoot", "Stem", "Upper-crown shoot", "Twig cambium", "Xylem", "Xylem parenchyma", "Shoot vascular tissue", "Pith", "Pith parenchyma", "Phloem", "Cambial meristem", "Cambium", "Cortex", "Branch", "branch")
  leaf <- c("Needle", "Inner crown leaf", "Leaf", "Upper crown leaf")
  root <- c("Root", "Root cambium")
  df.LT50 <- df.LT50 %>% filter(Organ %in% c(bud, flowerbud, branch, leaf, root)) 
  
  # classify in organ type
  df.LT50$OrganGroup <- NA
  df.LT50$OrganGroup[df.LT50$Organ %in% bud] <- "bud"
  df.LT50$OrganGroup[df.LT50$Organ %in% flowerbud] <- "flowerbud"
  df.LT50$OrganGroup[df.LT50$Organ %in% branch] <- "branch"
  df.LT50$OrganGroup[df.LT50$Organ %in% leaf] <- "leaf"
  df.LT50$OrganGroup[df.LT50$Organ %in% root] <- "root"
  
  rm(bud,flowerbud,branch,leaf,root)
  
  names(df.LT50) <- gsub(names(df.LT50), pattern = " ", replacement = "_")
  
  ### filter for analysis ###
  ###########################
  
  growthform.select <- c("tree","Tree") #"Small tree",
  plantage.select <- c("Mature",NA)
  organgroup.select <-c("branch","bud","flowerbud","leaf") #
  method.select <- c("EL", "Visual scoring") #
  period <- c("Winter")
  
  df.LT50 %>% 
    filter(Growth_form %in% growthform.select &
             Period %in% period &
             #Plant_age %in% plantage.select &
             OrganGroup %in% organgroup.select &
             Method %in% method.select 
    ) %>%
    filter(Temperature_of_resistance>(-100)) %>% 
    filter(!(Species %in% c("Abies balsamea",
                            "Betula platyphylla",
                            "Salix sachalinensis", 
                            "Larix laricina", 
                            "Pinus banksiana", 
                            "Pinus strobus", 
                            "Thuja occidentalis",
                            "Pinus resinosa") &
               Temperature_of_resistance == -196)) %>%
    filter(!((Species=="Larix decidua")& (Reference=="Charrier et al Tree Phys 2013"))) %>% # measure that do not exist
    filter(!(Reference=="Vitra et al New Physiol. 2017" & Date!="day 29")) %>% #remove measures coming from temporal dehardening series
    filter(!(Species %in% c("Juglans nigra") &
               Organ == "Cortex")) %>% 
    filter(Type == "LT50") %>% 
    filter(New_continent=="E_U") %>% 
    mutate(data.quality=NA) -> df.LT50.select
  
  ### Test of filtering ###
  #########################
  unique(df.LT50.select$Thawing_rate_.K.h.1.)
  df.LT50.select %>% 
    filter(Plant_age == "Mature" & 
             Fit %in% c("sigmoid","Sigmoid") &
             Method == "EL" &
             OrganGroup == "branch") %>%     
    group_by(Reference,Species) %>%
    slice(which.min(Temperature_of_resistance)) %>%
    ungroup()  -> df.LT50.select.strong
  
  
  thawing.rate <- c("2","3","5","7")
  #sp.deleted <- setdiff(unique(df.LT50.select$Species),unique(df.LT50.select.strong$Species))
  df.filtered <- as.data.frame(matrix(nrow = 0,
                                      ncol=dim(df.LT50.select)[2]))
  colnames(df.filtered) <- colnames(df.LT50.select)
  for (sp in unique(df.LT50.select$Species)){
    print(sp)
    print("step 1")
    #### All constraints
    df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" & 
                              df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                              df.LT50.select$Method == "EL" & 
                              df.LT50.select$OrganGroup == "branch" &
                              df.LT50.select$Fit %in% c("sigmoid","Sigmoid")
                              ,] %>%     
      group_by(Reference,Species) %>%
      slice(which.min(Temperature_of_resistance)) %>%
      ungroup() 
    if (dim(df.sp)[1]>0){
      df.sp$data.quality=1
      df.filtered <- rbind(df.filtered,df.sp)
    }else{
      #### Release "sigmoid" constraint
      print("step 2")
      df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" & 
                                df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                                df.LT50.select$Method == "EL" & 
                                df.LT50.select$OrganGroup == "branch" 
                              ,] %>%     
        group_by(Reference,Species) %>%
        slice(which.min(Temperature_of_resistance)) %>%
        ungroup() 
      if (dim(df.sp)[1]>0){
        df.sp$data.quality=2
        df.filtered <- rbind(df.filtered,df.sp)
      }else{
        #### Release "branch" only constraint t "branch" "bud"
          print("step 3")
          df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" & 
                                    df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                                    df.LT50.select$Method == "EL" & 
                                    df.LT50.select$OrganGroup %in% c("branch","bud")
                                    ,] %>%     
          group_by(Reference,Species) %>%
          slice(which.min(Temperature_of_resistance)) %>%
          ungroup() 
        if (dim(df.sp)[1]>0){
          df.sp$data.quality=3
          df.filtered <- rbind(df.filtered,df.sp)
        }else{
          #### Release "method" constraint
              print("step 4")
              df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" &
                                          df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                                          df.LT50.select$OrganGroup %in% c("branch","bud")
                                    ,]%>%     
              group_by(Reference,Species) %>%
              slice(which.min(Temperature_of_resistance)) %>%
              ungroup() 
            if (dim(df.sp)[1]>0){
              df.sp$data.quality=4
              df.filtered <- rbind(df.filtered,df.sp)
            }else{
              #### Release "thawing rate" constraint
              print("step 5")
              df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" &
                                        df.LT50.select$OrganGroup %in% c("branch","bud")
                                      ,]%>%     
                group_by(Reference,Species) %>%
                slice(which.min(Temperature_of_resistance)) %>%
                ungroup() 
              if (dim(df.sp)[1]>0){
                df.sp$data.quality=5
                df.filtered <- rbind(df.filtered,df.sp)
              }else{
              #### Release "age" constraint
                print("step 6")
                df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$OrganGroup %in% c("branch","bud")
                             ,] %>%     
                group_by(Reference,Species) %>%
                slice(which.min(Temperature_of_resistance)) %>%
                ungroup() 
              if (dim(df.sp)[1]>0){
                df.sp$data.quality=6
                df.filtered <- rbind(df.filtered,df.sp) 
              }else{
                print("step 7")
                df.filtered <- rbind(df.filtered,
                                     df.LT50.select[df.LT50.select$Species==sp
                                                    ,]%>%     
                                       group_by(Reference,Species) %>%
                                       slice(which.min(Temperature_of_resistance)) %>%
                                       ungroup() %>% 
                                       mutate(data.quality=7))
              }
            }
          }
      }
    }
    }
    rm(df.sp)}


  ### Summarise per species ###
  #############################
  df.species.traits <- df.filtered %>% 
    filter(New_continent=="E_U") %>%
    group_by(Species,Period,data.quality) %>% 
    summarise(LT50.mean=mean(Temperature_of_resistance,na.rm=TRUE),
              LT50.sd=sd(Temperature_of_resistance,na.rm=TRUE)) %>% 
    ungroup() %>% 
    arrange(LT50.mean) %>% 
    mutate(rank_filt=row_number())
  
  ## For raw data
  df.species.traits.raw <- df.LT50.select %>% 
    filter(New_continent=="E_U") %>%
    group_by(Species,Period) %>% 
    summarise(LT50.mean=mean(Temperature_of_resistance,na.rm=TRUE),
              LT50.sd=sd(Temperature_of_resistance,na.rm=TRUE)) %>% 
    ungroup() %>% 
    arrange(LT50.mean) %>% 
    mutate(rank=row_number()) %>% 
    left_join(df.species.traits[,c("Species","rank_filt")],by="Species")
  

  
    
 ### Few plots ###
  #################
  # ## Plot difference between mean computed over raw dataset and with filtered dataset
  # df.species.traits %>%
  #   filter(!is.na(LT50.sd)) %>%
  #   rename(`LT50.mean.filter`="LT50.mean",
  #          `LT50.sd.filter`="LT50.sd") %>%
  #   left_join(df.species.traits.raw %>%
  #               filter(!is.na(LT50.sd)) %>%
  #               rename(`LT50.mean.raw`="LT50.mean",
  #                      `LT50.sd.raw`="LT50.sd"),
  #             by="Species") %>%
  #   ggplot(aes(LT50.mean.filter-LT50.mean.raw,Species))+
  #   geom_point()
  # 
  # ## Boxplots for species with several measures
  # list=df.filtered  %>%
  #   mutate(data.quality=as.factor(data.quality)) %>% 
  #   filter(New_continent=="E_U") %>%
  #   # group_by(Species) %>%
  #   # filter(n()>1) %>%
  #   # ungroup() %>%
  #   ggplot(aes(Temperature_of_resistance,Species))+
  #   geom_line()+
  #   geom_point(aes(shape=data.quality))+
  #   theme_bw()+
  #   theme(axis.title=element_blank() )
  # ggsave(list,
  #        filename = "species_list.pdf",
  #        path="output/",
  #        device="pdf",
  #        height=12)
  # 
  # ## Plots with changes in species ranking before/after LT50 filtering
  # df.species.traits.raw %>% 
  #   filter(abs(rank-rank_filt)>5) %>% 
  #   ggplot(aes(rank-rank_filt,Species))+
  #   geom_point()
  # 
  # 
  # df.species.traits.raw %>%
  #   pivot_longer(cols = c("rank","rank_filt"),names_to = "rank.type",values_to = "rank") %>%
  #   mutate(rank.type=as.factor(rank.type)) %>%
  #   ggplot(aes(x = rank.type, y = rank, group = Species)) +
  #   geom_line(aes(color = Species, alpha = 1), size = 2) +
  #   geom_point(aes(color = Species, alpha = 1), size = 2) +
  #   geom_point(color = "#FFFFFF", size = 1) +
  #   scale_y_reverse(breaks = 1:nrow(df.species.traits.raw)) +
  #   scale_x_discrete(breaks = 1:10) +
  #   theme_bw()+
  #   theme(legend.position = 'none',
  #         axis.title = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.ticks.y= element_blank()) +
  #   geom_text(data = df.species.traits.raw %>% arrange(rank),
  #             aes(label = Species, x = .95) , hjust = .5,
  #             color = "#888888", size = 3) +
  #   labs(x = '', y = 'Rank', title = 'Changes in LT50 species ranking after filtering')
  # 
  # 
  
  
  write.csv(df.species.traits,"output/df_LT50_filtered.csv",row.names = FALSE)

  return(list(df.LT50.raw=df.LT50.select,
              df.LT50.cor=df.filtered,
              df.LT50sp.raw=df.species.traits.raw,
              df.LT50sp.cor=df.species.traits))
  
}




#' Get P50 database 
#' 
#' @description Load P50 from Max database
#' @return  Dataframe with P50 per European species 
get.P50 <- function(){
  species.list=get.specieslist()
  # df.P50.lit <- read.csv2("data/Species traits/BdD_full_by_species.csv")
  # df.P50.select <- df.P50.lit %>% 
  #   select(colnames(df.P50.lit)[!grepl("LT",colnames(df.P50.lit))]) %>% 
  #   select(-Dmean,-Psi_min_MD,-WD,-Plot_continent,-MHZ,-nplot,-ntree,-D99,-score,-X.gs90,-X.tlp,-densityMEANcor,-Exo,-genus,-GBIF_species,-class) %>% 
  #   filter(New_continent=="E_U") %>% 
  #   filter(!is.nan(P50tot))
  df.p50.msp<-read.csv2("data/Species traits/p50_nmsp.csv") |> 
    filter(species.binomial %in% species.list) |> 
    mutate(group=tolower(group),
           p88.mu=p50-50/slope,
           p50.mu=p50,
           bdd="martinsp") |> 
    select(group,species.binomial,p50.mu,p88.mu,bdd) |> 
    filter(!(group=="angiosperm"&is.na(p88.mu)))
  df.p50.lit <- read.csv2("data/Species traits/2022_10_20_cleaned_xft.csv") %>% 
    select(uniqueID,lat,long,XFT.database,Group,Family,Genus,Species,
           Developmental.stage,Growth.form,P50,P50.SD,P12,P88,Curve,Equation,
           psi.min.predawn..MPa.,psi..min.midday..MPa.,P50.number.of.samples) %>% 
    mutate(species.binomial=paste(Genus,Species),
           bdd="hammond") %>% 
    filter((species.binomial %in% species.list)&
             (!species.binomial %in% df.p50.msp$species.binomial )) |> 
    filter(Curve=="S") %>%
    filter(Growth.form=="T") %>% 
    filter(Developmental.stage =="A") |> 
    group_by(Group,Genus,Species,species.binomial,bdd) %>% 
    summarise(P50.mu=mean(as.numeric(P50),na.rm=TRUE),
              P50.sd=sd(as.numeric(P50),na.rm=TRUE),
              P88.mu=mean(as.numeric(P88),na.rm=TRUE),
              P88.sd=sd(as.numeric(P88),na.rm=TRUE)) %>% 
    ungroup() |> 
    rename_with(.cols=everything(),
                tolower) |> 
    select(group,species.binomial,p50.mu,p50.sd,p88.mu,p88.sd,bdd)
  df.p50=bind_rows(df.p50.lit,df.p50.msp)
  write.csv(df.p50,"output/df_P50_filtered.csv",row.names = FALSE)
  return(df.p50)
}

#' Get LT50/P50 database 
#' 
#' @description Load LT50/P50 for european species
#' @param dir.data data directory
#' @param dir.file directory of species file
#' @return dataframe of LT50/P50 traits per species
#' 

get.traits <- function(dir.distribution="data/chorological_maps_dataset",
                       df.P50, #=df.P50
                       df.LT50){ #=df.LT50$df.LT50sp.cor
  # df.P50 <- read.csv("output/df_P50_filtered.csv") 
  # df.LT50 <- read.csv("output/df_LT50_filtered.csv")
  
  df.traits.raw <- df.LT50 %>% 
    left_join(df.P50,by=c("Species"="species.binomial")) |> 
    rename_with(.cols=everything(),
                tolower)
  # df.traits.raw <- full_join(df.P50,df.LT50,by="Species")
  # df.traits.missing <- df.traits.raw %>% 
  #   filter(is.na(P50tot) | is.na(LT50.mean))
  df.traits <- df.traits.raw %>% 
    filter(!is.na(p50.mu)) %>% 
    filter(!is.na(lt50.mean)) %>% 
    select(species,group,p50.mu,p50.sd,p88.mu,p88.sd,lt50.mean,lt50.sd,data.quality) %>% #MATmean,MAPmean,Leaf_phenology,
    unique() %>% 
    rename(`species.binomial`="species") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2)),
           species.name=str_replace(species.binomial," ","_")) %>% 
    relocate(c("species.name","sp.ind"),.after="species") %>% 
    mutate(p.trait=case_when(group=="gymnosperm"~"p50",
                             group=="angiosperm"&!is.na(p88.mu)~"p88",
                             group=="angiosperm"&is.na(p88.mu)~"p50"),
           px.mu=case_when(p.trait=="p88"~p88.mu,
                           p.trait=="p50"~p50.mu),
           px.sd=case_when(p.trait=="p88"~p88.sd,
                           p.trait=="p50"~p50.sd)) %>% 
    filter(!(species.binomial=="Pinus contorta"&group=="angiosperm"))
  
  for (i in 1:dim(df.traits)[1]){
    species.files=list.files(file.path(dir.distribution,df.traits$species.binomial[i],"shapefiles"))
    if(length(species.files[grepl(paste0(df.traits$species.name[i],"_plg_clip"),species.files)])>1){
      df.traits$file[i]=paste0(df.traits$species.name[i],"_plg_clip")
    } else if(length(species.files[grepl(paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg_clip"),species.files)])>1){
      df.traits$file[i]=paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg_clip")
    } else if(length(species.files[grepl(paste0(df.traits$species.name[i],"_plg"),species.files)])>1) {
      df.traits$file[i]=paste0(df.traits$species.name[i],"_plg")
    } else if(length(species.files[grepl(paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg"),species.files)])>1) {
      df.traits$file[i]=paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg")
    } else {
      df.traits$file[i]=NA
    }
  }  
  write.csv(df.traits,"output/df_trait_filtered.csv",row.names = FALSE)
  return(df.traits)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Annexe - Useful chunks that were deleted  ####
#' @description few chunks that could be re-used
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# project raster
# SWC <- plyr::join_all(lapply(c("SWC_w_min","SWC_w_max","SWC_m_min","SWC_m_max"),
#                              function(x){
#                                assign(paste0(x,"_df"),
#                                       as.data.frame(project(get(x),texture),
#                                                     xy=TRUE))
#                              }
# ),by=c("x","y"))
# colnames(SWC) <- c("x","y","SWC_w_min","SWC_w_max","SWC_m_min","SWC_m_max")
# 
# 


#load rasters of variables of interest
# for (x in (1:length(Abv))){
#   assign(Abv[x],rast(list.files(file.path(dir.data,Variables[x]),full.names = TRUE)[1]))
# }


# psi_min <- as.data.frame(c(rast(texture,crs="epsg:4326"),SWC),xy=TRUE) %>% 
#   filter(clay!=0&sand!=0) %>% 
#   #categorize texture
#   mutate(texture=case_when(clay>60~1,
#                            clay>35&clay<=60~2,
#                            clay<=35&sand<15~3,
#                            clay>18&clay<=35&sand>15~4,
#                            clay<=18&sand>15&sand<65~4,
#                            clay<=18&sand>=65~5)) %>% 
#   left_join(texture_pot,by=c("texture")) %>% 
#   #left_join(SWC,by=c("x","y")) %>% 
#   mutate(Psi_w_min=Psi_e*(SWC_w_max/SWC_w_min)^(b),
#          SWC_w_ratio=SWC_w_max/SWC_w_min,
#          # Psi_m_min=Psi_e*(SWC_m_max/SWC_m_min)^(b),
#          # SWC_m_ratio=SWC_m_max/SWC_m_min,
#          Psi_w_med=texture_pot$Psi_e[3]*(SWC_w_max/SWC_w_min)^(texture_pot$b[3]))
# # Psi_m_med=texture_pot$Psi_e[3]*(SWC_m_max/SWC_m_min)^(texture_pot$b[3]))

# texture_pot <- data.frame(texture=c(5,4,3,2,1),
#                           Psi_e=c(-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
#                           b=c(2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933),
#                           teta_s=c(0.3845,0.4155,0.421,0.5005,0.576))




## gathering
# psi_min=cbind(psi_top %>% 
#                 rename(clay_t="clay",
#                        sand_t="sand",
#                        oc_t="oc",
#                        texture_t="texture"),
#               psi_sub[,c(3,4,5,8,17,18,19)] %>% 
#                 rename(clay_s="clay",
#                        sand_s="sand",
#                        oc_s="oc",
#                        texture_s="texture"))%>% 
#   filter(clay_t!=0&sand_t!=0&oc_t!=0) %>% 
#   relocate("texture_t",.after="oc_t") %>% 
#   relocate(c("clay_s","sand_s","oc_s","texture_s"),.after="texture_t") %>% 
#   relocate("SWC_ratio",.after="SWC_w_max")  %>% 
#   mutate(psi_t_VG=case_when(SWC_w_min>teta_s~-0.1,TRUE~psi_t_VG),
#          psi_s_VG=case_when(SWC_w_min>teta_s~-0.1,TRUE~psi_s_VG)) %>% 
#   rowwise() %>% 
#   mutate(psi_min_BG=max(psi_t_BG,psi_s_BG,na.rm=TRUE),
#          psi_min_CB=max(psi_t_CB,psi_s_CB,na.rm=TRUE),
#          psi_min_VG=max(psi_t_VG,psi_s_VG,na.rm=TRUE)) %>% 
#   ungroup()



# plot_psi %>% 
#   mutate(across(.cols=matches(c("psi")),
#                 ~ cut(.,
#                       breaks=c(-Inf, -20000,  -5000, -1000 , Inf),
#                       labels=c("<-20MPa", "-20<Psi<-5MPa","-5<Psi<-1MPa",
#                                "-1/0"))))  %>% 
#   select(matches(c("^x$","^y$","psi"))) %>% 
#   relocate(y,.after=x) %>% 
#   pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
#   ggplot()+
#   geom_tile(aes(x=x,y=y,fill=value)) +
#   theme_bw() +
#   theme(axis.title=element_blank(),
#         legend.key.size = unit(0.5,"cm"))+
#   labs(fill="Potentiel (MPa)")+
#   facet_wrap(~name)+
#   scale_fill_brewer(palette="RdYlBu")+
#   coord_quickmap()

#get_chelsa_wai
# chelsa.files=list.files(dir.chelsa)
# pet.files=chelsa.files[grepl("pet_penman_",chelsa.files)]
# prec.files=chelsa.files[grepl("pr_",chelsa.files)]
# test.files=chelsa.files[grepl("gdd",chelsa.files)]
# # Loop on all years
# for (i in substr(c(106:109), 2, 3)) {
#   df.loc=cbind(df.loc,
#                # extract(rast(file.path(dir.chelsa,
#                #                        chelsa.files[grepl(paste0("pet_penman_",i),
#                #                                           chelsa.files)])),
#                #         df.loc),
#                extract(rast(file.path(dir.chelsa,
#                                       chelsa.files[grepl(paste0("pr_",i),
#                                                          chelsa.files)])),
#                        df.loc[,c("x","y")])[,2])
#   # a modif avec pet
#   colnames(df.loc)=c(colnames(df.loc)[1:(dim(df.loc)[2]-1)],paste0("pr_",i))
# }
# db.wai=cbind(df.loc[,1:2],
#              apply(df.loc[,3:dim(df.loc)[2]],MARGIN = 1,mean))

# list=mget(ls()[grepl("psi_h",ls())])


#create a df for a year, with values of LT50 depending of budburst date

# Create LT50 time series
# LT50.dt=c(LT50+seq(0,day_budburst)*((-5-LT50)/day_budburst),
#           rep(-5,day_hardenning-1-day_budburst-1),
#           -5+seq(0,yday(as.Date(paste0(year.select,"-12-31")))-day_hardenning)*((LT50+5)/(yday(as.Date(paste0(year.select,"-12-31")))-day_hardenning)))

# Old filtering of LT50 data
# for (sp in sp.deleted){
#   print(sp)
#   df.sp <- df.LT50.select %>% 
#     filter(Species==sp)
#   
#   if(dim(df.sp)[1]>1&dim(df.sp)[1]<5){
#     tryCatch(
#       {
#         # Look how wide is the range
#         if(max(df.sp$Temperature_of_resistance,na.rm=TRUE)-min(df.sp$Temperature_of_resistance,na.rm=TRUE)>10){
#           # Check for plant age first
#           med.mat=median(df.sp[df.sp$Plant_age %in% c("Mature",NA),"Temperature_of_resistance"])
#           med.other=median(df.sp[! df.sp$Plant_age %in% c("Mature",NA),"Temperature_of_resistance"])
#           tryCatch({
#             if(abs(med.mat-med.other)>10){
#               df.sp <- df.sp[df.sp$Plant_age %in% c("Mature",NA),]
#               print("Some measures discarded")
#             }
#           },
#           error=function(e){print("One group not present")}
#           )
#           
#           # Check for fit quality 
#           med.sigmoid=median(df.sp[df.sp$Fit=="Sigmoid","Temperature_of_resistance"])
#           med.other=median(df.sp[df.sp$Fit!="Sigmoid","Temperature_of_resistance"])
#           tryCatch({
#             if(abs(med.sigmoid-med.other)>10){
#               df.sp <- df.sp[df.sp$Fit=="Sigmoid",]
#               print("Some measures discarded")
#             }
#           },
#           error=function(e){print("One group not present")}
#           )
#           
#           # Check for method 
#           med.el=median(df.sp[df.sp$Method=="EL","Temperature_of_resistance"])
#           med.vs=median(df.sp[df.sp$Method=="Visual scoring","Temperature_of_resistance"])
#           tryCatch(
#             {if(abs(med.el-med.vs)>10){
#               df.sp <- df.sp[df.sp$Method=="EL",]
#               print("Some measures discarded")
#             }},
#             error=function(e){print("One group not present")}
#           )
#           
#           # Check for organ group 
#           # performing t.test is not an option because groups are too small and therefore exclusion condition is not strong enough
#           med.br=median(df.sp[df.sp$OrganGroup=="branch","Temperature_of_resistance"])
#           med.other=median(df.sp[!df.sp$OrganGroup=="branch","Temperature_of_resistance"])
#           tryCatch(
#             {if(abs(med.br-med.other)>10){
#               df.sp <- df.sp[df.sp$OrganGroup=="branch",]
#               print("Some measures discarded")
#             }},
#             error=function(e){print("One group not present")}
#           )
#         }
#       },
#       error=function(e){print("Error for this species")}
#     )
#   } else if(dim(df.sp)[1]>5){
#     tryCatch(
#       {
#         # Look how wide is the range
#         if(max(df.sp$Temperature_of_resistance,na.rm=TRUE)-min(df.sp$Temperature_of_resistance,na.rm=TRUE)>10){
#           # Check for plant age first
#           tryCatch({
#             wilcox.age=wilcox.test(df.sp[df.sp$Plant_age %in% c("Mature",NA),"Temperature_of_resistance"],
#                                    df.sp[! df.sp$Plant_age %in% c("Mature",NA),"Temperature_of_resistance"])
#             if(wilcox.age$p.value<0.05){
#               df.sp <- df.sp[df.sp$Plant_age %in% c("Mature",NA),]
#               print("Some measures discarded")
#             }
#           },
#           error=function(e){print("Test not performed")}
#           )
#           
#           # Check for fit quality 
#           tryCatch({
#             wilcox.fit=wilcox.test(df.sp[df.sp$Fit=="Sigmoid","Temperature_of_resistance"],
#                                    df.sp[! df.sp$Plant_age %in% c("Mature",NA),"Temperature_of_resistance"])
#             if(wilcox.fit$p.value<0.05){
#               df.sp <- df.sp[df.sp$Fit=="Sigmoid",]
#               print("Some measures discarded")
#             }
#           },
#           error=function(e){print("One group not present")}
#           )
#           
#           # Check for method 
#           tryCatch({
#             wilcox.method=wilcox.test(df.sp[df.sp$Method=="EL","Temperature_of_resistance"],
#                                       df.sp[df.sp$Method=="Visual scoring","Temperature_of_resistance"])
#             if(wilcox.method$p.value<0.05){
#               df.sp <- df.sp[df.sp$Method=="EL",]
#               print("Some measures discarded")
#             }},
#             error=function(e){print("Test not performed")}
#           )
#           
#           # Check for organ group 
#           # performing t.test is not an option because groups are too small and therefore exclusion condition is not strong enough
#           tryCatch({
#             wilcox.organ=wilcox.test(df.sp[df.sp$OrganGroup=="branch","Temperature_of_resistance"],
#                                      df.sp[!df.sp$OrganGroup=="branch","Temperature_of_resistance"])
#             if(wilcox.organ$p.value<0.05){
#               df.sp <- df.sp[df.sp$OrganGroup=="branch",]
#               print("Some measures discarded")
#             }},
#             error=function(e){print("Test not performed")}
#           )
#         }
#       },
#       error=function(e){print("Error for this species")}
#     )
#   } 
#   df.filtered=rbind(df.filtered,df.sp)
# }