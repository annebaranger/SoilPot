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
  europe <-  raster::crop(europe, extent(-10, 46, 32, 72))
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


#' Data downloading from ERA5-land 
#' 
#' @description function used to download ERA5-land data using KrigR packages,
#'  and to save it in a folder
#' @note this function is obsolete as it takes a lot of time to download data 
#' from KrigR and this is much faster to do it directly from the website
#' @note Prior using KrigR function, it is necessary to sign in to 
#' "https://cds.climate.copernicus.eu" 
#' @param API_User API_user of user account
#' @param API_Key Key of user account
#' @param dir.data directory of project data
#' @param Dir.Shapes directory of countries masks
#' @param extent an object of SpatialExtent type that delimits the spatial zone
#' to extract. Default = NULL because countries name can be used to extract data
#' @param countries list of countries to extract ERA5-land data of
#' @param date_start "YYYY-MM-DD" beginning of the period to extract
#' @param date_end "YYYY-MM-DD" end of the period to extract
#' @param time_step "day", "month" or "year"
#' @output create folders with aggregated variables
era_data <- function(API_User="124078",
                     API_Key="ab883de0-39cd-4452-a1f0-7cd45ab252b0",
                     dir.data,
                     Dir.Shapes,
                     extent=NULL,
                     countries,
                     variables,
                     date_start,
                     date_end,
                     time_step){
  StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE)
  
  if(is.null(extent)){
    position <- unlist(lapply(countries,function(x) which(StateMask$geonunit==x)))
    extent <- StateMask[position,]
    # extent= extent(State_Shp)
  }
  lapply(variables,function(x){Dir.Var=file.path(dir.data,x)
                               dir.create(Dir.Var)
                               download_ERA(
                                Variable = x,
                                DataSet = "era5-land",
                                DateStart = date_start,
                                DateStop = date_end,
                                TResolution = time_step,
                                TStep = 1,
                                Extent = extent,
                                Dir = Dir.Var,
                                API_User = API_User,
                                API_Key = API_Key)})
}


#' Load data from TerraClimate
#'  
#' @description function used to load TerraClimate data into target envt
#' @note terraclimate data can be downloaded directly from website 
#' https://www.climatologylab.org
#' @note spatial extent is to be specified directly on the website when DLing
#' @return a dataframe of TerraClimate data
terraclimate_data <- function() {
  terraclimate <- as.data.frame(rast("data/agg_terraclimate_soil_1958_CurrentYear_GLOBE.nc"),
                             xy=TRUE)
  return(terraclimate)
}

#' Data downloading from DL soilgrid 
#' @note see Natheo script as function proposed by Soilgrid website is obsolete


#' Load data from Chelsa
#'  
#' @description function used to load Chelsa data into target envt. Bio6 variable
#' is loaded, and correspond to mean daily minimum air temperature of the coldest
#'  month during the period 1980-2010
#' @note Chelsa can be downloaded on
#'  https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
#' @note spatial extent is to be specified directly on the website when DLing
#' @return a dataframe of TerraClimate data
chelsabio6_data <- function(dir.data,dir.file,europe) {
  chelsabio6 <- crop(rast(file.path(dir.data,dir.file)),vect(europe),mask=TRUE)
  names(chelsabio6)="t.min"
  return(as.data.frame(chelsabio6,xy=TRUE))
}


#' Load FNI data from Franvr
#'  
#' @description function to load data from France National Inventory
#' @note Code from Matthieu Combaud
#' @param input.directory
#' @param zip.file
#' 
get.fni <- function(input.directory="data/IFN",
                    zip.file="export_dataifn_2005_2020.zip"){
  db.tree<-as.data.table(read.table(unz(paste0(input.directory,"/",zip.file), "ARBRE.csv"),
                                    header = TRUE, 
                                    sep = ";", 
                                    dec = ".", 
                                    encoding = "UTF-8"))
  db.cover<-as.data.table(read.table(unz(paste0(input.directory,"/",zip.file), "COUVERT.csv"), 
                                     header = TRUE,
                                     sep = ";",
                                     dec = ".",
                                     encoding = "UTF-8"))
  db.stand<-as.data.table(read.table(unz(paste0(input.directory,"/",zip.file), "PLACETTE.csv"), 
                                     header = TRUE,
                                     sep = ";",
                                     dec = ".", 
                                     encoding = "UTF-8"))
  return(list(tree=db.tree,
              cover=db.cover,
              stand=db.stand))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Computing essential variable over spatial area ####
#' @description Functions used to compute required variables to apply Campbell 
#' equations. Texture from ESDAC and ERA5. SWC daily or monthly timepath.
#' Parameters from Toth.
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Compute texture from ESDAC over a specified spatial extent
#' 
#' @description Function that load ESDAC texture data from directory
#' @note Texture data can be downloaded on ESDAC website. They encompass texture
#' from topsoil and subsoil at 1x1km resolution. Parameters are from Wosten 1999
#' @param dir.data directory of data of the project
#' @param dir.soil directory of soil data downloaded from ESDAC
#' @param europe SpatialPolygonsDataFrame of the spatial extent 
#' @return spatial dataframe (SpatRaster converted to dataframe) of clay and 
#' sand content average among topsoil and subsoil, cropped to spatial extent
get_textureESDAC <- function(dir.data,dir.soil,europe){
  # get spatial extent
  rast_model <- terra::vect(europe)
  
  # Clay 
  clay <- c(rast(file.path(dir.data,dir.soil,"STU_EU_T_CLAY.rst")),rast(file.path(dir.data,dir.soil,"STU_EU_S_CLAY.rst")))
  crs(clay) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  names(clay) <- c("topsoil","subsoil")
  clay <- terra::mask(project(clay,"epsg:4326"),rast_model)

    # Sand 
  sand <- c(rast(file.path(dir.data,dir.soil,"STU_EU_T_SAND.rst")),rast(file.path(dir.data,dir.soil,"STU_EU_S_SAND.rst")))
  crs(sand) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  names(sand) <- c("topsoil","subsoil")
  sand <- terra::mask(project(sand,"epsg:4326"),rast_model)

  # Organic content
  oc <- c(rast(file.path(dir.data,dir.soil,"STU_EU_T_OC.rst")),rast(file.path(dir.data,dir.soil,"STU_EU_S_OC.rst")))
  crs(oc) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  names(oc) <- c("topsoil","subsoil")
  oc <- terra::mask(project(oc,"epsg:4326"),rast_model)
  
  # Texture
  texture_raw <- c(clay,sand,oc)
  names(texture_raw) <- paste0(names(texture_raw),".",c("clay","clay","sand","sand","oc","oc"))
  texture_raw <- as.data.frame(texture_raw,xy=TRUE)
  
  # Pars from Wosten
  pars=list(topsoil=data.frame(texture=c(6,5,4,3,2,1),
                        psi_e=c(-0.790569415,-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
                        b=c(2.6411388301,2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933),
                        teta_r=c(0.010,0.0250,0.0100,0.0100,0.0100,0.0100),
                        teta_s=c(0.766,0.4030,0.4390,0.4300,0.5200,0.6140),
                        alpha=c(0.0130,0.0383,0.0314,0.0083,0.0367,0.0265),
                        n=c(1.2039,1.3774,1.1804,1.2539,1.1012,1.1033),
                        m=c(0.1694,0.2740,0.1528,0.2025,0.0919,0.0936)),
         subsoil=data.frame(texture=c(6,5,4,3,2,1),
                        psi_e=c(-0.790569415,-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
                        b=c(2.6411388301,2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933),
                        teta_r=c(0.010,0.0250,0.0100,0.0100,0.0100,0.0100),
                        teta_s=c(0.766,0.3660,0.3920,0.4120,0.4810,0.5380),
                        alpha=c(0.0130,0.0430,0.0249,0.0082,0.0198,0.0168),
                        n=c(1.2039,1.5206,1.1689,1.2179,1.0861,1.0730),
                        m=c(0.1694,0.3424,0.1445,0.1789,0.0793,0.0680)))
  
  texture <- sapply(names(pars),FUN=function(x){
    oc.tx <- paste0(x,".oc")
    clay.tx <- paste0(x,".clay")
    sand.tx <- paste0(x,".sand")
    assign(paste0("texture_",x),
           cbind(texture_raw[,c("x","y")],texture_raw[,grepl(x,colnames(texture_raw))]) %>%
             mutate(texture=case_when(.data[[clay.tx]]==0&.data[[sand.tx]]==0~0,
                                        .data[[oc.tx]]>15~6,
                                       .data[[clay.tx]]>60~1,
                                       .data[[clay.tx]]>35&.data[[clay.tx]]<=60~2,
                                       .data[[clay.tx]]<=35&.data[[sand.tx]]<15~3,
                                       .data[[clay.tx]]>18&.data[[clay.tx]]<=35&.data[[sand.tx]]>15~4,
                                       .data[[clay.tx]]<=18&.data[[sand.tx]]>15&.data[[sand.tx]]<65~4,
                                       .data[[clay.tx]]<=18&.data[[sand.tx]]>=65~5)) %>% 
             left_join(pars[[x]],by=c("texture")))
  },simplify=FALSE,USE.NAMES=TRUE)
  return(c(texture_raw,texture))
}

#' Compute texture from ERA5 over a specified spatial extent
#' 
#' @description Function that load ERA5 data from directory
#' @note Texture data can be downloaded on ERA5 website. They designate texture 
#' of subsoil at 10x10km resolution. Parameters can be chosen between Toth and Wosten.
#' @param dir.data directory of data of the project
#' @param dir.soil directory of soil data downloaded from ESDAC
#' @param europe SpatialPolygonsDataFrame of the spatial extent 
#' @return spatial dataframe (SpatRaster converted to dataframe) of texture class 
#' and associated parameters, cropped to spatial extent
#' 
# subsoil=data.frame(texture=c(6,5,4,3,2,1),
#                    psi_e=c(-0.790569415,-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
#                    b=c(2.6411388301,2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933),
#                    teta_r=c(0.010,0.0250,0.0100,0.0100,0.0100,0.0100),
#                    teta_s=c(0.766,0.3660,0.3920,0.4120,0.4810,0.5380),
#                    alpha=c(0.0130,0.0430,0.0249,0.0082,0.0198,0.0168),
#                    n=c(1.2039,1.5206,1.1689,1.2179,1.0861,1.0730),
#                    m=c(0.1694,0.3424,0.1445,0.1789,0.0793,0.0680))
get_textureERA5 <- function(dir.data,
                         dir.file,
                         europe,
                         pars){
  texture_ERA5=rast(file.path(dir.data,dir.file))
  texture_ERA5=crop(texture_ERA5,vect(st_as_sf(europe)),mask=TRUE)
  texture_ERA5=as.data.frame(texture_ERA5,xy=TRUE)  %>% 
    filter(slt!=0) %>% 
    mutate(slt=round(slt),
           slt=case_when(slt==1~"coarse",
                         slt==2~"medium",
                         slt==3~"mediumf",
                         slt==4~"fine",
                         slt==5~"vfine",
                         slt==6~"oc",
                         slt==7~"oc"),
           texture=case_when(slt=="coarse"~5,
                             slt=="medium"~4,
                             slt=="mediumf"~3,
                             slt=="fine"~2,
                             slt=="vfine"~1,
                             slt=="oc"~6)) %>% 
    left_join(pars,by=c("texture"))
  return(texture_ERA5[,c(1,2,4:11)])
  }

#' Soil volumetric content pre-treatment
#' 
#' @description crop era5land SWC to the good extent
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param dir.data directory of data of the project
#' @param texture spatial dataframe of soil textures in the area considered
#' @param europe SpatialPolygonsDataFrame of the spatial extent 
#' @return spatial dataframe (spatraster converted to dataframe) ofcropped soil 
#' water content layers. Each layer is a date. 
get_SWC <- function(dir.data,
                    dir.file,
                    europe){
  SWCtot <-rast(file.path(dir.data,dir.file))
  europe <- vect(st_as_sf(europe))
  SWCtot <- crop(SWCtot,europe,mask=TRUE)
  return(as.data.frame(SWCtot,xy=TRUE))
}

#' Compute weighted soil volumetric water content
#' 
#' @description Using ERA5-land SWC and a given depth, the function weight the
#' different horizons according to their thickness
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param dir.data directory of data of the project
#' @param SWCtot dataframe of SWCtot cropped to the accurate extent
#' @param depth numeric indicating to which depth swc is to be considered
#' @return dataframe with weighted swc over different horizons
#'  
weight_swc <- function(dir.data,
                       dir.file,
                       SWCtot,
                       depth){
  SWCtot=rast(SWCtot,crs="epsg:4326")
  
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
  SWC_t <- (x1*SWCtot[[grepl("swvl1", names(SWCtot))]]+
          x2*SWCtot[[grepl("swvl2", names(SWCtot))]]+
          x3*SWCtot[[grepl("swvl3", names(SWCtot))]]+
          x4*SWCtot[[grepl("swvl4", names(SWCtot))]])/(x1+x2+x3+x4)
  return(as.data.frame(SWC_t,xy=TRUE))
}

#' Compute min and max soil volumetric water content over a timeserie
#' 
#' @description Using ERA5-land SWC and a given depth, the function weight the
#' different horizons according to their thickness
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param SWC_t dataframe of timeserie of swc
#' @param timepath month or day
#' @param date_begin "YYYY-MM-DD" one unit timepath before the beginning of the
#' time serie
#' @return dataframe with min/max
#' date_begin='1949-12-01'
extr_swc <- function(SWC_t,timepath,date_begin) {
  SWC=setDT(as.data.frame(SWC_t,xy=TRUE)) #use data.table to deal with big db
  if (timepath=="month") {
    time=as.character(as.Date(date_begin) %m+% months(as.numeric(sub(".*_","",
                                                                       colnames(SWC)[c(-1,-2)]))))
  }
  if (timepath=="day"){
    time=as.character(as.Date(as.numeric(sub(".*_","",
                                             colnames(SWC)[c(-1,-2)])),
                              origin=paste0(date_begin)))
  }
  colnames(SWC)=c("x","y",time)
  # selection of max/min per year
  SWC=SWC[,loc_ID:=.I]
  SWC=melt.data.table(SWC,id.vars=c("loc_ID","x","y"))
  SWC[,year:=gsub("[-].*", "\\1",variable)]
  SWC[,month:=gsub(".*[-]([^.]+)[-].*", "\\1",variable)]
  SWC[,day:=gsub(".*[-]", "\\1",variable)]
  SWC=SWC[,.(Min=min(value),Max=max(value)),by="loc_ID,x,y,year"]
  SWC=SWC %>% 
    group_by(loc_ID,x,y) %>% 
    summarise(SWC_min=mean(head(sort(Min),5)),
              SWC_max=mean(tail(sort(Max),5))) %>%
    ungroup()
  SWC=rast(SWC[2:5],crs="epsg:4326")
  
  return(as.data.frame(SWC,xy=TRUE))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Computing Psi_min over the spatial area ####
#' @description Functions used to compute psi_min with different
#' methods, using weighted swc or wieghted psi over horizons
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Compute psi min with weighted swc over horizons
#' 
#' @description function applies Campbell equation to compute minimum soil 
#' potential over spatial extent selected earlier using swc weightedover soil 
#' depth
#' @param texture spatial dataframe to be converted into SpatRaster, containing
#' values of texture at 1x1km resolution
#' @param SWC spatial dataframe to be converted into SpatRaster, containing min
#' and max SWC, computed with 2 different methods, at 9x9km
#' @return Psi_min dataframe, that contains values of Psi_min computed with the
#' 2 different methods used for SWC (average/weighted)
#' 
#' Topsoils Coarse 0.025 0.403 0.0383 1.3774 0.2740 1.2500 60.000
compute_psiweighted <- function(texture,SWC){
  texture=rast(texture,crs="epsg:4326")
  SWC=rast(SWC,crs="epsg:4326")
  if(SWC@ptr$res[1]<texture@ptr$res[1]){
    texture=resample(texture,SWC,method="near")
  }
  if(SWC@ptr$res[1]>texture@ptr$res[1]){
    SWC=resample(SWC,texture)
  }
  psi_min=as.data.frame(c(texture,SWC),
                        xy=TRUE) %>% 
    # filter(clay!=0&sand!=0&oc!=0) %>% 
    mutate(SWC_ratio=SWC_max/SWC_min,
           psi_BG=psi_e*(SWC_max/SWC_min)^(b),
           psi_CB=psi_e*(teta_s/SWC_min)^(b),
           psi_VG=-((((teta_s-teta_r)/(SWC_min-teta_r))^(1/m)-1)^(1/n))*(1/alpha)*9.78*10^(-2)) %>%
    mutate(psi_VG=case_when(SWC_min>teta_s~-0.1,TRUE~psi_VG))
   return(psi_min)
}


#' Compute psi min weighted over horizons
#' 
#' @description function applies Campbell equation to compute minimum soil 
#' potential over spatial extent selected earlier using swc per horizon and 
#' associated parameters extracted from 3D soil database
#' @param SWCtot dataframe of SWCtot cropped to the accurate extent
#' @param hor number of horizons to account for
#' @param timepath month or day
#' @param date_begin "YYYY-MM-DD" one unit timepath before the beginning of the
#' time serie
#' @return Psi_min dataframe, that contains values of Psi_min computed with the
#' 2 different methods used for SWC (average/weighted)

compute_psihor <- function(SWCtot,
                           hor,
                           timepath,
                           begin_date,
                           dir.data="data",
                           dir.file="EU_SoilHydroGrids_1km"){
  
  ## Compute SWC extr per horizons ##

  SWCtot <-rast(SWCtot,crs="epsg:4326")
  for (i in 1:hor){
    assign(paste0("SWC_h",i),
           SWCtot[[grepl(paste0("swvl",i), names(SWCtot))]])
  }
  SWC=mget(ls()[grepl(paste0("SWC_h"),ls())])
  rm(list=ls()[grepl("SWC_h",ls())])
  
  # Compute minimum/max per horizon
  SWC_extr=lapply(SWC,function(x)extr_swc(x,timepath,begin_date))
  
  ## Compute SUREAU ##
  # necessary functions 
  kseriesum <- function(k1,k2) {return(1/(1/k1+1/k2))}
  compute.B_GC <- function (La, Lv, rootRadius) {
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
  return(plot_psi)
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 6 - Compute safety margins ####
#' @description Functions used to explore parametric pedotransfer functions
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute hydraulic and frost safety margins
#' 
#' @description Compute safety margin for selected species
#' @param psi.min df of psimin computed over europe
#' @param tmin df of min temperature over europe
#' @param dir.data data directory
#' @param dir.file directory of species file
#' @return two dataframe for each safety margin, columns are sm of each species
#' 
compute.sm <- function(psimin,
                       tmin,
                       dir.data="data/Species traits/",
                       dir.file="trait.select.csv"){
  species.traits=read.table(file=paste0(dir.data,dir.file),
                            header=TRUE,
                            sep=";",
                            dec=".") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2)))
  
  # hydraulic safety margin
  psimin <- psimin[,c("x","y","psi")]
  for (i in 1:dim(species.traits)[1]){
    psimin$hsm=psimin$psi-(species.traits$P50[i]*1000)
    names(psimin)[names(psimin) == "hsm"] <- species.traits$sp.ind[i]
  }
  
  # fros safety margin
  for (i in 1:dim(species.traits)[1]){
    if(!is.na(species.traits$LT50[i])){
      tmin$fsm=tmin$t.min-(species.traits$LT50[i])
      names(tmin)[names(tmin) == "fsm"] <- species.traits$sp.ind[i]
    }
  }
  return(list(psimin=psimin,tmin=tmin))
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
