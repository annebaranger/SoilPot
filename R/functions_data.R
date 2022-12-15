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
  terraclimate <- as.data.frame(rast("data/Terraclimate/agg_terraclimate_soil_1958_CurrentYear_GLOBE.nc"),
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
  
  return(list(texture_raw=texture_raw,
              texture=texture))
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
                    dir.file
                    ){
  SWCtot <-rast(file.path(dir.data,dir.file))
  #europe <- vect(st_as_sf(europe))
  #SWCtot <- crop(SWCtot,europe,mask=TRUE)
  return(as.data.frame(SWCtot,xy=TRUE))
}


# nc.data <-   read_ncdf("data/swc-1950-2021.nc",proxy=FALSE)
# sf_use_s2(FALSE)
# nc.stars.crop <-  nc.data[europe]
# plot(nc.data[,,,1:3])
# write_stars(nc.stars.crop,dsn="data/test.tif",layer="swvl1")
# # Input nc file
# nc.file <- file.path(dir.data,dir.file)
# # read nc data
# nc.data <- read_ncdf(nc.file)
# 
# # Read mask coordinates
# coordenades.poligon <- read_csv("coordenades_poligon.csv")
# colnames(coordenades.poligon) <- c("lon","lat")
# 
# # Build sf polygon to crop data
# europe <- st_as_sf(europe)
# bb_eu <- st_bbox(europe)
# 
# # Crop data
# nc.stars.crop <- st_crop(nc.data,europe,
#                          crop = TRUE,
#                          epsilon = sqrt(.Machine$double.eps),
#                          as_points = all(st_dimension(europe) == 2, na.rm = TRUE))
# 
# write_stars(nc.data[europe],dsn="data/test.nc")

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
weight_swc <- function(SWCtot,
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
#' @description Using monthly (or daily but not used here) time series of weighted
#' soil water content, the function compute the min and max of timeseries
#' @note ERA5-land data needs to be downloaded and weighted prior to applying 
#' the function
#' @param SWC_t dataframe of timeserie of swc
#' @param timepath month or day
#' @param date_begin "YYYY-MM-DD" one unit timepath before the beginning of the
#' time serie
#' @return dataframe with min/max
#' date_begin='1949-12-01'
extr_swc <- function(SWC_t,timepath,date_begin,europe) {
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
  europe <- vect(st_as_sf(europe))
  SWC <- crop(SWC,europe,mask=TRUE)
  return(as.data.frame(SWC,xy=TRUE))
}


#' Load min computed externally with python
#' 
#' @description Load and reshape SWC data computed with python program, per horizon
#' @note Python code need to be run prior executing this function
#' @param dir.data
#' @param dir.var
#' @param vars list of vars corresponding to each horizons computed
#' @return dataframe with min swc per horizon

load_swc <- function(dir.data="data",
                     dir.file="ERA5-land/swcd-1950-2021-",
                     vars=c("layer1","layer2","layer3","layer4")){
  rast.model=rast(paste0(dir.data,"/",dir.file,vars[1],".nc"))[[1]]
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
  names(rast.swc)=vars
  return(as.data.frame(rast.swc,xy=TRUE))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Computing Psi_min over the spatial area ####
#' @description Functions used to compute psi_min with different
#' methods, using weighted swc or wieghted psi over horizons
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Compute psi min with weighted swc over horizons
#' 
#' @description function applies CP and VG equation to compute minimum soil potential 
#' over weighted horizons. Using daily or monthly data, but not taking into account 
#' the different horizons
#' @param texture spatial dataframe to be converted into SpatRaster, containing
#' values of texture at 1x1km resolution
#' @param SWC spatial dataframe to be converted into SpatRaster, containing min
#' and max SWC, computed with 2 different methods, at 9x9km
#' @return Psi_min dataframe, that contains values of Psi_min computed with the
#' 2 different methods used for SWC (average/weighted)
#' 
#' Topsoils Coarse 0.025 0.403 0.0383 1.3774 0.2740 1.2500 60.000
compute_psiweighted <- function(texture,SWC,
                                file.output){
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
  write.csv(psi_min,file=file.output)
  
  return(psi_min)
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
  rast.temp=rast(rast.temp,crs="epsg:4326")
  rast.fdg=rast(rast.fdg,crs="epsg:4326")
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
  df.LT50.2022 <- read.csv2("data/Species traits/camp_LT50_DB.csv") # made by code in "Analysis new campaign data.Rmd" : output file
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
    rm(df.sp)
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
  }}


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
  # df.P50.lit <- read.csv2("data/Species traits/BdD_full_by_species.csv")
  # df.P50.select <- df.P50.lit %>% 
  #   select(colnames(df.P50.lit)[!grepl("LT",colnames(df.P50.lit))]) %>% 
  #   select(-Dmean,-Psi_min_MD,-WD,-Plot_continent,-MHZ,-nplot,-ntree,-D99,-score,-X.gs90,-X.tlp,-densityMEANcor,-Exo,-genus,-GBIF_species,-class) %>% 
  #   filter(New_continent=="E_U") %>% 
  #   filter(!is.nan(P50tot))
  df.p50.lit <- read.csv2("data/Species traits/2022_10_20_cleaned_xft.csv") %>% 
    select(uniqueID,lat,long,XFT.database,Group,Family,Genus,Species,
           Developmental.stage,Growth.form,P50,P50.SD,P12,P88,Curve,Equation,
           psi.min.predawn..MPa.,psi..min.midday..MPa.) %>% 
    mutate(species.binomial=paste(Genus,Species)) %>% 
    filter(Curve=="S") %>%
    filter(Growth.form=="T") %>% 
    group_by(Group,Genus,Species,species.binomial) %>% 
    summarise(P50.mu=mean(as.numeric(P50),na.rm=TRUE),
              P50.sd=sd(as.numeric(P50),na.rm=TRUE),
              P88.mu=mean(as.numeric(P88),na.rm=TRUE),
              P88.sd=sd(as.numeric(P88),na.rm=TRUE)) %>% 
    ungroup()
  write.csv(df.p50.lit,"output/df_P50_filtered.csv",row.names = FALSE)
  return(df.p50.lit)
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
    left_join(df.P50,by=c("Species"="species.binomial"))
  # df.traits.raw <- full_join(df.P50,df.LT50,by="Species")
  # df.traits.missing <- df.traits.raw %>% 
  #   filter(is.na(P50tot) | is.na(LT50.mean))
  df.traits <- df.traits.raw %>% 
    filter(!is.na(P50.mu)) %>% 
    filter(!is.na(LT50.mean)) %>% 
    select(Species,Group,P50.mu,P50.sd,P88.mu,P88.sd,LT50.mean,LT50.sd,data.quality) %>% #MATmean,MAPmean,Leaf_phenology,
    unique() %>% 
    rename(`species.binomial`="Species") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2)),
           species.name=str_replace(species.binomial," ","_")) %>% 
    relocate(c("species.name","sp.ind"),.after="species") %>% 
    mutate(p.trait=case_when(Group=="gymnosperm"~"P50",
                             Group=="angiosperm"&!is.na(P88.mu)~"P88",
                             Group=="angiosperm"&is.na(P88.mu)~"P50"),
           PX.mu=case_when(p.trait=="P88"~P88.mu,
                           p.trait=="P50"~P50.mu),
           PX.sd=case_when(p.trait=="P88"~P88.sd,
                           p.trait=="P50"~P50.sd))
  
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