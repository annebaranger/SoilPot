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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Downloading climatic and pedologic data ####
#' @description Functions used to download/load required climatic and pedologic
#'  data
#' @authors Anne Baranger (INRAE - LESSEM)
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
#' @param Dir.Data directory of project data
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
                     Dir.Data,
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
  lapply(variables,function(x){Dir.Var=file.path(Dir.Data,x)
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Computing Psi_min over the spatial area ####
#' @description Functions used to compute required variables to apply Campbell 
#' equations and get Psi_min 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute texture over a specified spatial extent
#' 
#' @description Function that load ESDAC texture data from directory
#' @note Texture data can be downloaded on ESDAC website. They encompass texture
#' from topsoil and subsoil at 1x1km resolution
#' @param Dir.Data directory of data of the project
#' @param Dir.Soil directory of soil data downloaded from ESDAC
#' @param europe SpatialPolygonsDataFrame of the spatial extent 
#' @return spatial dataframe (SpatRaster converted to dataframe) of clay and 
#' sand content average among topsoil and subsoil, cropped to spatial extent
texture_data <- function(Dir.Data,Dir.Soil,europe){
  # get spatial extent
  rast_model <- vect(europe)
  # Clay 
  clay <- c(rast(file.path(Dir.Data,Dir.Soil,"STU_EU_T_CLAY.rst")),rast(file.path(Dir.Data,Dir.Soil,"STU_EU_S_CLAY.rst")))
  crs(clay) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  clay <- crop(project(clay,"epsg:4326"),rast_model)
  clay <- app(clay,mean)
  # Sand 
  sand <- c(rast(file.path(Dir.Data,Dir.Soil,"STU_EU_T_SAND.rst")),rast(file.path(Dir.Data,Dir.Soil,"STU_EU_S_SAND.rst")))
  crs(sand) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  sand <- crop(project(sand,"epsg:4326"),rast_model)
  sand <- app(sand,mean) 
  # Texture
  texture <- c(clay,sand)
  names(texture) <- c("clay","sand")
  texture <- as.data.frame(texture,xy=TRUE)
  return(texture)
}


#' Compute soil volumetric content
#' 
#' @description Using ERA5-land SWC, the function compute minimum and maximum SWC
#' of each cell of the area, over the time period downloaded. It uses two 
#' computation method : an average over the 4 soil layers, or a weighted average.
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param Dir.Data directory of data of the project
#' @param texture spatial dataframe of soil textures in the area considered
#' @param europe SpatialPolygonsDataFrame of the spatial extent 
#' @return list of 3 spatial dataframe : first with min/max SWC computed with 2
#'  methods, second with SWC over the time period with weighted method, last is 
#'  the same with averaged method
#'  
# Dir.Data="data"
# Variables=c("Volumetric_soil_water_layer_1","Volumetric_soil_water_layer_2","Volumetric_soil_water_layer_3","Volumetric_soil_water_layer_4")
volumetric_content <- function(Dir.Data,texture,europe){ 
  # get texture raster
  texture <- rast(texture,crs="epsg:4326")

  # load SWC 
  SWCtot <- rast("data/swc-1950-2021.nc")
  SWC1 <- crop(subset(SWCtot,1:864),europe)
  SWC2 <- crop(subset(SWCtot,865:1728),europe)
  SWC3 <- crop(subset(SWCtot,1729:2592),europe)
  SWC4 <- crop(subset(SWCtot,2593:3456),europe)
  
  #weighted mean or absolute mean over the period downloaded
  if (sum(sapply(c("SWC1","SWC2","SWC3","SWC4"),function(x)exists(x)))==4){
      SWC_w<-(7*SWC1+21*SWC2+72*SWC3+189*SWC4)/289
      SWC_w_min<-app(SWC_w,function(x)mean(head(sort(x),5)))
      SWC_w_max<-app(SWC_w,function(x)mean(tail(sort(x),5)))
      SWC_m<-(SWC1+SWC2+SWC3+SWC4)/4
      SWC_m_min<-app(SWC_m,function(x)mean(head(sort(x),5)))
      SWC_m_max<-app(SWC_m,function(x)mean(tail(sort(x),5)))
      #Build SWC dataframe
      SWC <- plyr::join_all(lapply(c("SWC_w_min","SWC_w_max","SWC_m_min","SWC_m_max"),
                                function(x){
                                  assign(paste0(x,"_df"),
                                         as.data.frame(project(get(x),texture),
                                                       xy=TRUE))
                                }
      ),by=c("x","y"))
      colnames(SWC) <- c("x","y","SWC_w_min","SWC_w_max","SWC_m_min","SWC_m_max")
      
  }
  return(list(SWC,as.data.frame(SWC_w,xy=TRUE),as.data.frame(SWC_m,xy=TRUE)))
  
}


#' Compute minimum soil potential
#' 
#' @description function applies Campbell equation to compute minimum soil 
#' potential over spatial extent selected earlier
#' @param texture spatial dataframe to be converted into SpatRaster, containing
#' values of texture at 1x1km resolution
#' @param SWC spatial dataframe to be converted into SpatRaster, containing min
#' and max SWC, computed with 2 different methods, at 9x9km
#' @return Psi_min dataframe, that contains values of Psi_min computed with the
#' 2 different methods used for SWC (average/weighted)
soil_potential <- function(texture,SWC){
   # Format rasters
  texture <- rast(texture,crs="epsg:4326")
  SWC <- project(rast(SWC,crs="epsg:4326"),texture)
  #Conversion to potential
   texture_pot <- data.frame(texture=c(5,4,3,2,1),
                          Psi_e=c(-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
                          b=c(2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933))
   
   psi_min <- as.data.frame(texture,xy=TRUE) %>% 
     #categorize texture
     mutate(texture=case_when(clay>60~1,
                              clay<18&sand>65~5,
                              clay>18&clay<35&sand>15~4,
                              clay>18&sand>15&sand<65~4,
                              clay<35&sand>15~3,
                              clay>35&clay<60~2)) %>% 
     left_join(texture_pot,by=c("texture")) %>% 
     left_join(as.data.frame(SWC,xy=TRUE),by=c("x","y")) %>% 
     mutate(Psi_w_min=Psi_e*(SWC_w_max/SWC_w_min)^(b),
            SWC_w_ratio=SWC_w_max/SWC_w_min,
            Psi_m_min=Psi_e*(SWC_m_max/SWC_m_min)^(b),
            SWC_m_ratio=SWC_m_max/SWC_m_min,
            Psi_w_med=texture_pot$Psi_e[3]*(SWC_w_max/SWC_w_min)^(texture_pot$b[3]),
            Psi_m_med=texture_pot$Psi_e[3]*(SWC_m_max/SWC_m_min)^(texture_pot$b[3]))
   #psi_min=rast(psi_min)
   return(psi_min)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Annexe - Useful chunks that were deleted  ####
#' @description few chunks that could be re-used
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#load rasters of variables of interest
# for (x in (1:length(Abv))){
#   assign(Abv[x],rast(list.files(file.path(Dir.Data,Variables[x]),full.names = TRUE)[1]))
# }