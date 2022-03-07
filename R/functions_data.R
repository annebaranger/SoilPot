#### Data downloading from ERA5-land ####
#########################################


era_data <- function(API_User="124078",
                     API_Key="ab883de0-39cd-4452-a1f0-7cd45ab252b0",
                     Dir.Data,
                     Dir.Shapes, # directory of countries masks
                     extent=NULL,
                     countries,
                     variables,
                     date_start,
                     date_end,
                     time_step){
  StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE)
  
  if(is.null(extent)){
    position=unlist(lapply(countries,function(x) which(StateMask$geonunit==x)))
    extent=StateMask[position,]
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


#### DL soilgrid ####
#####################
# function not working with updated version of gdal
# library(rgdal)
# library(gdalUtils)
# 
# bb=st_as_sfc(st_bbox(SWC1,crs="+proj=longlat +datum=WGS84 +no_defs +type=crs"))
# bb=st_bbox(st_transform(bb,igh))
# igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection
# sg_url="/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"
# 
# gdal_translate(paste0(sg_url,'clay/clay_0-5cm_mean.vrt'),
#                "./crop_roi_igh_r.vrt",
#                of="VRT",tr=c(250,250),
#                projwin=bb,
#                projwin_srs =igh,
#                verbose=TRUE)



#### Function to compute a stack for texture ####
#################################################
# Dir.Data="data"
# Dir.Soil="STU_EU_Layers"
texture_data <- function(Dir.Data,Dir.Soil){
  # To get spatial extent
  rast_model=rast("data/Volumetric_soil_water_layer_1/Volumetric_soil_water_layer_1_2000-01-01_2015-01-01_month.nc")
  # Clay process
  clay=c(rast(file.path(Dir.Data,Dir.Soil,"STU_EU_T_CLAY.rst")),rast(file.path(Dir.Data,Dir.Soil,"STU_EU_S_CLAY.rst")))
  crs(clay)="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  clay=crop(project(clay,"epsg:4326"),rast_model)
  clay=app(clay,mean)
  # Sand process
  sand=c(rast(file.path(Dir.Data,Dir.Soil,"STU_EU_T_SAND.rst")),rast(file.path(Dir.Data,Dir.Soil,"STU_EU_S_SAND.rst")))
  crs(sand)="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
  sand=crop(project(sand,"epsg:4326"),rast_model)
  sand=app(sand,mean) 
  # Texture
  texture=c(clay,sand)
  names(texture)=c("clay","sand")
  texture=as.data.frame(texture,xy=TRUE)
  return(texture)
}



#### Function to compute soil volumetric content ####
#####################################################
##for each cell, using SWC from 4 soil layers over 15 years 
# Dir.Data="data"
# Variables=c("Volumetric_soil_water_layer_1","Volumetric_soil_water_layer_2","Volumetric_soil_water_layer_3","Volumetric_soil_water_layer_4")
# Abv=c("SWC1","SWC2","SWC3","SWC4")

volumetric_content <- function(Dir.Data,Variables,Abv,texture){
  # get texture raster
  texture=rast(texture,crs="epsg:4326")
  #load rasters of variables of interest
  for (x in (1:length(Abv))){
    assign(Abv[x],rast(list.files(file.path(Dir.Data,Variables[x]),full.names = TRUE)[1]))
  }
  #weighted mean or absolute mean over 15 years
  if (sum(sapply(c("SWC1","SWC2","SWC3","SWC4"),function(x)exists(x)))==4){
      SWC_w<<-(7*SWC1+21*SWC2+72*SWC3+189*SWC4)/289
      SWC_w_min<<-app(SWC_w,function(x)mean(head(sort(x),5)))
      SWC_w_max<<-app(SWC_w,function(x)mean(tail(sort(x),5)))
      SWC_m<<-(SWC1+SWC2+SWC3+SWC4)/4
      SWC_m_min<<-app(SWC_m,function(x)mean(head(sort(x),5)))
      SWC_m_max<<-app(SWC_m,function(x)mean(tail(sort(x),5)))
      #return(list(SWC_w,SWC_w_min,SWC_w_max,SWC_m,SWC_m_min,SWC_m_max))
  #Build SWC dataframe
      #SWC_time=apply(c(SWC_w,SWC_m),function(x)as.data.frame(x,xy=TRUE))
      SWC=plyr::join_all(lapply(c("SWC_w_min","SWC_w_max","SWC_m_min","SWC_m_max"),
                                function(x){
                                  assign(paste0(x,"_df"),
                                         as.data.frame(project(get(x),texture),
                                                       xy=TRUE))
                                }
      ),by=c("x","y"))
      colnames(SWC)=c("x","y","SWC_w_min","SWC_w_max","SWC_m_min","SWC_m_max")
      
  }
  return(list(SWC,as.data.frame(SWC_w,xy=TRUE),as.data.frame(SWC_m,xy=TRUE)))
  
}


#### Function to compute soil potential ####
############################################
soil_potential <- function(texture,SWC){
   # Format rasters
  texture=rast(texture,crs="epsg:4326")
  SWC=project(rast(SWC,crs="epsg:4326"),texture)
  #Conversion to potential
   texture_pot=data.frame(texture=c(5,4,3,2,1),
                          Psi_e=c(-0.790569415,-0.9128709202,-1.5811388301,-1.889822365,-5.9761430467),
                          b=c(2.6411388301,3.3057418584,4.3822776602,6.5796447301,14.9522860933))
   
   psi_min=as.data.frame(texture,xy=TRUE) %>% 
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
 