#### Packages required ####
###########################

package=c("KrigR", "ggplot2","tidyr","viridis","rgdal","raster","rosm","terra","dplyr","gdalUtils")
lapply(package, require, character.only = TRUE)


#### Data downloading from ERA5-land ####
#########################################

### Directories ###
Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data
Dir.Shapes <- file.path(Dir.Data, "Shapes") # folder path for shapefiles
Dirs <- sapply(c(Dir.Data, Dir.Shapes), function(x) if (!dir.exists(x)) dir.create(x))

### DL regions masks, function from KrigR package ###
if (!file.exists(file.path(Dir.Shapes, "StateMask.zip"))) { # if not downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip",
                destfile = file.path(Dir.Shapes, "StateMask.zip")
  ) # download cultural vector
  unzip(file.path(Dir.Shapes, "StateMask.zip"), exdir = Dir.Shapes) # unzip data
}

### DL ERA5 for a specific region ###
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
  ## Read regions masks
  StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE) 
  
  ## Find the extent of seleted countries
  if(is.null(extent)){
    position=unlist(lapply(countries,function(x) which(StateMask$geonunit==x)))
    extent=StateMask[position,]
    # extent= extent(State_Shp)
  }
  
  ## Dl ERA5 for each variable
  lapply(variables,function(x){
    Dir.Var=file.path(Dir.Data,x)
    dir.create(Dir.Var) # Create directory for each variable
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

## Dl monthly SWC from ERA5-land for France, between 2000 and 2015
era_data(Dir.Data="data",
         Dir.Shapes="data/Shapes",
         countries = "France",
         variables=c("Volumetric_soil_water_layer_1","Volumetric_soil_water_layer_2","Volumetric_soil_water_layer_3","Volumetric_soil_water_layer_4"),
         date_start="2000-01-01",
         date_end="2015-01-01",
         time_step="month")



#### Download Soildgrid texture ####
####################################
StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE) 
countries=c("France")
position=unlist(lapply(countries,function(x) which(StateMask$geonunit==x)))
extent=StateMask[position,]
igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection
sg_url="/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"

bb=st_as_sfc(st_bbox(extent,crs="+proj=longlat +datum=WGS84 +no_defs +type=crs"))
bb=st_bbox(st_transform(bb,igh))

### Download clay content for 0-5cm
gdal_translate(paste0(sg_url,'clay/clay_0-5cm_mean.vrt'),
               "./clay_0-5cm_mean.tif",
               of="VRT",tr=c(250,250),
               projwin=bb,
               projwin_srs =igh,
               verbose=TRUE)

