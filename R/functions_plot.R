#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing functions that enable to build
#' meaningful plots
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Exploration of climatic/environmental data ####
#' @description Function used to visualize the raw data that are imput into 
#' the computation of Psi_min
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot temporal evolution of SWC
#' 
#' @description Function that enables to plot temporal variation of SWC and to 
#' highlight minima and maxima that are averaged to compute SWC min/max
#' @param point_to_extrat a dataframe which rownames are location and that 
#' contains long/lat of points where SWC is to be explored. Default is set for 5
#' climate-contrasted location in France
#' @param moy numbers of maxima/minima points to be highlighted. Default is 5 as
#' set in the computation of SWC
#' @param SWC spatial dataframe containing SWC over a time period 
#' @param method method used to compute SWC over the different layers
#' @return a graph of SWC temporal evolution for each location
#' 
chronology_swc <- function(point_to_extract=data.frame(Luberon=c(5.387792,43.814892),
                                                     Landes= c(-0.438943,44.427567), #Landes
                                                     Maison= c(-1.844066,47.649440), #la maison
                                                     Chamrousse= c(5.860739,45.092439), #Chamrousse
                                                     GrandEst=c(5.569744,48.475034)), #Grand-Est
                         moy=5,
                         SWC,
                         method="Mean") {
  ## Format points to extract
  loc= colnames(point_to_extract)
  row.names(point_to_extract)=c("Lon","Lat")
  point_to_extract=t(point_to_extract)
  row.names(point_to_extract) <- NULL
  
  ## Format SWC_stack
  SWC_time=terra::rast(SWC,crs="epsg:4326")
  SWC_time=as.data.frame(extract(SWC_time,point_to_extract),
                           row.names = loc) %>% 
    tibble::rownames_to_column(var="Loc") %>% 
    pivot_longer(cols=colnames(.)[c(-1)]) %>% 
    separate(name,sep="_",c("drop","time")) %>% 
    select(-drop) %>% 
    mutate(time=as.numeric(time)) %>% 
    group_by(Loc) %>% 
    arrange(value) %>% 
    mutate(rnk=row_number()) %>% 
    mutate(maxima=dplyr::case_when(rnk<moy+1~"Min",
                                   rnk>(max(rnk)-moy)~"Max")) %>% 
    ungroup() %>% 
    arrange(Loc,time)
  SWC_time %>% 
    ggplot(aes(time,value,maxima))+
    geom_point(aes(color=maxima))+
    geom_line()+
    facet_wrap(~Loc)+
    labs(y=method)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Species distribution and hydraulic safety margin ####
#' @description Function used to compute HSM and to compare it with species
#' distribution
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Overlaying species distribution and hydraulic safety margin
#' @description functions that plots a polygon of approximate species distribution
#' extracted from EuForGen, and overlays HSM represented in a binary way
#' @param Dir.Data directory of data of the project
#' @param Dir.p50 files that contains mean P50 per species of interest, with files
#' extention
#' @param europe SpatialPolygonsDataFrame of the spatial extent, default is Europe
#' see before
#' @param psi_min spatial dataframe of psi_min
#' @return plots per species of interest
HSM_distribution <- function(Dir.Data="data",
                             Dir.p50="p50select.csv",
                             europe,
                             psi_min){
  P50_df=read.csv2(file.path(Dir.Data,Dir.p50))
  Species_distb=apply(P50_df,1,function(x){
    P50_sp=as.numeric(x[2])
    spdistrib=read_sf(dsn=file.path(Dir.Data,"chorological_maps_dataset",as.character(x[1]),"shapefiles"),
                      layer=as.character(x[3]))
    if(is.na(x[4])==FALSE){
      spdistrib_sup=read_sf(dsn=file.path(Dir.Data,"chorological_maps_dataset",x[1],"shapefiles"),
                            layer=as.character(x[4]))
      spdistrib=rbind(spdistrib,spdistrib_sup)
    }
    spdistrib=st_crop(spdistrib,sf::st_bbox(europe))
    
    graph=psi_min %>%
      select(x,y,Psi_w_min) %>%
      mutate(HSM=Psi_w_min-P50_sp) %>% #P50(F_syl)=-3.15
      mutate(HSM_bin=as.factor(case_when(HSM>0~1,
                                         HSM<=0~0))) %>% 
      select(x,y,HSM_bin) %>% 
      filter(is.na(HSM_bin)==FALSE) %>% 
      ggplot() + # create a plot
      geom_raster( aes(x = x, y = y, fill = HSM_bin)) + # plot the raw data
      geom_sf(data=spdistrib,fill=alpha("grey",0.6))+
      theme_bw() +
      scale_fill_discrete(na.value="transparent")+
      labs(x = "Longitude", y = "Latitude",fill=paste0(as.character(x[1])," presence")) 
    return(graph)
  })
  return(Species_distb)
  
}

#' Extract values from species distribution
#' @description function that extracts Psi_min or HSM into the distribution of a 
#' species and plot the distribution of the Psi_min, its quantiles and the P50
#' @return a distribution of HSM/Psi_min per species of interest 
#' 
Psi_min_distrib <- function(Dir.Data="data",
                            Dir.p50="p50select.csv",
                            europe,
                            psi_min){
  P50_df=read.csv2(file.path(Dir.Data,Dir.p50))
  Species_distb=apply(P50_df,1,
                      function(x){
                        P50_sp=as.numeric(x[2])
                        spdistrib=read_sf(dsn=file.path(Dir.Data,"chorological_maps_dataset",as.character(x[1]),"shapefiles"),
                                          layer=as.character(x[3]))
                        if(is.na(x[4])==FALSE){
                          spdistrib_sup=read_sf(dsn=file.path(Dir.Data,"chorological_maps_dataset",x[1],"shapefiles"),
                                                layer=as.character(x[4]))
                          spdistrib=rbind(spdistrib,spdistrib_sup)
                        }
                        return(vect(spdistrib))
                        }
                      )
  Distrib=sapply(1:length(Species_distb),function(x)
    terra::extract(rast(psi_min[,c(1,2,12)]),Species_distb[[x]],xy=TRUE)
    )
    
    
  
}
abies=vect(st_read(dsn="data/chorological_maps_dataset/Abies alba/shapefiles",layer = "Abies_alba_plg"))
distrib %>% 
  filter(Psi_w_min>-30) %>% 
  ggplot(aes(Psi_w_min))+
  geom_density()+geom_vline(xintercept=-3.79,colour="red") + 
  stat_summary(geom = "vline",
               orientation = "y",
               # y is a required aesthetic, so use a dummy value
               aes(y = 1, xintercept = after_stat(x)),
               fun = function(x) { 
                 quantile(x, probs = c(0.3, 0.9))
               }
  )







#### Comparison of Terraclimate dataset and ERA5-land ####
##########################################################

# terraclimate from 1958 to 2021
terraclimate_sub=rast(terraclimate[,c(1,2,3,3+120,3+240,3+360,3+480)],crs="epsg:4326")
SWC_sub=rast(SWC[,c(1,2,99,99+120,99+240,99+360,99+480)],crs="epsg:4326")
terraclimate_sub=crop(terraclimate_sub,SWC_sub)
SWC_sub=project(SWC_sub,terraclimate_sub)
comp=as.data.frame(c(SWC_sub,terraclimate_sub),xy=TRUE)
comp %>% 
  ggplot(aes(swvl1_577*289,soil_481))+geom_bin2d()+scale_fill_continuous(type = "viridis") +
  theme_bw()+xlab("ERA5")+ylab("TerraClimate")

