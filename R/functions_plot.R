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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Species distribution and both safety margins ####
#' @description Function used to compute HSM and to compare it with species
#' distribution
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Overlaying species distribution and hydraulic safety margin
#' @description functions that plots a polygon of approximate species distribution
#' extracted from EuForGen, and overlays HSM represented in a binary way
#' @param dir.data directory of data of the project
#' @param dir.p50 files that contains mean P50 per species of interest, with files
#' extention
#' @param europe SpatialPolygonsDataFrame of the spatial extent, default is Europe
#' see before
#' @param psi_min spatial dataframe of psi_min
#' @return plots per species of interest
plot.sftym <- function(dir.data="data",
                             dir.file="Species traits/trait.select.csv",
                             europe,
                             safety.margins.day){
  # Load and reshape species list and shapefiles list
  species.trait <- read.table(file=file.path(dir.data,dir.file),
                              header=TRUE,
                              sep=";",
                              dec=".") %>%
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2))) 

  # Aggregate safety margins raster for speed matters
  tmin.simpl=as.data.frame(aggregate(rast(safety.margins.day$tmin,crs="epsg:4326"),
                                   12,
                                   fun="mean"),
                         xy=TRUE)
  psimin=safety.margins.day$psimin
  rm(safety.margins.day)
  
  # Build and save graphs
  
  ## Load shapefiles
  for (i in 1:dim(species.trait)[1]){
    tryCatch(
      {
        spdistrib=read_sf(dsn=file.path(dir.data,
                                        "chorological_maps_dataset",
                                        species.trait$species.binomial[i],
                                        "shapefiles"),
                          layer=species.trait$file[i])
        if(is.na(species.trait$extra.file[i])==FALSE){
          spdistrib_sup=read_sf(dsn=file.path(dir.data,
                                              "chorological_maps_dataset",
                                              species.trait$species.binomial[i],
                                              "shapefiles"),
                                layer=species.trait$extra.file[i])
          spdistrib=st_union(rbind(spdistrib,spdistrib_sup))
        }
        spdistrib=st_crop(spdistrib,sf::st_bbox(europe))
        
        ## Build graphs
        HSM <- psimin %>% 
          pivot_longer(cols=colnames(.)[c(-1,-2)],
                       names_to = "Species",
                       values_to = "HSM") %>% 
          filter(Species==species.trait$sp.ind[i]) %>% 
          mutate(HSM=cut(HSM,
                         breaks=c(-Inf, -10000, -5000,-3000,0,500,1000,5000,Inf),
                         labels=c("<-10MPa","-10<HSM<-5MPa","-5<HSM<-3MPa",
                                  "-3<HSM<0MPa","0/0.5","0.5/1","1/5",">5MPa"))) %>% 
          ggplot() +
          geom_tile(aes(x=x,y=y,fill=HSM))+
          geom_sf(data=spdistrib,
                  fill=alpha("grey",0.6),
                  lwd = 0)+
          coord_sf()+
          theme_bw() +
          theme(axis.title=element_blank(),
                legend.key.size = unit(0.5,"cm"),
                legend.text = element_text(size=8),
                legend.title = element_text(size=9))+
          labs(fill=paste0("HSM (MPa), \nSpecies: ",species.trait$species.binomial[i],
                           "\nP50=",species.trait$P50[i]))+
          scale_fill_brewer(palette="RdYlBu")
        ggsave(HSM,
               filename = paste0(species.trait$sp.ind[i],"HSM.png"),
               path="figs/",
               device="png",
               scale=2)
        
        FSM <- tmin.simpl %>% 
          pivot_longer(cols=colnames(.)[c(-1,-2)],
                       names_to = "Species",
                       values_to = "FSM") %>% 
          filter(Species==species.trait$sp.ind[i]) %>% 
          mutate(FSM=cut(FSM,
                         breaks=c(-Inf, -15,-10, -5,0,5,10,15,20,25,30,35,Inf),
                         labels=c("<-15", "-15<FSM<-10", "-10/-5","-5/0","0/5","5/10","10/15","15/20","20/25","25/30","30/35",">35"))) %>% 
          ggplot() +
          geom_tile(aes(x=x,y=y,fill=FSM))+
          geom_sf(data=spdistrib,
                  fill=alpha("grey",0.6),
                  lwd = 0)+
          coord_sf()+
          theme_bw() +
          theme(axis.title=element_blank(),
                legend.key.size = unit(0.5,"cm"),
                legend.text = element_text(size=8),
                legend.title = element_text(size=9))+
          labs(fill=paste0("FSM (MPa), \nSpecies: ",species.trait$species.binomial[i],
                           "\nLT50=",round(species.trait$LT50[i],digit=1)))+
          scale_fill_brewer(palette="RdYlBu")
        ggsave(FSM,
               filename = paste0(species.trait$sp.ind[i],"FSM.png"),
               path="figs/",
               device="png",
               scale=2)
      },
      error=function(e){print("File not loaded")}
    )
  }
  }



#' Extract values from species distribution
#' @description function that extracts Psi_min or HSM into the distribution of a 
#' species and plot the distribution of the Psi_min, its quantiles and the P50
#' @return a distribution of HSM/Psi_min per species of interest 
#' 
Psi_min_distrib <- function(dir.data="data",
                            dir.p50="p50select.csv",
                            europe,
                            psi_min){
  P50_df=read.csv2(file.path(dir.data,dir.p50))
  Species_distb=apply(P50_df,1,
                      function(x){
                        P50_sp=as.numeric(x[2])
                        spdistrib=read_sf(dsn=file.path(dir.data,"chorological_maps_dataset",as.character(x[1]),"shapefiles"),
                                          layer=as.character(x[3]))
                        if(is.na(x[4])==FALSE){
                          spdistrib_sup=read_sf(dsn=file.path(dir.data,"chorological_maps_dataset",x[1],"shapefiles"),
                                                layer=as.character(x[4]))
                          spdistrib=rbind(spdistrib,spdistrib_sup)
                        }
                        return(terra::vect(spdistrib))
                        }
                      )
  Distrib=sapply(1:length(Species_distb),function(x){
    dist=terra::extract(rast(psi_min[,c(1,2,12)]),Species_distb[[x]],xy=TRUE) %>% 
      select(-ID) %>% 
      relocate(c("x","y"))
    return(dist)
    }
  )
    
    
  
}
# abies=terra::vect(st_read(dsn="data/chorological_maps_dataset/Abies alba/shapefiles",layer = "Abies_alba_plg"))
# distrib %>% 
#   filter(Psi_w_min>-30) %>% 
#   ggplot(aes(Psi_w_min))+
#   geom_density()+geom_vline(xintercept=-3.79,colour="red") + 
#   stat_summary(geom = "vline",
#                orientation = "y",
#                # y is a required aesthetic, so use a dummy value
#                aes(y = 1, xintercept = after_stat(x)),
#                fun = function(x) { 
#                  quantile(x, probs = c(0.3, 0.9))
#                }
#   )
# 
# 
# 
# 
# 
# 
# 
# #### Comparison of Terraclimate dataset and ERA5-land ####
# ##########################################################
# 
# # terraclimate from 1958 to 2021
# terraclimate_sub=rast(terraclimate[,c(1,2,3,3+120,3+240,3+360,3+480)],crs="epsg:4326")
# SWC_sub=rast(SWC[,c(1,2,99,99+120,99+240,99+360,99+480)],crs="epsg:4326")
# terraclimate_sub=crop(terraclimate_sub,SWC_sub)
# SWC_sub=project(SWC_sub,terraclimate_sub)
# comp=as.data.frame(c(SWC_sub,terraclimate_sub),xy=TRUE)
# comp %>% 
#   ggplot(aes(swvl1_577*289,soil_481))+geom_bin2d()+scale_fill_continuous(type = "viridis") +
#   theme_bw()+xlab("ERA5")+ylab("TerraClimate")
# 
