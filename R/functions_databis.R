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
  
  return(as.data.frame(chelsabio6,xy=TRUE))
}



#' Soil volumetric content scaled for targetted area, daily data
#' 
#' @description crop era5land SWC to the good extent
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param dir.data directory of data of the project
#' @param dir.file directory of data file
#' @param europe SpatialPolygonsDataFrame of the spatial extent 
#' @return dataframe with temporal series of SWC
#' 
#' ext_fresp=ext(-9.95, 10.05, 34.95, 52.05)

get_SWCd <- function(dir.data="data",
                    dir.file="swcd1-3_2011-2021_FrEsp.nc",
                    europe,
                    ext=NA,
                    crop=FALSE){
  SWC_dtot=rast(file.path(dir.data,dir.file))
  if (crop==TRUE){
    if(is.na(ext)==TRUE){
      europe <- vect(st_as_sf(europe))
      SWC_dtot=crop(SWC_dtot,europe,mask=TRUE)
    }
    else{
      SWC_dtot=crop(SWC_dtot,ext,mask=TRUE)
    }
  }
  return(as.data.frame(SWC_dtot,xy=TRUE))
}



#' Compute minimum soil volumetric content per horizon with daily timeserie
#' 
#' @description Using ERA5-land SWC, the function compute minimum and maximum SWC
#' of each cell of the area, and for each horizon over the time period downloaded.
#' It computes min/max for each of the 3 first horizons
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param SWC_dtot dataframe with timeseries per horizons
#' @param horizons SpatialPolygonsDataFrame of the spatial extent 
#' @return list of dataframes of min/max per horizons
#'  
# dir.data="data"

volumetric_content <- function(SWC_dtot,
                               horizons=3){
  # create a list of time series per horizons 
  SWC_dtot <-rast(SWC_dtot,crs="epsg:4326")
  for (i in 1:horizons){
    assign(paste0("SWC_dh",i),
           SWC_dtot[[grepl(paste0("swvl",i), names(SWC_dtot))]])
  }
  SWCdh=mget(ls()[grepl(paste0("SWC_dh"),ls())])
  rm(list=ls()[grepl("SWC_dh",ls())])
  
  # Compute minimum/max per horizon
  SWC_mx=lapply(SWCdh,function(x){
    SWC_df=setDT(as.data.frame(x,xy=TRUE)) 
    time=as.character(as.Date(as.numeric(sub(".*_","",
                                             colnames(SWC_df)[c(-1,-2)])),
                              origin=paste0('2010-12-31')))
    colnames(SWC_df)=c("x","y",time)
    # selection of max/min per year
    SWC_df=SWC_df[,loc_ID:=.I]
    SWC_df=melt.data.table(SWC_df,id.vars=c("loc_ID","x","y"))
    SWC_df[,year:=gsub("[-].*", "\\1",variable)]
    SWC_df[,month:=gsub(".*[-]([^.]+)[-].*", "\\1",variable)]
    SWC_df[,day:=gsub(".*[-]", "\\1",variable)]
    SWC_df=SWC_df[,.(Min=min(value),Max=max(value)),by="loc_ID,x,y,year"]
    SWC_df=SWC_df %>% 
      group_by(loc_ID,x,y) %>% 
      summarise(SWC_min=mean(head(sort(Min),5)),
                SWC_max=mean(tail(sort(Max),5))) %>%
      ungroup()
    return(SWC_df)
  })
  return(SWC_mx)
}




#' Compute minimum soil potential per horizons
#' 
#' @description function applies VG equation to compute minimum soil 
#' potential per horizons using parameters extracted from 3D hydrosoil
#' @param SWCdh list of min/max per horizon
#' @param horizons number of horizons for SWC
#' @param dir.data 
#' @param dir.file 
#' @return Psi_min dataframe, that contains values of Psi_min computed with the

soil_potential <- function(SWCdh,
                           horizons=3,
                           dir.data="data",
                           dir.file="EU_SoilHydroGrids_1km"){
  # necessary functions 
  kseriesum <- function(k1,k2) {return(1/(1/k1+1/k2))}
  compute.B_GC <- function (La, Lv, rootRadius) {
    b <- 1 / sqrt(pi * Lv)
    return(La * 2 * pi / (log(b/rootRadius)))
  }
  
    # compute root characteristics
  betaRootProfile =0.97 
  depth <- c(0.05,0.15,0.3,0.6,1)
  SoilVolume <- depth 
  rootDistribution <-numeric(length(depth)) #Three soil layers
  rootDistribution[1] = 1-betaRootProfile^(depth[1]*100) # conversion of depth to cm 
  rootDistribution[2] = (1-betaRootProfile^(depth[2]*100))-rootDistribution[1]
  rootDistribution[3] = (1-betaRootProfile^(depth[3]*100))-(rootDistribution[1]+rootDistribution[2])
  rootDistribution[4] = (1-betaRootProfile^(depth[4]*100))-(rootDistribution[1]+rootDistribution[2]+rootDistribution[3])
  rootDistribution[5] = 1-(rootDistribution[1]+rootDistribution[2]+rootDistribution[3]+rootDistribution[4])
  k_RootToStem   = 1 * rootDistribution #Kroot to stem for eachh root
  
  fRootToLeaf=1
  LAImax = 5
  RAI = LAImax*fRootToLeaf
  rootRadius=0.0004
  La = RAI*rootDistribution / (2*pi*rootRadius)
  Lv = La/(SoilVolume)
  compute.B_GC(La,Lv,rootRadius)

  #Compute psi and ksoil and ksoiltostem
  
  pars.files=list.files(file.path(dir.data,dir.file))
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  hor_correspondance=matrix(data=c(1,2,2,3,3,4,4,1,2,3,4,5,6,7),ncol=2)
  for (i in 1:5){
    assign(paste0("psi_dh",i),
           c(rast(SWCdh[[hor_correspondance[i,1]]][,2:5]),
             resample(rast(file.path(dir.data,dir.file,mrc.files[endsWith(mrc.files,paste0(i,".tif"))])),
                      rast(SWCdh[[hor_correspondance[i,1]]][,2:5]),
                      method="near"),
             resample(rast(file.path(dir.data,dir.file,hcc.files[endsWith(hcc.files,paste0(i,".tif"))])),
                      rast(SWCdh[[hor_correspondance[i,1]]][,2:5]),
                      method="near")
           )
           )
  }
  
  list=mget(ls()[grepl("psi_dh",ls())])
  psi_min=lapply(seq_along(list),function(i){
    psi_df=as.data.frame(list[[i]],xy=TRUE)
    REW_mrc=(psi_df[,3]-psi_df[,8]*10^(-4))/((psi_df[,9]-psi_df[,8])*10^(-4))
    REW_hcc=(psi_df[,3]-psi_df[,15]*10^(-4))/((psi_df[,16]-psi_df[,15])*10^(-4))
    psi_df$psi_min=-(((1/REW_mrc)
                      ^(1/(psi_df[,6]*10^(-4)))
                      -1)
                     ^(1/(psi_df[,7]*10^(-4)))
                     *(1/(psi_df[,5]*10^(-4)))
                     *9.78*10^(-2))
    psi_df$ksoil=(psi_df[,11]
                  *(REW_hcc^(psi_df[,12]*10^(-4)))
                  *(1-(1-REW_hcc^(1/(psi_df[,13]*10^(-4))))
                    ^(psi_df[,13]*10^(-4)))^2)
    psi_df$ksoilGC=(psi_df[,11]
                  *(REW_hcc^(psi_df[,12]*10^(-4)))
                  *(1-(1-REW_hcc^(1/(psi_df[,13]*10^(-4))))
                    ^(psi_df[,13]*10^(-4)))^2 *1000 *compute.B_GC(La,Lv,rootRadius)[i])
    psi_df$k_soiltostem=kseriesum(psi_df[,19],k_RootToStem[i])
    psi_df$psi_w=psi_df[,19]*psi_df[,17]
    return(psi_df)
  }
  )
  rm(list)
  return(psi_min)
}


plot_psi=psi_min[[1]][,1:2]
for(i in 1:length(psi_min)){
  plot_psi=cbind(plot_psi,psi_min[[i]][,c(3,4,17)])
  names(plot_psi)=c(names(plot_psi)[1:(2+3*(i-1))],paste0(names(psi_min[[i]][,c(3,4,17)]),i))
}
plot_psi %>% 
  mutate(across(.cols=matches(c("psi")),
                ~ cut(.,
                      breaks=c(-Inf, -20000, -10000, -5000, -3000, -1000, -500, -300, -200, -100, Inf),
                      labels=c("<-20MPa", "-20<Psi<-10MPa", "-10<Psi<-5MPa", "-5<Psi<-3MPa","-3<Psi<-1MPa",
                               "-1<Psi<-0.5MPa","-0.5/-0.3","-0.3/-0.2","-0.2/-0.1","-0.1/0"))))  %>% 
  select(matches(c("^x$","^y$","psi"))) %>% 
  relocate(y,.after=x) %>% 
  pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=value)) +
  theme_bw() +
  theme(axis.title=element_blank(),
        legend.key.size = unit(0.5,"cm"))+
  labs(fill="Potentiel (MPa)")+
  facet_wrap(~name)+
  scale_fill_brewer(palette="RdYlBu")+
  coord_quickmap()


sum_psi=numeric(dim(psi_min[[i]])[1])
sum_k=numeric(dim(psi_min[[i]])[1])
for (i in 1:length(psi_min)){
  sum_psi=sum_psi+psi_min[[i]][,21]
  sum_k=sum_k+psi_min[[i]][,19]
}
psi=sum_psi/sum_k 

plot_psi=cbind(plot_psi,psi)
plot_psi


