#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analyses.R  
#' @description R script containing all functions relative to data
#               analyses
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Presence/absence dataset ####
#' @description Download pres/abs dataset
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

get.mauri <- function(dir.occ="data/EUForestsMauri/EUForestspecies.csv",
                      df.traits){
  #load species
  #df.traits=read.csv("output/df_trait_filtered.csv")
  
  df.Euforest=read.table("data/EUForestsMauri/EUForestspecies.csv",
                         header=TRUE,
                         sep=",",
                         dec=".")%>% 
    #remove duplicated location_tree 
    filter(!duplicated(paste(X,Y, SPECIES.NAME))) %>%
    #create plot id
    mutate(Pres = rowSums(cbind(NFI==1,FF==1,BS==1)),
           Node = unclass(factor(paste(X, Y))))
  #filter targetted species
  # filter(SPECIES.NAME %in% as.character(df.species$species.binomial))
  
  db.cont=df.Euforest %>% 
    group_by(Node,SPECIES.NAME) %>%
    summarise(n=n())%>%
    spread(SPECIES.NAME, n) %>% 
    ungroup() %>%
    pivot_longer(colnames(.)[-1],
                 names_to = "species.binomial",
                 values_to = "presence") %>% 
    filter(species.binomial %in% as.character(df.traits$species.binomial)) %>% 
    left_join(df.Euforest[!duplicated(df.Euforest$Node), c("X","Y","COUNTRY","EEO","Node")],
              by = "Node") %>% 
    left_join(df.traits[,c("species.binomial","sp.ind","p.trait","PX.mu","PX.sd","LT50.mean","LT50.sd","data.quality")],#df.traits[,c("species.binomial","sp.ind","MATmean","MAPmean","P50tot","LT50.mean","LT50.sd")],
              by="species.binomial" )
  #Change coordinates
  coordinates(db.cont) <- c("X",  "Y")
  proj4string(db.cont) <- CRS("+init=epsg:3035")
  db.cont <- as.data.frame(spTransform(db.cont, CRS("+proj=longlat +datum=WGS84")))
  
  
  
  ## Load climate
  # fdg
  rast.fdg=rast("data/eea_2000-2016/SOS_2000_2016_DOY.BIL")
  rast.fdg.mean <- mean(rast.fdg,na.rm=TRUE)
  rast.fdg.mean <- terra::project(rast.fdg.mean,"epsg:4326")
  rast.fdg.mean[rast.fdg.mean<0] <- 30
  names(rast.fdg.mean)="fdg"
  
  rm(rast.fdg)
  # function to get sgdd and wai
  get_waisgdd <- function(dir.chelsa="data/CHELSA/",
                          df.loc){
    annual_temp=extract(rast("data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif"),
                        df.loc)[2]
    pr_chelsa=extract(rast("data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif"),
                      df.loc)[2]
    pet_chelsa=extract(rast("data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif"),
                       df.loc)[2]
    sgdd_chelsa=extract(rast("data/CHELSA/CHELSA_gdd5_1981-2010_V.2.1.tif"),
                        df.loc)[2]
    df.loc=cbind(df.loc,
                 temp.mean=annual_temp,
                 pr=pr_chelsa,
                 pet=pet_chelsa,
                 sgdd=sgdd_chelsa)
    colnames(df.loc)=c("x","y","temp.mean","pr","pet","sgdd")
    df.loc <- df.loc %>% 
      mutate(wai=(pr-12*pet)/pet)
    return(df.loc)
  }
  
  # get safety margins
  psi_min=rast(fread("output/psihorday_real.csv")[,c("x","y","psi")],crs="epsg:4326")
  frost.index.spring=rast(fread("output/budburst_tquant.csv")[,-1],crs="epsg:4326")
  names(frost.index.spring)="frost.spring"
  frost.index.winter=min(rast("data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc"),na.rm=FALSE)
  frost.index.winter=classify(frost.index.winter, cbind(6553500, NA)) #set as NA default value
  names(frost.index.winter)="frost.winter"
  
  
  
  #Extract climate
  for (rast in list(psi_min,frost.index.spring,frost.index.winter,rast.fdg.mean)){
    print(rast)
    db.cont <- cbind(db.cont,
                     extract(rast,
                             data.frame(x=db.cont$X,y=db.cont$Y))[-1])
  }
  db.cont$presence[is.na(db.cont$presence)] <- 0
  date.dehardening=0
  df.fdg.sp <- db.cont %>%
    filter(presence==1) %>% 
    select(species.binomial,sp.ind,LT50.mean,fdg) %>% 
    group_by(species.binomial,sp.ind,LT50.mean) %>% 
    summarise(fdg.sp.mean=mean(fdg,na.rm=TRUE),
              fdg.sp.sd=sd(fdg,na.rm=TRUE)) %>% 
    ungroup() %>% 
    mutate(LT50.sp.spring=-5+(30/2)*((5+LT50.mean)/(fdg.sp.mean-date.dehardening)))
  
  db.cont <- db.cont %>% 
    left_join(df.fdg.sp[,c("species.binomial","LT50.sp.spring")],by="species.binomial") %>% 
    mutate(psi_old=psi,
           psi=case_when(psi<(-20000)~(-20000),
                         TRUE~psi)) %>% 
    mutate(HSM=psi-PX.mu*1000,
           FSMwinter=frost.winter-LT50.mean,
           FSMspring=frost.spring-(LT50.sp.spring),
           species.binomial=as.factor(species.binomial))
  db.cont$HSM[is.nan(db.cont$HSM)] <- NA
  db.cont$FSMwinter[is.nan(db.cont$FSMwinter)] <- NA
  db.cont$FSMspring[is.nan(db.cont$FSMspring)] <- NA
  
  # Extract climatic variables
  db.clim <- cbind(
    db.cont,
    get_waisgdd(df.loc=data.frame(x=db.cont$X,
                                  y=db.cont$Y))[,3:7])
  write.csv(db.clim,"output/db_EuForest.csv")
  write.csv(df.fdg.sp,"output/LT50_spring.csv")
  
  return(db.clim)
}
