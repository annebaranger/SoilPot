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
                      dir.species="data/Species traits/trait.select.csv"){
  #load species
  df.traits=read.csv("output/df_trait_filtered.csv")
  
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



#' resampling absence data into ecoregions
#' @description Using maps of ecorgeions in Europe, identify areas without forest
#' to sample new absence data for nothern species
#' @param europe Europe boarders
#' @param dir.ecoregions directory of ecoregions shapefiles
#' @note ecoregions are downloaded from the wwf
#' @param db.clim pres/abs dataset with climatic variables
#' @param safety.margins
#' @param n_sample numbers of resampled pseudo-abs
#' @return a map of occ data overlayed with distribution polygon
#'

sample.absence <- function(europe,
                           dir.ecoregions="data/WWF/official",
                           db.clim,
                           safety.margins,
                           n_sample,
                           dir.species="data/Species traits/trait.select.csv"){
  #Load species selection
  df.species=read.table(file=dir.species,
                        header=TRUE,
                        sep=";",
                        dec=".") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2)),
           species.name=str_replace(species.binomial," ",".")) %>% 
    filter(!is.na(LT50))
  
  # Load Ecoregions to european scale
  europe=st_as_sf(europe)
  ecoregions=read_sf(dsn=dir.ecoregions,
                     layer="wwf_terr_ecos") 
  sf_use_s2(FALSE)
  ecoregions=st_crop(ecoregions,sf::st_bbox(europe))
  
  # Draw buffer around abs+pres data using only 20000 points (for speed matter)
  hull <- sample_n(db.clim,
                   20000)  %>%
    st_as_sf( coords = c( "X", "Y" ), crs = 4326 ) %>% 
    # summarise( geometry = st_combine( geometry ) ) %>% 
    st_buffer(dist=units::as_units(0.25, 'degree')) %>% 
    st_union()
  
  abs.sampling=hull %>% 
    st_difference(x=ecoregions[ecoregions$ECO_NAME=="Scandinavian Montane Birch forest and grasslands",
                               c("ECO_NAME","geometry")],
                  y=.) %>% 
    st_sample(n_sample,
              type="regular")
  
  
  abs.plots=data.frame(st_coordinates(abs.sampling)) %>% 
    crossing(df.species[,c(1,4,5)]) %>% 
    cbind(.,
          extract(rast(safety.margins$psimin,crs="epsg:4326")[["psi"]],
                  data.frame(x=.$X,y=.$Y))[-1],
          extract(rast(safety.margins$tmin,crs="epsg:4326")[["t.min"]],
                  data.frame(x=.$X,y=.$Y))[-1]) %>% 
    mutate(Node= max(db.clim$Node)+ unclass(factor(paste(X, Y))),
           fsm=t.min-LT50,
           hsm=psi-P50*1000,
           species.binomial=as.factor(species.binomial))
  abs.plots$presence <- 0
  abs.plots$fsm[is.nan(abs.plots$fsm)] <- NA
  abs.plots$hsm[is.nan(abs.plots$hsm)] <- NA
  
  
  abs.plots <- cbind(
    abs.plots,
    get_waisgdd(df.loc=data.frame(x=abs.plots$X,
                                  y=abs.plots$Y))[,3:6])
  
  db.clim=db.clim %>% 
    bind_rows(abs.plots)
  
  return(db.clim)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Generalized linear mixte model ####
#' @description 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Glmm IFN France extent
#' 
#' @description Fit a GLMM using France IFN presnece/absence data
#' @param db.fni Fni dataset
#' @param safety.margins Frost and drought safety margins computed over a 
#' selection of species
#' @return a glmm fitted
#'

sm.glmm.ifn <- function(db.fni,
                    safety.margins,
                    dir.data="data/Species traits/",
                    dir.file="trait.select.csv"){
  df.species=read.table(file=paste0(dir.data,dir.file),
                        header=TRUE,
                        sep=";",
                        dec=".") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2))) %>% 
    filter(!is.na(LT50))
  
  # Extract values per sampled stand
  stand.sample=sample(unique(db.fni$stand$IDP), 
                      dim(db.fni$stand)[1]/5, 
                      replace = FALSE,
                      prob = NULL)
  db.stand.sample=db.fni$stand[IDP %in% stand.sample]
  db.stand.sample=unique(db.stand.sample[,c("IDP","XL","YL")] )
  db.stand.sample=st_as_sf(db.stand.sample, coords = c("XL", "YL"), crs = 2154)
  db.stand.sample=st_transform(db.stand.sample,crs="epsg:4326")

  db.stand.sample.hsm=cbind(db.stand.sample,
                        extract(rast(safety.margins$psimin[,-3],
                                     crs="epsg:4326"),
                                vect(db.stand.sample))) %>%
    dplyr::select(-ID) %>% 
    pivot_longer(cols=setdiff(colnames(.),c("IDP","geometry")),values_to = "hsm") %>% 
    filter(name %in% df.species$sp.ind) %>% 
    left_join(df.species[,c("fni.code","sp.ind")],by=c("name"="sp.ind")) 
    
  
  db.stand.sample.fsm=cbind(db.stand.sample,
                             extract(rast(safety.margins$tmin[,-3],
                                          crs="epsg:4326"),
                                     vect(db.stand.sample))) %>%
    dplyr::select(-ID) %>% 
    pivot_longer(cols=setdiff(colnames(.),c("IDP","geometry")),values_to = "fsm") %>% 
    filter(name %in% df.species$sp.ind) %>% 
    left_join(df.species[,c("fni.code","sp.ind")],by=c("name"="sp.ind"))
  
  db.stand.sample=cbind(db.stand.sample.hsm,
                        fsm=db.stand.sample.fsm$fsm)
  rm(db.stand.sample.fsm,
     db.stand.sample.hsm)
  # db.stand.sample=cbind(db.stand.sample,
  #                       extract(rast(safety.margins$tmin,
  #                                    crs="epsg:4326"),
  #                               vect(db.stand.sample)))
   
  #  Look presence of trees per stand
  db.tree.sample=db.fni$tree[IDP %in% stand.sample]
  db.tree.sample=unique(db.tree.sample[-which(db.tree.sample$ESPAR == ""), .(IDP,ESPAR)])
  db.tree.sample[, presence:=1]
  db.tree.sample=db.tree.sample[ESPAR %in% df.species$fni.code]
  
  # Build datatable with only presence absence of each species in each stand
  db.stand.sample=db.stand.sample %>%
    filter(hsm>(-25000)) %>% 
    rename(`ESPAR`="fni.code") %>% 
    left_join(db.tree.sample,by=c("IDP","ESPAR")) %>% 
    mutate(ESPAR=as.factor(ESPAR),
           IDP=as.factor(IDP))
  db.stand.sample[is.na(db.stand.sample)] <- 0
  
  db.stand.sample=db.stand.sample %>%
    mutate(hsm=scale(hsm,scale=TRUE,center=FALSE),
           fsm=scale(fsm,scale=TRUE,center=FALSE))
  # Fit glmm
  sm.glmm=glmer(presence ~ fsm + hsm + (1|ESPAR), 
                data = db.stand.sample,
                family=binomial(link = "logit"),
                nAGQ=20
                )
}

#' Glmm IFN Europe
#' @description Fit a GLMM using Mauri presence/absence data
#' @param safety.margins Frost and drought safety margins computed over a 
#' selection of species
#' @return a glmm fitted
#'

sm.glmm.mauri <- function(db.clim,
                          safety.margins,
                          dir.occ="data/EUForestsMauri/EUForestspecies.csv",
                          dir.species="data/Species traits/trait.select.csv"){
  db.clim=fread("output/db_EuForest.csv")
    # rename_with(.cols=everything(),
    #             tolower) %>% 
    # filter(!is.na(hsm)) %>% 
    # filter(!is.na(fsmwinter)) %>% 
    # filter(!is.na(fsmspring)) %>% 
    # filter(!is.na(sgdd)) %>% 
    # filter(!is.na(wai)) 
  #Load species selection
  df.traits=read.csv("output/df_trait_filtered.csv") 
  
  # Pres/abs tidy dataframe :
  ## filter extreme HSM
  ## scale & center safety margins per species -> relevant?
  hsm.q01=quantile(db.clim$HSM,0.01,na.rm=TRUE)
  db.pres <- db.clim %>%
    rename_with(.cols=everything(),
                tolower) %>% 
    filter(species.binomial%in%c("Pinus sylvestris","Fagus sylvatica","Fraxinus excelsior","Quercus ilex")) %>% 
    filter(hsm>hsm.q01) %>% 
    filter(!is.na(hsm)) %>% 
    filter(!is.na(fsmwinter)) %>% 
    filter(!is.na(fsmspring)) %>% 
    filter(!is.na(sgdd)) %>% 
    filter(!is.na(wai)) 
  hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
  hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
  fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
  fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)
  # group_by(species.binomial) %>% 
  db.pres <- db.pres %>% 
    mutate(hsm=scale(hsm,scale=TRUE,center=TRUE),
           fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE),
           pet=scale(pet,scale=TRUE,center=TRUE),
           sgdd=scale(sgdd,scale=TRUE,center=TRUE)) 
  #   ungroup() %>% 

  
  # How centering and scaling predictors
  # db.pres %>%
  #   filter(hsm>-10000) %>%
  #   sample_frac(0.3) %>%
  #  # group_by(species.binomial) %>%
  #   # mutate(hsm=scale(hsm,scale=TRUE,center=TRUE),
  #   #        fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)) %>%
  #   #ungroup() %>%
  #   ggplot()+
  #   geom_density(aes(hsm,color=species.binomial))+
  #   geom_density(aes(hsm),color="red")
  
  ## Fit glmer with predictors centered scaled overall 
  glm.winter=glmer(presence ~  hsm + fsmwinter + (1|species.binomial), 
                data = db.pres,
                family=binomial(link = "logit")) 
  #save(glm.winter,file="output/model/glm.winter.RData")
  glm.spring=glmer(presence ~  hsm + fsmspring + (1|species.binomial), 
                   data = db.pres,
                   family=binomial(link = "logit"))  
  #save(glm.spring,file="output/model/glm.spring.RData")

  ## Fit glm per species
  summary.sp <- db.pres %>% 
    filter(presence==1) %>% 
    group_by(species.binomial) %>% 
    summarise(pr_mean=mean(pr,na.rm=TRUE),
              temp_mean=mean(temp.mean,na.rm=TRUE),
              sgdd_mean=mean(sgdd,na.rm=TRUE),
              pet_mean=mean(pet,na.rm=TRUE)) %>% 
    mutate(hsm.coef=NA,
           fsm.coef=NA)
  fsm.param=c()
  hsm.param=c()
  for (species in unique(db.pres$species.binomial)){
    #hsm.q01=quantile(db.pres$HSM,0.01,na.rm=TRUE)
    db.pres.sp<- db.pres %>%
      filter(species.binomial==species) 
      # filter(HSM>hsm.q01) %>% 
      # mutate(hsm=scale(hsm,scale=TRUE,center=TRUE),
      #        fsm=scale(fsm,scale=TRUE,center=TRUE),
      #        pet=scale(pet,scale=TRUE,center=TRUE),
      #        sgdd=scale(sgdd,scale=TRUE,center=TRUE))
    glm.winter.sp=glm(presence ~  hsm + fsmwinter, 
                    data = db.pres.sp,
                    family=binomial(link = "logit"))
    hsm.param=c(hsm.param,glm.winter.sp$coefficients[["hsm"]][[1]])
    fsm.param=c(fsm.param,glm.winter.sp$coefficients[["fsmwinter"]][[1]])
    # glm.spring.sp=glm(presence ~  hsm + fsmspring, 
    #                  data = db.pres.sp,
    #                  family=binomial(link = "logit"))  
    print(species)
    print(summary(glm.winter.sp))
    # print(summary(glm.spring.sp))
  }
  summary.sp$hsm.coef <- hsm.param
  summary.sp$fsm.coef <- fsm.param
  
  
  summary.sp %>% 
    ggplot(aes(temp_mean,hsm.coef))+
    geom_point()+
    theme_bw()+
    xlab("Mean temperature of species distribution")+
    ylab("Coefficient of HSM")
  
  summary.sp %>% 
    ggplot(aes(temp_mean,fsm.coef))+
    geom_point()+
    theme_bw()+
    xlab("Mean temperature of species distribution")+
    ylab("Coefficient of FSM")
  
  
  ## sous sampling dans les données d'absences
  db.pres.pres <- db.pres %>% 
    filter(presence==1) 
  db.pres.abs <- db.pres %>% 
    filter(presence==0) %>% 
    sample_frac(0.4)
  db.pres.sample <- rbind(db.pres.pres,
                          db.pres.abs)
  rm(db.pres.abs,db.pres.pres)
  glm.winter.sample=glmer(presence ~  hsm + fsmwinter + (1|species.binomial), 
                   data = db.pres,
                   family=binomial(link = "logit")) 
  #save(glm.winter.sample,file="output/model/glm.winter.sample.RData")
  glm.spring.sample=glmer(presence ~  hsm + fsmspring + (1|species.binomial), 
                   data = db.pres,
                   family=binomial(link = "logit"))  
  #save(glm.spring.sample,file="output/model/glm.spring.sample.RData")


  ## conditionnal model on temperature zone
  species.zone <- db.pres %>% 
    filter(presence==1) %>% 
    mutate(species.binomial=as.factor(species.binomial)) %>% 
    group_by(species.binomial,presence) %>% 
    summarise(temp.quart.sup=quantile(temp.mean,probs=0.6)[[1]],
              temp.quart.inf=quantile(temp.mean,probs=0.4)[[1]]) %>% 
    ungroup() %>% 
    select(-presence)
  
  db.pres.south <- db.pres %>% 
    left_join(species.zone, by="species.binomial") %>% 
    filter(temp.mean>temp.quart.inf)
  db.pres.north <- db.pres %>% 
    left_join(species.zone, by="species.binomial") %>% 
    filter(temp.mean<temp.quart.sup) %>% 
    mutate(fsm=fsmwinter + fsmspring)
  
  glm.south=glmer(presence ~  hsm + (1|species.binomial), 
                   data = db.pres.south,
                   family=binomial(link = "logit")) 
  #save(glm.south,file="output/model/glm.south.RData")
  glm.north=glmer(presence ~ fsmwinter + (1|species.binomial), 
                   data = db.pres.north,
                   family=binomial(link = "logit"))  
  #save(glm.north,file="output/model/glm.north.RData")
  
  # test for nothern part
  glm.north.winter=glmer(presence ~ fsmwinter + (1|species.binomial), 
                         data = db.pres.north,
                         family=binomial(link = "logit"))
  glm.north.spring=glmer(presence ~ fsmspring + (1|species.binomial), 
                         data = db.pres.north,
                         family=binomial(link = "logit"))
  glm.north.wisp=glmer(presence ~ fsm + (1|species.binomial), 
                       data = db.pres.north,
                       family=binomial(link = "logit"))
  
  ## segmented regression
  hsm.bkpt=-hsm.mu/hsm.sd
  fsm.bkpt=-fsm.mu/fsm.sd
  glm.winter=glmmTMB(presence ~  hsm + fsmwinter + (1|species.binomial), 
                     data = db.pres,
                     family=binomial(link = "logit"))
  glm.segmented=segmented(obj=glm.winter,seg.Z=~hsm+fsmwinter) #,psi=list(hsm=hsm.bkpt,fsmwinter=fsm.bkpt)
  
  ## safety margins per species
  df.shadetol=read.csv2("data/Species traits/data_Niinemets&Valladares_2006.csv")
  db.clim %>% 
    filter(presence==1) %>% 
    group_by(species.binomial) %>% 
    mutate(n=n()) %>% 
    filter(n>1000) %>% 
    ggplot(aes(hsm,species.binomial))+
    geom_boxplot()+
    scale_x_continuous(trans = ssqrt_trans)+
    geom_vline(xintercept=0,color="red")

  #plot percentage of presence that fall into a compatible zone of the species regarding safety margins
db.clim %>% 
    mutate(compatibility=case_when(hsm>0&fsmwinter>0~1,
                                   TRUE~0)) %>%
  group_by(species.binomial) %>% 
  filter(presence==1) %>% 
  mutate(n=n()) %>% 
  filter(n>1000) %>% 
  ungroup() %>% 
  group_by(species.binomial,n,compatibility) %>% 
  summarise(n.cat=n()) %>% 
  mutate(occupancy=n.cat/n,
         order=case_when(compatibility==0~1-occupancy,
                         TRUE~occupancy)) %>% 
  # filter(compatibility==1) %>% 
  # arrange(desc(occupancy))
  ggplot(aes(n.cat,reorder(species.binomial,order),
             fill=as.factor(compatibility)))+
  geom_bar(position = "fill",stat = "identity")

}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Exploratory Plots ####
#' @description 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Presence data comparison with predicted niche
#' @description Overlay species distribution
#' @param dir.occ Mauri occ data
#' @param dir.species list of species of interest
#' @return a map of occ data overlayed with distribution polygon
#'

overlay.abs <- function(dir.species="data/Species traits/trait.select.csv",
                        db.pres,
                        europe){
  #load europe extent in sf format
  europe=st_as_sf(europe)
  
  #load species list and files
  df.species=read.table(file=dir.species,
                        header=TRUE,
                        sep=";",
                        dec=".") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2)),
           species.name=str_replace(species.binomial," ",".")) %>% 
    filter(!is.na(LT50))
  
  #plot presence data and euforgen distribution
  for (i in 1:dim(df.species)[1]){
    tryCatch(
      {
        spdistrib=read_sf(dsn=file.path("data",
                                  "chorological_maps_dataset",
                                  df.species$species.binomial[i],
                                  "shapefiles"),
                          layer=df.species$file[i])
        if(is.na(df.species$extra.file[i])==FALSE){
          spdistrib_sup=read_sf(dsn=file.path("data",
                                              "chorological_maps_dataset",
                                              df.species$species.binomial[i],
                                              "shapefiles"),
                                layer=df.species$extra.file[i])
          spdistrib=st_union(rbind(spdistrib,spdistrib_sup))
          }
        spdistrib=st_crop(spdistrib,sf::st_bbox(europe))
        
        ## Build graphs
        presence= db.pres %>%
          filter(species.binomial==df.species$species.binomial[i]) %>%
          filter(presence==1) %>% 
          ggplot()+
          geom_sf(data=spdistrib,
                  fill=alpha("grey",0.6),
                  lwd = 0)+
          geom_sf(data=europe,
                  fill=NA)+
          geom_point(aes(X,Y),size=0.5)+
          theme_bw()
        ggsave(presence,
               filename = paste0(df.species$sp.ind[i],"presence.png"),
               path="figs/",
               device="png",
               scale=2)
        absence=db.pres %>%
          filter(species.binomial==df.species$species.binomial[i]) %>%
          filter(presence!=1) %>% 
          ggplot()+
          geom_sf(data=spdistrib,
                  fill=alpha("grey",0.6),
                  lwd = 0)+
          geom_sf(data=europe,
                  fill=NA)+
          geom_point(aes(X,Y),
                     size=0.05)+
          theme_bw()
        ggsave(absence,
               filename = paste0(df.species$sp.ind[i],"absence.png"),
               path="figs/",
               device="png",
               scale=2)
      },
      error=function(e){print("File not loaded")}
    )
  }
}


#' Random plots for exploratory analysis
#' @description 1 Ecoregions & absence of a species
#' @param dir.occ Mauri occ data
#' @param dir.species list of species of interest
#' @return a map of occ data overlayed with distribution polygon
#'

explots <- function(db.clim,
                    europe,
                    dir.ecoregions="data/WWF/official"){
  # Load Ecoregions to european scale
  europe=st_as_sf(europe)
  ecoregions=read_sf(dsn=dir.ecoregions,
                     layer="wwf_terr_ecos") 
  sf_use_s2(FALSE)
  ecoregions=st_crop(ecoregions,sf::st_bbox(europe))
  
  # Plot of ecoregions + absence of a particular species
  db.clim %>%
    filter(species.binomial=="Picea abies") %>%
    filter(presence!=1) %>% 
    ggplot()+
    geom_sf(data=ecoregions,aes(fill=ECO_NAME),
            show.legend = FALSE)+
    geom_point(aes(X,Y),
               size=0.05)+
    theme_bw()+
    coord_sf()
  
  # Plot buffer and newly sampled pseudo-abs over europe
  db.clim %>% 
    dplyr::select(Node,X,Y,t.min,psi) %>% 
    distinct() %>% 
    mutate(pseudo.abs=case_when(eval(max(db.clim$Node))-Node<10008~"Added", ## warning , change the formula according to number of pseudo-abs
                                TRUE~"Mauri")) %>%
    filter(pseudo.abs=="Added") %>% 
    ggplot()+
    geom_point(aes(x=X,y=Y),size=0.01)+
    geom_sf(data=ecoregions,fill=NA)#+
  # geom_sf(data=hull,fill="red")

  # Plot proportion of pseudo absence in each bin of t.min 
  db.clim %>% 
    dplyr::select(Node,X,Y,t.min,psi) %>% 
    distinct() %>% 
    mutate(pseudo.abs=case_when(eval(max(db.clim$Node))-Node<10008~"Added", ## warning , change the formula according to number of pseudo-abs
                                TRUE~"Mauri")) %>% 
    ggplot(aes(t.min,fill=pseudo.abs))+
    geom_histogram(position = "fill")
  
  
}


# specific case of abies alba
db.clim %>% 
  mutate(pseudo.abs=case_when(eval(max(db.clim$Node))-Node<10008~"Added", ## warning , change the formula according to number of pseudo-abs
                              TRUE~"Mauri")) %>%
  filter(species.binomial=="Abies alba") %>% 
  filter(presence!=1) %>% 
  filter(pseudo.abs=="Mauri") %>% 
  ggplot(aes(fsm))+
  geom_density()

db.pres %>% 
  mutate(pseudo.abs=case_when(eval(max(db.clim$Node))-Node<10008~"Added", ## warning , change the formula according to number of pseudo-abs
                              TRUE~"Mauri")) %>%
  filter(species.binomial=="Abies alba") %>% 
  ggplot(aes(hsm,presence))+
  geom_point()+
  geom_smooth()

## density plot of presence according to hsm/fsm
db.clim %>% 
  filter(hsm>hsm.q05) %>% 
  mutate(pseudo.abs=case_when(eval(max(db.clim$Node))-Node<10008~"Added", ## warning , change the formula according to number of pseudo-abs
                              TRUE~"Mauri")) %>%
  filter(species.binomial=="Abies alba") %>% 
  filter(presence!=0) %>% 
  ggplot(aes(fsm,hsm))+
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

#subset data in the alps
db.alps=
    st_filter(
      st_as_sf(db.clim,
               coords=c('X','Y'),
               crs="epsg:4326"),
      ecoregions[ecoregions$ECO_NAME=="Alps conifer and mixed forests",]) %>% 
    dplyr::mutate(X = sf::st_coordinates(.)[,1],
                  Y = sf::st_coordinates(.)[,2])
db.alps %>% 
  filter(hsm>hsm.q05) %>% 
  filter(species.binomial=="Picea abies") %>%
  ggplot(aes(fsm,hsm))+
  geom_point(aes(colour=presence)) +
  theme_bw()
#%%%%%%%%%%%%%%
#### explo ####
#%%%%%%%%%%%%%%


# Look for relation between presence and all explicative variables
db.pres %>% 
  select(species.binomial,presence,fsmwinter,hsm,pet,sgdd) %>% 
  pivot_longer(cols = colnames(.)[c(-1,-2)],
               names_to = "variables",
               values_to = "value") %>% 
  sample_frac(0.3) %>% 
  ggplot(aes(value,presence))+ 
  geom_point()+
  geom_smooth(method="gam")+
  facet_grid(variables~species.binomial,scales = "free_x")
db.pres %>% 
  ggplot(aes(hsm,presence))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~species.binomial)
db.pres %>% 
  mutate(presence=as.factor(presence)) %>% 
  ggplot(aes(presence,hsm))+
  geom_boxplot()

db.pres %>% 
  ggplot(aes(-psi))+geom_density()+scale_x_log10()

db.pres %>% 
  #filter(species.binomial=="Fagus sylvatica") %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(fsm,presence)) +
  geom_point()+
  geom_smooth()

db.pres %>% 
  filter(presence==1) %>% 
  ggplot(aes(hsm,wai))+
  geom_point()


# draw corplots
corrplot <- db.clim %>% 
  select(-P50,-NFI,-FF,-BS,-EEO,-LT50,-COUNTRY) %>% 
  filter(!is.na(psi)) %>% 
  filter(!is.na(t.min)) %>% 
  filter(species.binomial%in%c("Picea abies")) %>% 
  filter(psi>quantile(db.clim$psi,probs=0.05,na.rm=TRUE)) %>% 
  group_by(species.binomial) %>% 
  mutate(hsm=hsm/1000,
         fsm.class=cut(fsm,
                       breaks=10),
         hsm.class=cut(hsm,
                       breaks=10)) %>% 
  ungroup() %>% 
  group_by(species.binomial,fsm.class,hsm.class) %>% 
  mutate(ntot=n(), #number of obs per species in each class
         npres=sum(presence>0),
         prop=npres/ntot+0.0001) %>% 
  ungroup() %>% 
  select(species.binomial,fsm.class,hsm.class,prop,npres,ntot) %>% 
  distinct() %>% 
  ggplot(aes(x=fsm.class,y=hsm.class))+
  geom_tile(aes(fill=prop))+
  scale_fill_distiller(palette="YlGnBu",direction=1)+
  geom_point(aes(size=ntot))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.6))+
  xlab("FSM class (celsius)")+
  ylab("HSM class (MPa)")+
  labs(fill="Presence proportion",
       size="Number of observations")+
  facet_wrap(~species.binomial,scales="free")

#%%%%%%%%
#### end
#%%%%%%%%
for (species in df.species$species.binomial){
  print(species)
  threshold.cold=quantile(db.pres[c(db.pres$species.binomial==species & db.pres$presence==1),]$t.min)[4]
  threshold.hot=quantile(db.pres[c(db.pres$species.binomial==species & db.pres$presence==1),]$t.min)[2]
  glm.cold=glm(formula = presence ~ hsm + fsm,
                     data=db.pres %>% 
                       filter(t.min<threshold.cold) %>% 
                       filter(species.binomial==species),
                     family=binomial(link = "logit"))
  glm.hot=glm(formula = presence ~ hsm + fsm,
               data=db.pres %>% 
                 filter(t.min>threshold.hot) %>% 
                 filter(species.binomial==species),
               family=binomial(link = "logit"))
  print("cold")
  print(summary(glm.cold))
  print("hot")
  print(summary(glm.hot))
  
}
threshold.cold=quantile(db.pres[c(db.pres$species.binomial=="Pinus sylvestris"& db.pres$presence==1),]$t.min)[2]
threshold.hot=quantile(db.pres[c(db.pres$species.binomial=="Abies alba"& db.pres$presence==1),]$t.min)[2]
glm.abies.cold=glm(formula = presence ~ fsm,
                  data=db.pres %>% 
                    filter(t.min<threshold.cold) %>% 
                    filter(species.binomial=="Pinus sylvestris"),
                  family=binomial(link = "logit"))

summary(glm.abies.cold)


#%%%%%%%%%%%%%
#### explo glm
#%%%%%%%%%%%%%
# step AIC using db.pres
fullmod=glmer(presence ~  hsm + fsm + sgdd + wai + (1|species.binomial), 
                data = db.pres,
                family=binomial(link = "logit"),
              na.action = "na.omit")
 nullmod=glmer(presence ~  (1|species.binomial), 
              data = db.pres,
              family=binomial(link = "logit"))
forwards=step(nullmod,
              scope=list(lower=formula(nullmod),
                         upper=formula(fullmod)),
              direction="forward",
              )
dredge(fullmod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Recompute vdp ####
#' @description 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compute.vpd <- function(europe,
                        dir.file="data/CHELSA/CHELSA_vpd_max_1981-2010_V.2.1.tif"){
  europe=st_as_sf(europe)
  vpd.max=crop(rast(dir.file),
               vect(europe),
               mask=TRUE)
  names(vpd.max)="vpd.max"
  psi.contrib=-5*(vpd.max/101325)
  vpd.psi <- as.data.frame(vpd.max,xy=TRUE) %>% 
    mutate(psi.contrib=-5*(vpd.max/101325))
  vpd.psi %>% 
    ggplot()+
    geom_tile(aes(x=x,y=y,fill=psi.contrib))+
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") + 
    scale_fill_gradientn(name = "potential min",colours = viridis(80),na.value="transparent")+
    coord_quickmap()
  
  return(vpd.psi)
}


