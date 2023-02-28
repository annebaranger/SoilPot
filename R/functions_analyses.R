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


#' Extract from Chelsea and compute WAI/SGDD
#' 
#' @description extract required variables from chelsea and computed sgdd and wai
#' on points entered as argument
#' @note Chelsea dwonloaded from website
#' @param dir.chelsa
#' @param df.loc data.frame of x and y location
#' @return df.loc + extracted variables + computed variables

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


#' Load Mauri pres/abs and format dataset
#' 
#' @description load Mauri dataset, format presence/absence and load safety margins
#' and others climatics variables associated to pres/abs points
#' @note Chelsea dwonloaded from website
#' @param dir.occ directory of mauri abs.pres dataset
#' @param df.traits dataframe of traits associated with eahc species
#' @param psi_min psimin computed over europe
#' @param frost.index frost index computed over europe
#' @return dataframe, where for each location, abs/pres of all selected species
#' is mentionned, and associated predictors

get.mauri <- function(dir.occ="data/EUForestsMauri/EUForestspecies.csv",
                      df.traits,
                      psi_min,
                      frost.index){
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

  # get safety margins
  psi_min=rast(fread(psi_min)[,c("x","y","psi")],crs="epsg:4326")
  frost.index.spring=rast(fread(frost.index)[,-1],crs="epsg:4326")
  names(frost.index.spring)="frost.spring"
  frost.index.winter=min(rast("data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc"),na.rm=FALSE)
  frost.index.winter=classify(frost.index.winter, cbind(6553500, NA)) #set as NA default value
  names(frost.index.winter)="frost.winter"
  
  
  
  #Extract safety margins on eahc abs.pres points
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
                                  y=db.cont$Y))[,3:7]) %>% 
    rename_with(.cols=everything(),
                tolower) %>% 
    filter(!is.na(hsm)) %>% 
    filter(!is.na(fsmwinter)) 
  write.csv(db.clim,"output/db_EuForest.csv")
  write.csv(df.fdg.sp,"output/LT50_spring.csv")
  
  return(db.clim)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Compute species characteristics ####
#' @description Compute index per species 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Compute prevalence
#' 
#' @description function that compute for each species in a dataframe its prevalence
#' as the percentage of presence in its presumed distribution
#' @note Distribution from EuForgen
#' @param df.traits df of traits for each species
#' @param db.clim database of pres/abs
#' @return df.loc + extracted variables + computed variables
#' 

get.prevalence <- function(df.traits,
                           db.clim){
  # create an empty column for prevalence
  df.traits$prevalence <- NA
  
  # for loop occuring along all species
  for (sp in unique(df.traits$species.binomial)){
    print(sp)
    # look for euforgen distribution if it exists for the targetted species and
    # compute the prevalelnce
    if (!is.na(df.traits[df.traits$species.binomial==sp,"file"][[1]])){
      # filter Mauri db with only points of sp
      db.pres <- db.clim %>% 
        filter(species.binomial==sp)
      #load euforgen distrib
      file.sp=df.traits[df.traits$species.binomial==sp,"file"][[1]]
      spdistrib=read_sf(dsn=file.path("data",
                                      "chorological_maps_dataset",
                                      sp,
                                      "shapefiles"),
                        layer=file.sp)
      # create spatial points with mauri db 
      db.pres.geo <- st_as_sf(db.pres,
                              coords=c("x","y"),
                              crs="epsg:4326") %>% 
        st_join(spdistrib, join = st_within,left=FALSE) %>% # select only points falling in euforgen distrib
        as.data.frame() 
      df.traits[df.traits$species.binomial==sp,"prevalence"] <-  sum(db.pres.geo$presence==1)/dim(db.pres.geo)[1]
    } else { # if euforgen distrib do not exist, prevalence default set to 0.1
      df.traits[df.traits$species.binomial==sp,"prevalence"] <- 0.1
    }
  }
  return(df.traits[,c("species.binomial","prevalence")])
}



#' Compute niche caracteristics
#' 
#' @description compute mean lat/long and continentality of species niche
#' @note Continentality from https://doi.org/10.1038/s41597-020-0464-0
#' @param df.traits df of traits for each species
#' @param db.clim database of pres/abs
#' @return df.loc + extracted variables + computed variables
#' 

get.niche <- function(df.traits,
                           db.clim){
  df.niche <- db.clim %>% 
    filter(presence==1) %>% 
    group_by(species.binomial) %>% 
    summarise(lat.mean=mean(y),
              lat.sd=sd(y),
              long.mean=mean(x),
              long.sd=sd(x),
              lat.q05=quantile(y,prob=0.05)[[1]],
              lat.q95=quantile(y,prob=0.95)[[1]],
              psi.q05=quantile(psi,prob=0.05)[[1]],
              psi.q95=quantile(psi,prob=0.95)[[1]],
              frost.q05=quantile(frost.winter,prob=0.05)[[1]],
              frost.q95=quantile(frost.winter,prob=0.95)[[1]]
    )
  
  
  rast.cont=as.data.frame(mean(rast("jci_year.nc")),xy=TRUE) %>% 
    mutate(z=x,
           x=y,
           y=z) %>% 
    select(-z) %>% 
    rast(crs="epsg:4326")
  df.traits$jci=NA
  for (sp in df.traits$species.binomial){
    print(sp)
    db.pres <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      filter(presence==1)
    db.pres <- cbind(db.pres,
                     jci=extract(rast.cont,db.pres[,c("x","y")])[["mean"]])
    jci=mean(db.pres$jci,na.rm=TRUE)
    df.traits[df.traits$species.binomial==sp,"jci"]=jci
  }
  
  df.overlap <- db.clim %>% 
    filter(presence==1) %>% 
    group_by(node) %>% 
    mutate(nb.sp=sum(presence==1),
           overlap.hsm=sum(hsm>0),
           overlap.fsm=sum(fsmwinter>0)) %>% 
    ungroup() %>% 
    group_by(species.binomial) %>% 
    summarise(overlap=mean(nb.sp),
              overlap.hsm=mean(overlap.hsm),
              overlap.fsm=mean(overlap.fsm))
  
  
  df.traits <- df.traits %>% 
    left_join(df.niche,by="species.binomial") %>% 
    left_join(df.overlap,by="species.binomial")
  return(df.traits[,c("species.binomial","lat.mean","lat.sd","long.mean","long.sd",
                      "lat.q05","lat.q95","psi.q05","psi.q95","frost.q05","frost.q95",
                      "jci","overlap","overlap.hsm","overlap.fsm")])
}

#' Get shade tol
#' 
#' @description compute mean lat/long and continentality of species niche
#' @note Continentality from https://doi.org/10.1038/s41597-020-0464-0
#' @param df.traits df of traits for each species
#' @return df.loc + extracted variables + computed variables
#' 

get.shadetol <- function(df.traits){
  df.shadetol=read.csv2("data/Species traits/data_Niinemets&Valladares_2006.csv")
  df.traits=df.traits %>% 
    left_join(df.shadetol,by=c("species.binomial"="Species"))
  return(df.traits[,c("species.binomial","shade_tolerance.mean","drought_tolerance.mean","waterlogging_tolerance.mean")]) 
}

#' Get traits and caract
#' 
#' @description compute niche traits and caract
#' @param df.traits df of traits for each species
#' @return df.loc + extracted variables + computed variables
#' 

get.species <- function(df.traits,
                        df.preval,
                        df.shadetol,
                        df.niche){
  df.species <- df.traits %>% 
    left_join(df.shadetol,by="species.binomial")%>% 
    left_join(df.niche,by="species.binomial")%>% 
    left_join(df.preval,by="species.binomial")
  write.csv(df.species,"output/df.species.csv")
  return(df.species) 
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Fit logistic regression ####
#' @description Compute index per species 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Fit model
#' 
#' @description function that fit 4 models for each species and save it
#' @param df.species df of traits for each species
#' @return no return, fit_mod6 is the folder with all fitted models
#' 

fit.logistic <- function(db.clim.file,
                         df.species){
  # read db.clim and filter some species
  db.clim=fread(db.clim.file) %>% 
      filter(!species.binomial %in% c("Juniperus virginiana",
                                      "Prunus persica",
                                      "Quercus palustris",
                                      "Cedrus libani",
                                      "Liriodendron tulipifera",
                                      "Pinus ponderosa",
                                      "Alnus cordata",
                                      "Acer negundo",
                                      "Prunus cerasifera",
                                      "Tsuga heterophylla",
                                      "Prunus serotina",
                                      "Juglans nigra",
                                      "Quercus rubra",
                                      "Pinus contorta",
                                      "Abies grandis",
                                      "Juglans regia",
                                      "Abies nordmanniana",
                                      "Pseudotsuga menziesii",
                                      "Robinia pseudoacacia"
      ))
  
  
  # set species list, as sp present in occ data
  species.list=unique(db.clim$species.binomial)
  
  # define a dataframe to 
  df.output=df.species %>% 
    filter(species.binomial %in% species.list) %>% 
    select(species.binomial,Group,p.trait,PX.mu,PX.sd,LT50.mean,LT50.sd,prevalence) %>% 
    mutate(K_int=NA,
           r_fsm=NA,
           r_hsm=NA,
           t_hsm=NA,
           t_fsm=NA,
           Rhat=NA,
           lp__=NA) %>% 
    crossing(mod=c("2sm","hsm","fsm","none"))

    for (sp in species.list){
    print(sp)
    db.pres <- db.clim %>%
      filter(psi>(-10000)) %>% #remove very low value of psi
      filter(species.binomial==sp) %>% 
      mutate(hsm=hsm/1000,
             fsm=fsmwinter) 
    
    # 2 var
    print("2sm")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      fsm=as.numeric(db.pres$fsm),
                      hsm=as.numeric(db.pres$hsm),
                      prior_K=df.species[df.species$species.binomial==sp,
                                        "prevalence"],
                      NULL)
    fit.2var <- stan(file = "glm_log_1sp_III.stan",
                     data=data.list,
                     iter=1000,
                     chains=3,
                     core=3,
                     include=FALSE,
                     pars=c("proba","K_vect"))
    save(fit.2var,file=paste0("fit_mod6/",sp,"_fitIII.RData"))
    summary=summary(fit.2var)$summary
    df.output[df.output$species.binomial==sp
              & df.output$mod=="2sm",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__")]=
      list(summary["K_int","mean"],
           summary["r_fsm","mean"],
           summary["r_hsm","mean"],
           summary["t_hsm","mean"],
           summary["t_fsm","mean"],
           max(summary[,"Rhat"]),
           summary["lp__","mean"]
      )
    
    # HSM
    print("hsm")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      hsm=as.numeric(db.pres$hsm),
                      prior_K=df.species[df.species$species.binomial==sp,
                                        "prevalence"],
                      NULL)
    fit.hsm <- stan(file = "glm_log_1sp_III_1varh.stan",
                    data=data.list,
                    iter=1000,
                    chains=3,
                    core=3,
                    include=FALSE,
                    pars=c("proba","K_vect"))
    save(fit.hsm,file=paste0("fit_mod6/",sp,"_fitIII_hsm.RData"))
    summary=summary(fit.hsm)$summary
    df.output[df.output$species.binomial==sp
              & df.output$mod=="hsm",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__")]=
      list(summary["K_int","mean"],
           0,
           summary["r_hsm","mean"],
           summary["t_hsm","mean"],
           0,
           max(summary[,"Rhat"]),
           summary["lp__","mean"]
      )
    
    # FSM
    print("fsm")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      fsm=as.numeric(db.pres$fsm),
                      prior_K=df.species[df.species$species.binomial==sp,
                                        "prevalence"],
                      NULL)
    fit.fsm <- stan(file = "glm_log_1sp_III_1var.stan",
                    data=data.list,
                    iter=1000,
                    chains=3,
                    core=3,
                    include=FALSE,
                    pars=c("proba","K_vect"))
    save(fit.fsm,file=paste0("fit_mod6/",sp,"_fitIII_fsm.RData"))
    summary=summary(fit.fsm)$summary
    df.output[df.output$species.binomial==sp
              & df.output$mod=="fsm",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__")]=
      list(summary["K_int","mean"],
           summary["r_fsm","mean"],
           0,
           0,
           summary["t_fsm","mean"],
           max(summary[,"Rhat"]),
           summary["lp__","mean"]
      )
    
    # None
    print("none")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      prior_K=df.species[df.species$species.binomial==sp,
                                        "prevalence"],
                      NULL)
    fit.none <- stan(file = "glm_log_1sp_III_0var.stan",
                     data=data.list,
                     iter=1000,
                     chains=3,
                     core=3,
                     include=FALSE,
                     pars=c("proba","K_vect"))
    save(fit.none,file=paste0("fit_mod6/",sp,"_fitIII_none.RData"))
    summary=summary(fit.none)$summary
    df.output[df.output$species.binomial==sp
              & df.output$mod=="none",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__")]=
      list(summary["K_int","mean"],
           0,
           0,
           0,
           0,
           max(summary[,"Rhat"]),
           summary["lp__","mean"]
      )
    
  }
}

#' Build output
#' 
#' @description load each model and extract relevant parameters
#' @param df.species df of traits for each species
#' @param output.path
#' @return df.output
#'

get.output <- function(db.clim.file,
                       df.species,
                       output.path){
  # read db.clim and filter some species
  db.clim=fread(db.clim.file) %>% 
    filter(!species.binomial %in% c("Juniperus virginiana",
                                    "Prunus persica",
                                    "Quercus palustris",
                                    "Cedrus libani",
                                    "Liriodendron tulipifera",
                                    "Pinus ponderosa",
                                    "Alnus cordata",
                                    "Acer negundo",
                                    "Prunus cerasifera",
                                    "Tsuga heterophylla",
                                    "Prunus serotina",
                                    "Juglans nigra",
                                    "Quercus rubra",
                                    "Pinus contorta",
                                    "Abies grandis",
                                    "Juglans regia",
                                    "Abies nordmanniana",
                                    "Pseudotsuga menziesii",
                                    "Robinia pseudoacacia"
    ))
  
  
  # set species list, as sp present in occ data
  species.list=unique(db.clim$species.binomial)
  
  # define a dataframe to 
  df.output=df.species %>% 
    filter(species.binomial %in% species.list) %>% 
    select(species.binomial,Group,p.trait,PX.mu,PX.sd,LT50.mean,LT50.sd,prevalence) %>% 
    mutate(K_int=NA,
           r_fsm=NA,
           r_hsm=NA,
           t_hsm=NA,
           t_fsm=NA,
           Rhat=NA,
           lp__=NA,
           divergence=NA) %>% 
    crossing(mod=c("2sm","hsm","fsm","none"))
  
  for (sp in species.list){
    print(sp)

    # 2 var
    print("2sm")
    load(file=paste0(output.path,sp,"_fitIII.RData"))
    summary=summary(fit.2var)$summary
    par_mod=get_sampler_params(fit.2var)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="2sm",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           summary["r_fsm","mean"],
           summary["r_hsm","mean"],
           summary["t_hsm","mean"],
           summary["t_fsm","mean"],
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
    
    
    # HSM
    print("hsm")
    load(file=paste0(output.path,sp,"_fitIII_hsm.RData"))
    summary=summary(fit.hsm)$summary
    par_mod=get_sampler_params(fit.hsm)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="hsm",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           0,
           summary["r_hsm","mean"],
           summary["t_hsm","mean"],
           0,
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
    
    # FSM
    print("fsm")
    load(file=paste0(output.path,sp,"_fitIII_fsm.RData"))
    summary=summary(fit.fsm)$summary
    par_mod=get_sampler_params(fit.fsm)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="fsm",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           summary["r_fsm","mean"],
           0,
           0,
           summary["t_fsm","mean"],
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
    
    # None
    print("none")
    load(file=paste0(output.path,sp,"_fitIII_none.RData"))
    summary=summary(fit.none)$summary
    par_mod=get_sampler_params(fit.none)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="none",
              c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","Rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           0,
           0,
           0,
           0,
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
    
  }
  
  #' Compute BIC
  df.output <- df.output %>% 
    mutate(nb.par=case_when(mod=="2sm"~5,
                            mod=="hsm"|mod=="fsm"~3,
                            mod=="none"~1),
           bic=2*nb.par-2*lp__)
  write.csv(df.output,file="output/df.output.csv")
  return(df.output)
}


#' Compute auc
#' 
#' @description compute auc for each model
#' @param df.species df of traits for each species
#' @param output.path
#' @return df.output
#'

compute.auc <- function(df.output,
                        db.clim.file){
  # read db.clim and filter some species
  db.clim=fread(db.clim.file) %>% 
    filter(!species.binomial %in% c("Juniperus virginiana",
                                    "Prunus persica",
                                    "Quercus palustris",
                                    "Cedrus libani",
                                    "Liriodendron tulipifera",
                                    "Pinus ponderosa",
                                    "Alnus cordata",
                                    "Acer negundo",
                                    "Prunus cerasifera",
                                    "Tsuga heterophylla",
                                    "Prunus serotina",
                                    "Juglans nigra",
                                    "Quercus rubra",
                                    "Pinus contorta",
                                    "Abies grandis",
                                    "Juglans regia",
                                    "Abies nordmanniana",
                                    "Pseudotsuga menziesii",
                                    "Robinia pseudoacacia"
    ))
  
  
  # set species list, as sp present in occ data
  species.list=unique(db.clim$species.binomial)
  
  df.output$auc=NA
  df.output$threshold=NA
  
  for (sp in species.list){
    print(sp)
    
    ### 2 var
    print("2sm")
    K_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="2sm",
                               "K_int"])
    r_fsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="2sm",
                               "r_fsm"])
    t_fsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="2sm",
                               "t_fsm"])
    r_hsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="2sm",
                               "r_hsm"])
    t_hsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="2sm",
                               "t_hsm"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=K_int/
               ((1+exp(-r_fsm*(fsmwinter-t_fsm)))*
                  (1+exp(-r_hsm*(hsm-t_hsm)))),
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="2sm",
              c("threshold","auc")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"])
    rm(db.pred,db.pred.auc,df.auc)
    
    # HSM
    print("hsm")
    K_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="hsm",
                               "K_int"])
    r_hsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="hsm",
                               "r_hsm"])
    t_hsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="hsm",
                               "t_hsm"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=K_int/
               ((1+exp(-r_hsm*(hsm-t_hsm)))),
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="hsm",
              c("threshold","auc")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"])
    rm(db.pred,db.pred.auc,df.auc)
    
    
    # FSM
    print("fsm")
    K_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="fsm",
                               "K_int"])
    r_fsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="fsm",
                               "r_fsm"])
    t_fsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="fsm",
                               "t_fsm"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=K_int/
               ((1+exp(-r_fsm*(fsmwinter-t_fsm)))),
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="fsm",
              c("threshold","auc")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"])
    rm(db.pred,db.pred.auc,df.auc)
    
    
    # None
    print("none")
    K_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="none",
                               "K_int"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=K_int,
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="none",
              c("threshold","auc")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"])
    rm(db.pred,db.pred.auc,df.auc)
    
    
  }
  
  return(df.output)
}


#' Model selection
#' 
#' @description compute auc for each model
#' @param df.output
#' @return df.output
#'

select.model <- function(df.output){
  #' Model selection
  df.mod.select <- df.output %>% 
    filter(Rhat<1.2) %>% 
    filter(divergence <0.1) %>% 
    group_by(species.binomial) %>% 
    slice(which.min(bic)) %>% 
    ungroup()
  
  #' Index of inflexion
  hsm.95=5.1#quantile(db.clim$hsm,prob=0.95)[[1]]/1000
  fsm.95=21.3#quantile(db.clim$fsmwinter,prob=0.95)[[1]]
  df.mod.select <- df.mod.select %>% 
    mutate(inflex_fsm=100*((16*exp(-r_fsm*(fsm.95-t_fsm)))/(2+2*exp(-r_fsm*(fsm.95-t_fsm)))^2),#+
           #(K_int*r_fsm*exp(-r_fsm*(fsm.05-t_fsm)))/(1+exp(-r_fsm*(fsm.05-t_fsm)))^2,
           #K_int*r_fsm/4,
           inflex_hsm=100*((16*exp(-r_hsm*(hsm.95-t_hsm)))/(2+2*exp(-r_hsm*(hsm.95-t_hsm)))^2)#+
           #(K_int*r_hsm*exp(-r_hsm*(hsm.05-t_hsm)))/(1+exp(-r_hsm*(hsm.05-t_hsm)))^2
           #K_int*r_hsm/4
    )
  return(df.mod.select)
}