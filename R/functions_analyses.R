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

get.mauri <- function(dir.occ){
  # load european species
  species.list=list.files("data/chorological_maps_dataset/")
  
  df.Euforest=read.table(dir.occ,
                         header=TRUE,
                         sep=",",
                         dec=".")%>% 
    #remove duplicated location_tree 
    filter(!duplicated(paste(X,Y, SPECIES.NAME))) %>%
    #create plot id
    mutate(Pres = rowSums(cbind(NFI==1,FF==1,BS==1)),
           Node = unclass(factor(paste(X, Y))))

  db.cont=df.Euforest %>% 
    group_by(Node,SPECIES.NAME) %>%
    summarise(n=n())%>%
    spread(SPECIES.NAME, n) %>% 
    ungroup() %>%
    pivot_longer(colnames(.)[-1],
                 names_to = "species",
                 values_to = "presence") %>% 
    filter(species %in% species.list) %>% 
    left_join(df.Euforest[!duplicated(df.Euforest$Node), c("X","Y","COUNTRY","EEO","Node")],
              by = "Node") 
  
  species.list.eff = db.cont |> 
    filter(presence==1) |> 
    group_by(species) |> 
    summarise(n=n()) |> 
    filter(n>400)
  
  db.cont<-db.cont |> 
    filter(species %in% species.list.eff$species) |> 
    rename_with(.cols=everything(),
                tolower)
  
  #Change coordinates
  coordinates(db.cont) <- c("x",  "y")
  proj4string(db.cont) <- CRS("+init=epsg:3035")
  db.cont <- as.data.frame(spTransform(db.cont, CRS("+proj=longlat +datum=WGS84")))
  
  return(db.cont)
}


get.occ.clim<-function(db.mauri,
                       clim.list,
                       file.path){
  ## Load safe ty margin indicators
  clim.list=lapply(clim.list,rast)
  df.clim=lapply(seq_along(clim.list),function(x) {
    df=extract(clim.list[[x]],data.frame(x=db.mauri$x,y=db.mauri$y))[-1]
    colnames(df)=names(clim.list)[x]
    return(df)}
    )
  
  db.mauri<-cbind(db.mauri,
                 as.data.frame(df.clim))
  
  db.mauri$presence[is.na(db.mauri$presence)] <- 0

  fwrite(db.mauri,file.path)
  # write.csv(df.fdg.sp,"output/LT50_spring.csv")
  
  return(db.mauri)
}


#' Compute safety margins on abs/pres dataset
#' 
#' @description merge with traits and compute safety margins
#' @param db.clim formated dataset
#' @param df.traits dataframe of traits associated with eahc species
#' @return dataframe, where for each location, abs/pres of all selected species
#' is mentionned, and associated predictors

get.occurence <- function(db.cont,
                          df.traits,
                          file.path){
  db.cont <- db.cont %>% 
    left_join(df.traits) 
  
  # compute LT50.sp.spring, averaged by species /// IF SPRING FROST NEEDED
  # date.dehardening=0
  # df.fdg.sp <- db.cont %>%
  #   filter(presence==1) %>% 
  #   select(species.binomial,sp.ind,lt50.mean,fdg) %>% 
  #   group_by(species.binomial,sp.ind,lt50.mean) %>% 
  #   summarise(fdg.sp.mean=mean(fdg,na.rm=TRUE),
  #             fdg.sp.sd=sd(fdg,na.rm=TRUE)) %>% 
  #   ungroup() %>% 
  #   mutate(lt50.sp.spring=-5+(30/2)*((5+lt50.mean)/(fdg.sp.mean-date.dehardening)))
  
  
  db.cont <- db.cont %>% 
    # left_join(df.fdg.sp[,c("species.binomial","lt50.sp.spring")],by="species.binomial") %>% 
    mutate(psi_old=psi,
           psi=case_when(psi<(-20000)~(-20000),
                         TRUE~psi),
           psi.100_old=psi.100,
           psi.100=case_when(psi.100<(-20000)~(-20000),
                             TRUE~psi.100)) %>% 
    mutate(hsm=psi-px.mu*1000,
           hsm.100=psi.100-px.mu*1000,
           fsm.winter=frost.winter-lt50.mean,
           # fsm.spring=frost.spring-(lt50.sp.spring),
           species=as.factor(species))
  db.cont$hsm[is.nan(db.cont$hsm)] <- NA
  db.cont$hsm.100[is.nan(db.cont$hsm.100)] <- NA
  db.cont$fsm.winter[is.nan(db.cont$fsm.winter)] <- NA
  # db.cont$fsm.spring[is.nan(db.cont$fsm.spring)] <- NA
  
  
  fwrite(db.cont,file.path)
  
  return(db.cont)
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

get.prevalence <- function(species.list,
                           db.clim){
  df.traits=data.frame(species=species.list)
  # create an empty column for prevalence
  df.traits$prevalence <- NA
  
  # for loop occuring along all species
  
  for (sp in unique(df.traits$species)){
    print(sp)
    # look for euforgen distribution if it exists for the targetted species and
    # compute the prevalelnce
    path=file.path("data",
                   "chorological_maps_dataset",
                   sp,
                   "shapefiles")
    tryCatch(
      {
        if (file.exists(path)){
          db.pres <- db.clim %>%
            filter(species==sp)
          if (length(grep("plg_clip\\.",
                        list.files(path),
                        value = TRUE))!=0){
          #load euforgen distrib
          file.list=grep("plg_clip\\.",
                         list.files(path),
                         value = TRUE)
          file.sp=unique(vapply(strsplit(file.list,"\\."), `[`, 1, FUN.VALUE=character(1)))
          } else if (length(grep("plg\\.",
                           list.files(path),
                           value = TRUE))!=0){
            #load euforgen distrib
          file.list=grep("plg\\.",
                         list.files(path),
                         value = TRUE)
          file.sp=unique(vapply(strsplit(file.list,"\\."), `[`, 1, FUN.VALUE=character(1)))
          }
          # filter Mauri db with only points of sp
          spdistrib=do.call(rbind,lapply(file.sp,function(x)read_sf(dsn=path,x))) |> 
            summarise(geometry = sf::st_union(geometry)) 
          spdistrib=st_crop(spdistrib,
                            xmin=-10,
                            xmax=28.231836,
                            ymin=st_bbox(spdistrib)$ymin[[1]],
                            ymax=st_bbox(spdistrib)$ymax[[1]])
          # create spatial points with mauri db 
          db.pres.geo <- st_as_sf(db.pres,
                                  coords=c("x","y"),
                                  crs="epsg:4326") %>% 
            st_join(spdistrib, join = st_within,left=FALSE) %>% # select only points falling in euforgen distrib
            as.data.frame() 
          df.traits[df.traits$species==sp,"prevalence"] <-  sum(db.pres.geo$presence==1)/dim(db.pres.geo)[1]
        } else { # if euforgen distrib do not exist, prevalence default set to 0.1
          df.traits[df.traits$species==sp,"prevalence"] <- 0.1
        }
        },
        error=function(e){print(paste0("error for ",sp))}
        )
  } 
  df.traits$prevalence[is.nan(df.traits$prevalence)] <- mean(df.traits$prevalence,na.rm=TRUE) #because a value is always needed for mod prior
  df.traits$prevalence[is.na(df.traits$prevalence)] <- mean(df.traits$prevalence,na.rm=TRUE)
  return(df.traits[,c("species","prevalence")])
}



#' Compute niche caracteristics
#' 
#' @description compute mean lat/long and continentality of species niche
#' @note Continentality from https://doi.org/10.1038/s41597-020-0464-0
#' @param df.traits df of traits for each species
#' @param db.clim database of pres/abs
#' @return df.loc + extracted variables + computed variables
#' 

get.niche <- function(species.list,
                      db.clim,
                      psi="psi_cerraday_real",
                      tmin="tmin_cerra"){
  df.traits=data.frame(species=species.list)
  
  df.niche <- db.clim %>% 
    filter(presence==1) %>% 
    group_by(species) %>% 
    summarise(lat.mean=mean(y),
              lat.sd=sd(y),
              long.mean=mean(x),
              long.sd=sd(x),
              lat.q05=quantile(y,prob=0.05)[[1]],
              lat.q95=quantile(y,prob=0.95)[[1]],
              psi.q05=quantile(eval(parse(text=psi)),prob=0.05,na.rm=TRUE)[[1]],
              psi.q95=quantile(eval(parse(text=psi)),prob=0.95,na.rm=TRUE)[[1]],
              tmin.q05=quantile(eval(parse(text=tmin)),prob=0.05,na.rm=TRUE)[[1]],
              tmin.q95=quantile(eval(parse(text=tmin)),prob=0.95,na.rm=TRUE)[[1]]
    )
  
  
  # rast.cont=as.data.frame(mean(rast("data/jci_year.nc")),xy=TRUE) %>% 
  #   mutate(z=x,
  #          x=y,
  #          y=z) %>% 
  #   select(-z) %>% 
  #   rast(crs="epsg:4326")
  # df.traits$jci=NA
  # for (sp in df.traits$species.binomial){
  #   print(sp)
  #   db.pres <- db.clim %>% 
  #     filter(species.binomial==sp) %>% 
  #     filter(presence==1)
  #   db.pres <- cbind(db.pres,
  #                    jci=extract(rast.cont,db.pres[,c("x","y")])[["mean"]])
  #   jci=mean(db.pres$jci,na.rm=TRUE)
  #   df.traits[df.traits$species.binomial==sp,"jci"]=jci
  # }
  
  # df.overlap <- db.clim %>%
  #   filter(presence==1) %>%
  #   group_by(node) %>%
  #   mutate(nb.sp=sum(presence==1),
  #          overlap.hsm=sum(hsm>0),
  #          overlap.fsm=sum(fsm.winter>0)) %>%
  #   ungroup() %>%
  #   group_by(species.binomial) %>%
  #   summarise(overlap=mean(nb.sp),
  #             overlap.hsm=mean(overlap.hsm,na.rm=TRUE),
  #             overlap.fsm=mean(overlap.fsm,na.rm=TRUE))

  df.traits <- df.traits %>%
    left_join(df.niche,by="species") #%>%
    # left_join(df.overlap,by="species.binomial")
  return(df.traits[,c("species","lat.mean","lat.sd","long.mean","long.sd",
                      "lat.q05","lat.q95","psi.q05","psi.q95","tmin.q05","tmin.q95")])#,"jci","overlap","overlap.hsm","overlap.fsm"
}

#' Get shade tol
#' 
#' @description compute mean lat/long and continentality of species niche
#' @note Continentality from https://doi.org/10.1038/s41597-020-0464-0
#' @param df.traits df of traits for each species
#' @return df.loc + extracted variables + computed variables
#' 

get.shadetol <- function(species.list,
                         db.clim){
  df.traits=data.frame(species=species.list)
  df.shadetol=read.csv2("data/Species traits/data_Niinemets&Valladares_2006.csv")
  df.traits=df.traits %>% 
    left_join(df.shadetol,by=c("species"="Species")) |> 
    mutate(across(c(shade_tolerance.mean,drought_tolerance.mean,waterlogging_tolerance.mean),
                  as.numeric))
  
  df.overlapshade=db.clim |>
    filter(presence==1) |> 
    select(node,species,presence) |> 
    left_join(df.traits[,c("species","shade_tolerance.mean","drought_tolerance.mean")], by="species") |> 
    filter(!is.na(shade_tolerance.mean)&
             !is.na(drought_tolerance.mean)) |> 
    group_by(node) |> 
    mutate(overlap_plot_shade=colSums(outer(shade_tolerance.mean,shade_tolerance.mean,">")),
           overlap_plot_drought=colSums(outer(drought_tolerance.mean,drought_tolerance.mean,">"))) |> #shade_overlap(shade_tolerance.mean,presence) 
    ungroup() |> 
    group_by(species) |> 
    summarise(overlap_shade=mean(overlap_plot_shade),
              overlap_drought=mean(overlap_plot_drought))
  
  df.traits=df.traits %>% 
    left_join(df.overlapshade,by=c("species"))
  return(df.traits[,c("species","shade_tolerance.mean","drought_tolerance.mean","waterlogging_tolerance.mean","overlap_shade","overlap_drought")]) 
}

#' Get traits and caract
#' 
#' @description compute niche traits and caract
#' @param df.traits df of traits for each species
#' @return df.loc + extracted variables + computed variables
#' 

get.species <- function(species.list,
                        df.preval,
                        df.shadetol,
                        df.niche,
                        file.output){
  df.species <- data.frame(species=species.list) %>% 
    left_join(df.shadetol,by="species")%>% 
    left_join(df.niche,by="species")%>% 
    left_join(df.preval,by="species")
  fwrite(df.species,file=file.output)
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
#' @param soil.depth "real" or "constant"
#' @return no return, fit_mod6 is the folder with all fitted models
#' 

fit.logistic <- function(db.clim.file,
                         df.species,
                         soil.depth="real",
                         output){
  
  # read db.clim and filter some species
  db.clim=fread(db.clim.file) %>% 
      # filter(!species.binomial %in% c("Juniperus virginiana",
      #                                 "Prunus persica",
      #                                 "Quercus palustris",
      #                                 "Cedrus libani",
      #                                 "Liriodendron tulipifera",
      #                                 "Pinus ponderosa",
      #                                 "Alnus cordata",
      #                                 "Acer negundo",
      #                                 "Prunus cerasifera",
      #                                 "Tsuga heterophylla",
      #                                 "Prunus serotina",
      #                                 "Juglans nigra",
      #                                 "Quercus rubra",
      #                                 "Pinus contorta",
      #                                 "Abies grandis",
      #                                 "Juglans regia",
      #                                 "Abies nordmanniana",
      #                                 "Pseudotsuga menziesii",
      #                                 "Robinia pseudoacacia"
      # )) |> 
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm.winter)) |> 
    filter(psi>(-10000)) |>  #remove very low value of psi
    mutate(hsm=hsm/1000,
           hsm.100=hsm.100/1000,
           fsm=fsm.winter) 
  
  if (soil.depth=="constant"){
    db.clim <- db.clim |> 
      select(-hsm) |> 
      rename(`hsm`="hsm.100")
  }
  # set species list, as sp present in occ data
  species.list=unique(db.clim$species.binomial)
  

  for (sp in species.list){
    print(sp)
    db.pres <- db.clim %>%
      filter(species.binomial==sp) 
    
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
    save(fit.2var,file=paste0(output,sp,"_fitIII.RData"))
    
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
    save(fit.hsm,file=paste0(output,sp,"_fitIII_hsm.RData"))
    
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
    save(fit.fsm,file=paste0(output,sp,"_fitIII_fsm.RData"))
    
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
    save(fit.none,file=paste0(output,sp,"_fitIII_none.RData"))
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
                       output.path,
                       file.path){
  # read db.clim 
  db.clim=fread(db.clim.file) %>% 
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm.winter)) |> 
    mutate(hsm=hsm/1000,
           hsm.100=hsm.100/1000,
           fsm=fsm.winter) 

  
  # set species list, as sp present in occ data
  species.list=unique(db.clim$species.binomial)
  
  # define a dataframe to 
  df.output=df.species %>% 
    filter(species.binomial %in% species.list) %>% 
    select(species.binomial,prevalence) %>% 
    mutate(k_int=NA,
           r_fsm=NA,
           r_hsm=NA,
           t_hsm=NA,
           t_fsm=NA,
           rhat=NA,
           lp__=NA,
           divergence=NA) %>% 
    crossing(mod=c("2sm","hsm","fsm","none"))
  
  for (sp in species.list){
    print(sp)
    tryCatch(
      {# 2 var
      print("2sm")
      load(file=paste0(output.path,sp,"_fitIII.RData"))
      summary=summary(fit.2var)$summary
      par_mod=get_sampler_params(fit.2var)
      df.output[df.output$species.binomial==sp
                & df.output$mod=="2sm",
                c("k_int","r_fsm","r_hsm","t_hsm","t_fsm","rhat","lp__","divergence")]=
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
                c("k_int","r_fsm","r_hsm","t_hsm","t_fsm","rhat","lp__","divergence")]=
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
                c("k_int","r_fsm","r_hsm","t_hsm","t_fsm","rhat","lp__","divergence")]=
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
                c("k_int","r_fsm","r_hsm","t_hsm","t_fsm","rhat","lp__","divergence")]=
        list(summary["K_int","mean"],
             0,
             0,
             0,
             0,
             max(summary[,"Rhat"]),
             summary["lp__","mean"],
             mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
        )
      },
      error=function(e){print(paste0("error for ",sp))}
    )
  }
  
  #' Compute BIC
  df.output <- df.output %>% 
    mutate(nb.par=case_when(mod=="2sm"~5,
                            mod=="hsm"|mod=="fsm"~3,
                            mod=="none"~1),
           bic=2*nb.par-2*lp__)
  fwrite(df.output,file=file.path)
  return(df.output)
}


#' Build output with uncertainties
#' 
#' @description load each model and extract relevant parameters
#' @param df.species df of traits for each species
#' @param output.path
#' @return df.output
#'

get.output.uncertainties <- function(output.path,
                                     file.path){
  # Get a list of all .RData files in the directory
  files <- list.files(path = output.path, pattern = "*.RData")
  
  # Initialize an empty list to store the data frames
  list_of_tables <- list()
  
  # Loop over the files
  for (i in 133:140) {
    print(i)
    tryCatch({
      new_env <- new.env()
      load(file=paste0(output.path, files[i]),envir=new_env)
      obj_name <- ls(new_env)[1]
      model <- get(obj_name, envir = new_env)
      par_mod=rstan::get_sampler_params(model)
      rm(new_env)
      
      summary=as.data.frame(summary(model)$summary)
      summary$div=mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
          
      # Add columns for species and model
      summary$species <- strsplit(files[i], "_")[[1]][1]
      summary$model <- if_else(is.na(strsplit(files[i], "_")[[1]][3]),"2sm",strsplit(files[i], "_")[[1]][3])
      summary$model <- gsub(".RData", "", summary$model) # Remove the .RData from the model name
      
      # Add the modified table to the list
      list_of_tables[[i]] <- summary
      rm(model,par_mod,summary)
    },
    error=function(e){print(paste0("error for ",i))}
    )
    
  }
  
  # Use lapply to apply the same function to each table
  reshaped_tables <- lapply(list_of_tables, function(df) {
    tryCatch({
      df <- as.data.frame(df) |> tibble::rownames_to_column(var="parameter")
      df <- df %>% pivot_longer(cols = c(-parameter,-species,-model), names_to = "posterior", values_to = "quantile")
    },
    error=function(e){print(paste0("error"))}
    )
    return(df)
  })
  
  # Bind all the reshaped tables together
  df.mod.output <- bind_rows(reshaped_tables) |> 
    rename(mod=model)
  bic <- df.mod.output %>% 
    filter(posterior=="mean"&parameter=="lp__") |> 
    mutate(nb.par=case_when(mod=="2sm"~5,
                            mod=="hsm"|mod=="fsm"~3,
                            mod=="none"~1),
           bic=2*nb.par-2*quantile) |> 
    select(species,mod,bic)
  df.mod.output<-df.mod.output |> 
    left_join(bic,by=c("species","mod"))
  
  
  df.mod.output |> filter(posterior%in% c("Rhat","div")&parameter=="lp__") |> 
    pivot_wider(names_from = posterior,
                values_from = quantile) |> 
    filter(Rhat<1.2) %>% 
    filter(div <0.1) %>% 
    group_by(species) %>% 
    slice(which.min(bic)) %>% 
    ungroup() -> df.output.select
  return(df.mod.output)
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
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm.winter)) |> 
    filter(psi>(-10000)) |>  #remove very low value of psi
    mutate(hsm=hsm/1000,
           hsm.100=hsm.100/1000,
           fsm=fsm.winter) 
  
  
  # set species list, as sp present in occ data
  species.list=unique(df.output$species.binomial)
  
  df.output$auc=NA
  df.output$threshold_auc=NA
  df.output$tss=NA
  df.output$threshold_tss=NA
  
  
  for (sp in species.list){
    print(sp)
    
    ### 2 var
    print("2sm")
    k_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="2sm",
                               "k_int"])
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
      mutate(proba=k_int/
               ((1+exp(-r_fsm*(fsm.winter-t_fsm)))*
                  (1+exp(-r_hsm*(hsm-t_hsm)))),
             pred=as.integer(case_when(proba>proba_thresh~1,
                            TRUE~0)))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA,
                      tss=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
      df.auc[df.auc$threshold==i,"tss"]=confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Sensitivity"]] +
                                        confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Specificity"]] -
                                        1
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="2sm",
              c("threshold_auc","auc","threshold_tss","tss")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"],
           df.auc[which.max(df.auc$tss),"threshold"],
           df.auc[which.max(df.auc$tss),"tss"])
    rm(db.pred,db.pred.auc,df.auc)
    
    # HSM
    print("hsm")
    k_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="hsm",
                               "k_int"])
    r_hsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="hsm",
                               "r_hsm"])
    t_hsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="hsm",
                               "t_hsm"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=k_int/
               ((1+exp(-r_hsm*(hsm-t_hsm)))),
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA,
                      tss=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
      df.auc[df.auc$threshold==i,"tss"]=confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Sensitivity"]] +
        confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Specificity"]] -
        1
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="hsm",
              c("threshold_auc","auc","threshold_tss","tss")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"],
           df.auc[which.max(df.auc$tss),"threshold"],
           df.auc[which.max(df.auc$tss),"tss"])
    rm(db.pred,db.pred.auc,df.auc)
    
    
    # FSM
    print("fsm")
    k_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="fsm",
                               "k_int"])
    r_fsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="fsm",
                               "r_fsm"])
    t_fsm=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="fsm",
                               "t_fsm"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=k_int/
               ((1+exp(-r_fsm*(fsm.winter-t_fsm)))),
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA,
                      tss=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
      df.auc[df.auc$threshold==i,"tss"]=confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Sensitivity"]] +
        confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Specificity"]] -
        1
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="fsm",
              c("threshold_auc","auc","threshold_tss","tss")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"],
           df.auc[which.max(df.auc$tss),"threshold"],
           df.auc[which.max(df.auc$tss),"tss"])
    rm(db.pred,db.pred.auc,df.auc)
    
    
    # None
    print("none")
    k_int=as.numeric(df.output[df.output$species.binomial==sp
                               & df.output$mod=="none",
                               "k_int"])
    db.pred <- db.clim %>% 
      filter(species.binomial==sp) %>% 
      crossing(proba_thresh=seq(0.01,1,by=0.01)) %>% 
      mutate(proba=k_int,
             pred=case_when(proba>proba_thresh~1,
                            TRUE~0))
    df.auc=data.frame(threshold=seq(0.01,1,by=0.01),
                      auc=NA,
                      tss=NA)
    for (i in seq(0.01,1,by=0.01)){
      db.pred.auc <- db.pred %>% 
        filter(proba_thresh==i)
      df.auc[df.auc$threshold==i,"auc"]=as.numeric(pROC::auc(db.pred.auc$presence,db.pred.auc$pred))
      df.auc[df.auc$threshold==i,"tss"]=confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Sensitivity"]] +
        confusionMatrix(as.factor(db.pred.auc$pred),as.factor(db.pred.auc$presence))$byClass[["Specificity"]] -
        1
    }
    df.output[df.output$species.binomial==sp
              & df.output$mod=="none",
              c("threshold_auc","auc","threshold_tss","tss")]=
      list(df.auc[which.max(df.auc$auc),"threshold"],
           df.auc[which.max(df.auc$auc),"auc"],
           df.auc[which.max(df.auc$tss),"threshold"],
           df.auc[which.max(df.auc$tss),"tss"])
    rm(db.pred,db.pred.auc,df.auc)
    
    
  }
  
  return(df.output)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Fit model with climate ####
#' @description Fit model with classical climatic variable
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Fit model
#' 
#' @description function that fit 4 models for each species and save it
#' @param df.species df of traits for each species
#' @param soil.depth "real" or "constant"
#' @return no return, fit_mod6 is the folder with all fitted models
#' 

fit.logistic.clim <- function(db.clim.file,
                         df.species,
                         output,
                         file.path){
  # check that ouput directory exists
  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
    cat("Folder created at:", output, "\n")
  } else {
    cat("Folder already exists at:", output, "\n")
  }
  
  
  
  # read db.clim and filter some species
  db.clim=fread(db.clim.file) %>% 
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm.winter)) |> 
    filter(psi>(-10000)) |>  #remove very low value of psi
    mutate(hsm=hsm/1000,
           hsm.100=hsm.100/1000,
           fsm=fsm.winter) |> 
    filter(!is.na(wai)) |> 
    filter(!is.na(temp.mean))
    
  
  
  # set species list, as sp present in occ data
  species.list=unique(db.clim$species.binomial)
  
  # prepare dataframe output
  df.output=df.species %>% 
    filter(species.binomial %in% species.list) %>% 
    select(species.binomial,prevalence) %>% 
    mutate(k_int=NA,
           r_mat=NA,
           r_wai=NA,
           t_mat=NA,
           t_wai=NA,
           rhat=NA,
           lp__=NA,
           divergence=NA) %>% 
    crossing(mod=c("2var","wai","mat"))
  
  
  # for loop for model fitting, no none model because already fitted with safety
  # margins
  for (sp in species.list){
    print(sp)
    
    #initialize
    db.pres <- db.clim %>%
      filter(species.binomial==sp) 
    mean.mat=mean(db.pres$temp.mean)
    mean.wai=mean(db.pres$wai)
    preval=df.species[df.species$species.binomial==sp,
                      "prevalence"]
    
    # 2 var
    print("2var")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      mat=as.numeric(db.pres$temp.mean),
                      wai=as.numeric(db.pres$wai),
                      prior_K=preval,
                      mean_mat=mean.mat,
                      mean_wai=mean.wai,
                      NULL)
    fit.2var <- stan(file = "glm_log_1sp_2clim.stan",
                     data=data.list,
                     iter=1000,
                     chains=3,
                     core=3,
                     include=FALSE,
                     pars=c("proba","K_vect"))
    save(fit.2var,file=paste0(output,sp,"_fit_2clim.RData"))
    
    summary=summary(fit.2var)$summary
    par_mod=get_sampler_params(fit.2var)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="2var",
              c("k_int","r_mat","r_wai","t_wai","t_mat","rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           summary["r_mat","mean"],
           summary["r_wai","mean"],
           summary["t_wai","mean"],
           summary["t_mat","mean"],
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
    
    
    # WAi (equivalent HSM)
    print("wai")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      wai=as.numeric(db.pres$wai),
                      prior_K=preval,
                      mean_wai=mean.wai,
                      NULL)
    fit.wai <- stan(file = "glm_log_1sp_wai.stan",
                    data=data.list,
                    iter=1000,
                    chains=3,
                    core=3,
                    include=FALSE,
                    pars=c("proba","K_vect"))
    save(fit.wai,file=paste0(output,sp,"_fit_wai.RData"))
    
    summary=summary(fit.wai)$summary
    par_mod=get_sampler_params(fit.wai)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="wai",
              c("k_int","r_mat","r_wai","t_wai","t_mat","rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           0,
           summary["r_wai","mean"],
           summary["t_wai","mean"],
           0,
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
    
    
    # MAT (equivalent FSM)
    print("mat")
    data.list <- list(N=dim(db.pres)[1],
                      presence=db.pres$presence,
                      mat=as.numeric(db.pres$temp.mean),
                      prior_K=preval,
                      mean_mat=mean.mat,
                      NULL)
    fit.mat <- stan(file = "glm_log_1sp_mat.stan",
                    data=data.list,
                    iter=1000,
                    chains=3,
                    core=3,
                    include=FALSE,
                    pars=c("proba","K_vect"))
    save(fit.mat,file=paste0(output,sp,"_fit_mat.RData"))
    
    summary=summary(fit.mat)$summary
    par_mod=get_sampler_params(fit.mat)
    df.output[df.output$species.binomial==sp
              & df.output$mod=="mat",
              c("k_int","r_mat","r_wai","t_wai","t_mat","rhat","lp__","divergence")]=
      list(summary["K_int","mean"],
           summary["r_mat","mean"],
           0,
           0,
           summary["t_mat","mean"],
           max(summary[,"Rhat"]),
           summary["lp__","mean"],
           mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
      )
  }
  
  
  #' Compute BIC & inflexion index
  wai.95=quantile(db.clim$wai,prob=0.95)[[1]]
  mat.95=quantile(db.clim$temp.mean,prob=0.95)[[1]]
  df.output <- df.output %>% 
    mutate(nb.par=case_when(mod=="2var"~5,
                            mod=="wai"|mod=="mat"~3,
                            mod=="none"~1),
           bic=2*nb.par-2*lp__) |> 
    mutate(inflex_mat=100-100*((16*exp(-r_mat*(mat.95-t_mat)))/(2+2*exp(-r_mat*(mat.95-t_mat)))^2),
           inflex_wai=100-100*((16*exp(-r_wai*(wai.95-t_wai)))/(2+2*exp(-r_wai*(wai.95-t_wai)))^2))
  fwrite(df.output,file=file.path)
  
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 5 - Model selection ####
#' @description Select best model for each species 
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Model selection
#' 
#' @description compute auc for each model
#' @param df.output
#' @return df.output
#'

select.model <- function(df.output){
  #' Model selection
  df.mod.select <- df.output %>% 
    filter(rhat<1.2) %>% 
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


fit.allspecies<- function(db.clim.file,
                                  output.clim.file,
                                  output,
                                  mod.folder){
  # check that ouput directory exists
  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
    cat("Folder created at:", output, "\n")
  } else {
    cat("Folder already exists at:", output, "\n")
  }
  
  if (!dir.exists(mod.folder)) {
    dir.create(mod.folder, recursive = TRUE)
    cat("Folder created at:", mod.folder, "\n")
  } else {
    cat("Folder already exists at:", mod.folder, "\n")
  }
  
  # read db.clim and filter some species
  db.clim=fread(db.clim.file) %>% 
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm.winter)) |> 
    filter(psi>(-10000)) |>  #remove very low value of psi
    mutate(hsm=hsm/1000,
           hsm.100=hsm.100/1000,
           fsm=fsm.winter) |> 
    filter(!is.na(wai)) |> 
    filter(!is.na(temp.mean))
  
  
  
  # set species list, as sp present in occ data
  # species.list=db.clim |> 
  #   filter(presence==1) |> 
  #   group_by(species.binomial) |> 
  #   summarise(n=n()) |> 
  #   arrange(n)
  # n.pres=quantile(species.list$n,probs=0.1)
  # 
  # 
  # db.clim.sub=db.clim |>
  #   sample_n(dim(db.clim)[1]/100,replace=TRUE)
  if(!file.exists(paths=paste0(mod.folder,"fit_random_allsp.RData"))){
    data.list<-list(N=dim(db.clim)[1],
                  S=nlevels(as.factor(db.clim$species.binomial)),
                  presence=db.clim$presence,
                  species=as.numeric(as.factor(db.clim$species.binomial)),
                  fsm=db.clim$fsm,
                  hsm=db.clim$hsm)
  
    fit.allsp <- stan(file = "glm_log_all.stan",
                      data=data.list,
                      iter=1000,
                      chains=3,
                      core=3,
                      include=FALSE,
                      pars=c("proba","K_vect"))
    save(fit.allsp,file=paste0(mod.folder,"fit_random_allsp.RData"))
  } else{
    print("discarded safmarg model fitting")
    load(paste0(mod.folder,"fit_random_allsp.RData"))
  }
    
  if(!file.exists(paths=paste0(mod.folder,"fit_random_allspClim.RData"))){
    outputclim=fread(output.clim.file) |>  #output.clim.file="output/df.outputClim.csv"
      filter(rhat<1.2) %>% 
      filter(divergence <0.1) %>% 
      group_by(species.binomial) %>% 
      slice(which.min(bic)) %>% 
      ungroup() |> 
      mutate(across(c("r_mat","r_wai","t_wai","t_mat"),
                    ~na_if(.,0)))|>
      mutate(across(where(is.numeric),
                    ~mean(.,na.rm=TRUE))) |> 
      select(-species.binomial,-mod) |> unique()
    # db.clim<-db.clim |> filter(species.binomial %in% c("Abies alba","Fagus sylvatica","Picea abies")) |> sample_n(70000)
    data.list<-list(N=dim(db.clim)[1],
                    S=nlevels(as.factor(db.clim$species.binomial)),
                    presence=db.clim$presence,
                    species=as.numeric(as.factor(db.clim$species.binomial)),
                    mat=db.clim$temp.mean,
                    wai=db.clim$wai,
                    prior=as.numeric(outputclim[,c("r_mat","r_wai","t_mat","t_wai")])
                    # mean_mat=mean(db.clim.sub$temp.mean),
                    # mean_wai=mean(db.clim.sub$wai)
    )
    fit.allspClim <- stan(file = "glm_log_allClim.stan",
                          data=data.list,
                          iter=1000,
                          chains=3,
                          core=3,
                          include=FALSE,
                          pars=c("proba","K_vect"))
    save(fit.allspClim,file=paste0(mod.folder,"fit_random_allspClim.RData"))
  } else{
    print("discarded clim model fitting")
    load(paste0(mod.folder,"fit_random_allspClim.RData"))
  }
  
  posteriors_mean_sfm<-as.data.frame(t(summary(fit.allsp)$summary)) |> 
    select(!matches("K_sp"))
  posteriors_mean_clim<-as.data.frame(t(summary(fit.allspClim)$summary)) |> 
    select(!matches("K_sp"))
  
  db.clim.auc<-db.clim |> 
    # filter(species.binomial==sp) |> 
    select(species.binomial,presence,x,y,hsm,fsm,temp.mean,wai) |> 
    mutate(pred_sfm=posteriors_mean_sfm$K_int[1]/
             ((1+exp(-posteriors_mean_sfm$r_fsm[1]*(fsm-posteriors_mean_sfm$t_fsm[1])))*
                (1+exp(-posteriors_mean_sfm$r_hsm[1]*(hsm-posteriors_mean_sfm$t_hsm[1])))),
           pred_clim=posteriors_mean_clim$K_int[1]/
             ((1+exp(-posteriors_mean_clim$r_mat[1]*(temp.mean-posteriors_mean_clim$t_mat[1])))*
                (1+exp(-posteriors_mean_clim$r_wai[1]*(wai-posteriors_mean_clim$t_wai[1]))))
    ) |> 
    group_by(species.binomial) |> 
    summarize(auc_sfm=as.numeric(auc(presence,pred_sfm)),
              auc_clim=as.numeric(auc(presence,pred_clim)))
  # prior sur le range des aprametres du mod par esp
  fwrite(db.clim.auc,file=paste0(output,"db.clim.auc.csv"))
  return(db.clim.auc)
}
  

