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
#### Section 1 - Generalized linear mixte model ####
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

sm.glmm.mauri <- function(safety.margins,
                          dir.occ="data/EUForestsMauri/EUForestspecies.csv",
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
  
  #Load Mauri occurence data
  db.tree.mauri=read.table(file=dir.occ,
                           header=TRUE,
                           sep=",",
                           dec=".") %>%
    #remove duplicated location_tree 
    filter(!duplicated(paste(X,Y, SPECIES.NAME))) %>%
    #create plot id
    mutate(Pres = rowSums(cbind(NFI==1,FF==1,EEO==1)),
           Node = unclass(factor(paste(X, Y)))) %>% 
    #filter targetted species
    filter(SPECIES.NAME %in% as.character(df.species$species.binomial))

  # Pres/abs dataset
  db.cont=db.tree.mauri %>% 
    group_by(Node,SPECIES.NAME) %>%
    summarise(n=n())%>%
    spread(SPECIES.NAME, n) %>% 
    ungroup() %>%
    pivot_longer(colnames(.)[-1],
                 names_to = "species.binomial",
                 values_to = "presence") %>% 
    left_join(db.tree.mauri[!duplicated(db.tree.mauri$Node), c("X","Y","COUNTRY","NFI","FF","BS","EEO","Node")],
              by = "Node") %>% 
    left_join(df.species[,c(1,4,5)],
              by="species.binomial" )
  #Change coordinates
  coordinates(db.cont) <- c("X",  "Y")
  proj4string(db.cont) <- CRS("+init=epsg:3035")
  db.cont <- as.data.frame(spTransform(db.cont, CRS("+proj=longlat +datum=WGS84")))
  
  #Extract safety margins
  db.cont <- cbind(db.cont,
                   extract(rast(safety.margins$psimin,crs="epsg:4326")[["psi"]],
                           data.frame(x=db.cont$X,y=db.cont$Y))[-1],
                   extract(rast(safety.margins$tmin,crs="epsg:4326")[["t.min"]],
                           data.frame(x=db.cont$X,y=db.cont$Y))[-1]) %>% 
    mutate(fsm=t.min-LT50,
           hsm=psi-P50*1000,
           species.binomial=as.factor(species.binomial))
  db.cont$presence[is.na(db.cont$presence)] <- 0
  db.cont$fsm[is.nan(db.cont$fsm)] <- NA
  db.cont$hsm[is.nan(db.cont$hsm)] <- NA
  
  # Extract climatic variables
  db.clim <- cbind(
    db.cont,
    get_waisgdd(df.loc=data.frame(x=db.cont$X,
                                  y=db.cont$Y)[,3:5])
  )
  
  
  # Pres/abs filtered
  hsm.q001=quantile(db.cont$hsm,0.001,na.rm=TRUE)
  db.pres <- db.cont %>%
    filter(species.binomial%in%c("Acer campestre","Acer monspessulanum","Fagus sylvatica","Fraxinus excelsior","Quercus ilex")) %>% 
    filter(hsm>hsm.q001) %>% 
    mutate(hsm=scale(hsm,scale=TRUE,center=TRUE),
           fsm=scale(fsm,scale=TRUE,center=TRUE))

  
#%%%%%%%%%%%%%%
#### explo ####
#%%%%%%%%%%%%%%

  db.pres %>% 
    ggplot(aes(fsm,presence))+ 
    geom_point()+
    geom_smooth()+
    facet_wrap(~species.binomial)
  db.pres %>% 
    ggplot(aes(hsm,presence))+ 
    geom_point()+
    geom_smooth()+
    facet_wrap(~species.binomial)
  db.pres %>% 
    mutate(presence=as.factor(presence)) %>% 
    ggplot(aes(presence,hsm))+
    geom_boxplot()
  
  db.cont %>% 
    ggplot(aes(psi))+geom_density()
  
corrplot <- db.cont %>% 
    select(-P50,-NFI,-FF,-BS,-EEO,-LT50,-COUNTRY) %>% 
    filter(!is.na(psi)) %>% 
    filter(!is.na(t.min)) %>% 
    # filter(species.binomial%in%c("Fagus sylvatica","Pinus sylvestris","Quercus ilex")) %>% 
    filter(psi>quantile(db.cont$psi,probs=0.01,na.rm=TRUE)) %>% 
    mutate(hsm=hsm/1000,
           fsm.class=cut(fsm,
                         breaks=c(-Inf,0,5,10,12,14,16,18,20,25,30, Inf)),
           hsm.class=cut(hsm,
                         breaks=c(-Inf,-20,-15,-10,-7.5,-5,-2.5,0,1,2,3,4,5,7.5,10,15,Inf))) %>%
    group_by(species.binomial,fsm.class,hsm.class) %>% 
    mutate(ntot=n()) %>% 
    filter(presence!=0) %>% 
    mutate(npres=n(),
           prop=npres/ntot) %>% 
    select(species.binomial,fsm.class,hsm.class,prop,npres) %>% 
    ungroup() %>% 
    distinct() %>% 
  ggplot(aes(x=fsm.class,y=hsm.class))+
    geom_tile(aes(fill=prop))+
    scale_fill_distiller(palette="YlGnBu",direction=1)+
    geom_point(aes(size=npres))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 0.6))+
  xlab("FSM class (celsius)")+
  ylab("HSM class (MPa)")+
  labs(fill="Presence proportion",
       size="Number of observations")+
  facet_wrap(~species.binomial)
corrplot
  ggsave(corrplot,
       filename = paste0("corrplot_all.png"),
       path="figs/",
       device="png",
       scale=2)
  
#%%%%%%%%
#### end
#%%%%%%%%
  
  # Fit glmm
  sm.glmm=glmer(presence ~  hsm + fsm + (1|species.binomial), 
                data = db.pres,
                family=binomial(link = "logit"))
  for (species in unique(db.cont$species.binomial)){
    hsm.q001=quantile(db.cont$hsm,0.001,na.rm=TRUE)
    db.pres.sp<- db.cont %>%
      filter(species.binomial==species) %>% 
      filter(hsm>hsm.q001) %>% 
      mutate(hsm=scale(hsm,scale=TRUE,center=TRUE),
             fsm=scale(fsm,scale=TRUE,center=TRUE))
    smsp.glmm=glm(presence ~  hsm + fsm, 
                    data = db.pres.sp,
                    family=binomial(link = "logit"))
    print(species)
    print(summary(smsp.glmm))
  }

  
  
}
