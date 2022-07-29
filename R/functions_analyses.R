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

sm.glmm <- function(db.ifn,
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
