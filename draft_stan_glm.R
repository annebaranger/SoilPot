#%%%%%%%%%%%%%%%%%%
#### Section 1 ####
#'@description load packages needed


# Load dataset
library("data.table")
library("tidyr")
library(dplyr)
library(ggplot2)
library(mcp)
library(rstan)
library("boot")
library(stringr)

#%%%%%%%%%%%%%%%%%%
#### Section 2 ####
#'@description load datasets 
db.clim=fread("output/db_EuForest.csv") %>% 
  rename_with(.cols=everything(),
              tolower) %>% 
  filter(!is.na(hsm)) %>% 
  filter(!is.na(fsmwinter)) %>% 
  filter(!is.na(fsmspring)) %>% 
  filter(!is.na(sgdd)) %>% 
  filter(!is.na(wai)) %>% 
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
df.traits=read.csv("output/df_trait_filtered.csv") %>% 
  slice(-5)

#%%%%%%%%%%%%%%%
#### Section 2.1 ####
#'@description computing prevalence data
#' example for betula pendula
library(sf)
df.traits$prevalence <- NA
for (sp in unique(df.traits$species.binomial)){
  if (!is.na(df.traits[df.traits$species.binomial==sp,"file"])){
    db.pres <- db.clim %>% 
      filter(species.binomial==sp)
    file.sp=df.traits[df.traits$species.binomial==sp,"file"]
    spdistrib=read_sf(dsn=file.path("data",
                                    "chorological_maps_dataset",
                                    sp,
                                    "shapefiles"),
                      layer=file.sp)
    db.pres.geo <- st_as_sf(db.pres,
                            coords=c("x","y"),
                            crs="epsg:4326") %>% 
      st_join(spdistrib, join = st_within,left=FALSE) %>% 
      as.data.frame() 
    df.traits[df.traits$species.binomial==sp,"prevalence"] <-  sum(db.pres.geo$presence==1)/dim(db.pres.geo)[1]
  } else {
    df.traits[df.traits$species.binomial==sp,"prevalence"] <- 0.1
  }
}

#%%%%%%%%%%%%%%%%%%
#### Section 3 ####
#'@description format dataset
#'@output db.pres cleaned dataframe without NA and subset from original dataset

# Function pseudolog
pseudoLog10 <- function(x) {asinh(x/2)/log(10)}
#hsm.q01=quantile(db.clim$HSM,0.01,na.rm=TRUE) #because hsm ceilled at -20MPa

# Filter targetted species adn subset 1% 
set.seed(4554)
db.pres <- db.clim %>%
  # filter(species.binomial%in%c("Fagus sylvatica","Pinus sylvestris","Fraxinus excelsior","Quercus ilex")) %>% #"
  filter(species.binomial %in% c("Abies alba","Quercus ilex","Olea europaea",
                                 "Pinus pinaster","Quercus robur","Fagus sylvatica",
                                 "Betula pendula","Picea abies")) %>%
  # filter(species.binomial=="Fagus sylvatica") %>% 
  #filter(hsm>hsm.q01) %>% 
  sample_frac(0.05)
# rm(db.clim)

hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
#hsm.plog.mu=mean(pseudoLog10(db.pres$hsm),na.rm=TRUE)
hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
#hsm.plog.sd=sd(pseudoLog10(db.pres$hsm),na.rm=TRUE)
fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)

# Scale variables
db.pres <- db.pres %>% 
  mutate(hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
         # hsm.pl10=pseudoLog10(hsm),
         fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1],
         pet=scale(pet,scale=TRUE,center=TRUE)[, 1],
         sgdd=scale(sgdd,scale=TRUE,center=TRUE)[, 1]) 


#%%%%%%%%%%%%%%%%%%
#### Section 4 ####
#'@description few exploratory plots

# db.pres %>% 
#   mutate(hsm.ps.sc=pseudoLog10(scale(hsm,scale=TRUE,center=TRUE)),
#          hsm.sc.ps=scale(pseudoLog10(hsm),scale=TRUE,center=TRUE)) %>%
#   ggplot()+
#   geom_point(aes(hsm,hsm.ps.sc),color="green")+
#   geom_point(aes(hsm,hsm.sc.ps),color="red")+
#   geom_density(aes(hsm.ps.sc),color="red") + 
#   geom_hline(yintercept = -0.07497452,color="green")+
#   geom_hline(yintercept = -hsm.plog.mu/hsm.plog.sd,color="red")+
#   geom_vline(xintercept=0)
# 
# db.pres %>%
#   filter(species.binomial=="Fagus sylvatica") %>% 
#   sample_frac(0.1) %>% 
#   ggplot(aes(fsmwinter,presence,color=species.binomial))+
#   geom_point()+
#   geom_vline(xintercept = 0)

# parameters of scaling
hsm.bkpt=-hsm.mu/hsm.sd
fsm.bkpt=-fsm.mu/fsm.sd
hsm.plog.bkpt=pseudoLog10(-hsm.mu/hsm.sd)


#%%%%%%%%%%%%%%%%%%
#### Section 5 ####
#'@description fit segmented regression with "segmented" package

# glm fitter
glm.winter=glmmTMB(presence ~  hsm , 
                   data = db.pres.sb,
                   family=binomial(link = "logit")) ## format not working

glm.winter.bis=glm(presence ~  hsm + species.binomial , 
                   data = db.pres.sb,
                   family=binomial(link = "logit"))
# segmented regressions
glm.segmented=segmented(obj=glm.winter.bis,seg.Z=~hsm,npsi=1) 


#%%%%%%%%%%%%%%%%%%
#### Section 6 ####
#'@description fit segmented regression with "mcp" package

## models on FSM
# MCP model for one species, one safety margins
library(mcp)
model = list(
  presence ~ fsm,  # changing rate
  ~ 0 
)
prior=list(
  cp_1="dnorm(-0.9199044,1)"
)

data.mcp=data.frame(presence= db.pres.sb$presence,
                    fsm=db.pres.sb$fsmwinter)
fit.mcp.raw= mcp(model, family = bernoulli(),data=data.mcp)
fit.mcp = mcp(model, family = bernoulli(),data=data.mcp,sample = "both", adapt = 10000, prior = prior)
summary(fit.mcp)


#MCP model for several species
model = list(
  presence ~ fsm,  # changing rate
  1+ (1|species.binomial)~ 0 
)
prior=list(
  cp_1="dnorm(-0.9199044,1)"
)

data.mcp=data.frame(presence= db.pres.sb$presence,
                    fsm=db.pres.sb$fsmwinter,
                    species.binomial=as.factor(db.pres.sb$species.binomial))
fit.mcp.raw= mcp(model, family = bernoulli(),data=data.mcp)
fit.mcp = mcp(model, family = bernoulli(),data=data.mcp,sample = "both", adapt = 10000, prior = prior)
summary(fit.mcp)


#%%%%%%%%%%%%%%%%%%
#### Section 7 ####
#'@description fit segmented regression with "stan" package

# df for species code
species <- data.frame(species=unique(factor(db.pres$species.binomial)),
                   species.nb=as.character(as.numeric(unique(factor(db.pres$species.binomial)))))

# stan model, S species, FSM only
data.list <- list(N=dim(db.pres)[1],
                  S=nlevels(factor(db.pres$species.binomial)),
                  presence=db.pres$presence,
                  #hsm=as.numeric(db.pres.sb$hsm),
                  species=as.numeric(factor(db.pres$species.binomial)),
                  fsm=as.numeric(db.pres$fsmwinter),
                  prior_breakpoint=-fsm.mu/fsm.sd,
                  NULL)
fit <- stan(file = "glm_seg_hsm.stan",
            data=data.list,
            iter=1000,
            chains=2)
sum.fit=summary(fit)

##plot results
breakpoint=sum.fit$summary["breakpoint","mean"]
b_fsm=sum.fit$summary["b_fsm","mean"]
plateau=sum.fit$summary["plateau","mean"]
x.fsm=seq(min(data.list$fsm),max(data.list$fsm),len=1000)
p.fsm=inv.logit((x.fsm>breakpoint)*plateau +
  (x.fsm<breakpoint)*(b_fsm * (x.fsm - breakpoint) + plateau))
fit.mean=data.frame(x.fsm,p.fsm)
db.pres %>% 
  ggplot()+
  geom_point(aes(fsmwinter,presence),color="red")+
  geom_line(data=fit.mean,aes(x.fsm,p.fsm))+
  geom_vline(xintercept = data.list$prior_breakpoint)


# stan model, S species, FSM & HSM
data.list <- list(N=dim(db.pres)[1],
                  S=nlevels(factor(db.pres$species.binomial)),
                  presence=db.pres$presence,
                  #hsm=as.numeric(db.pres.sb$hsm),
                  species=as.numeric(factor(db.pres$species.binomial)),
                  fsm=db.pres$fsmwinter,
                  hsm=db.pres$hsm,
                  prior_threshold_fsm=-fsm.mu/fsm.sd,
                  prior_threshold_hsm=-hsm.mu/hsm.sd,
                  NULL)
fit.2sm <- stan(file = "glm_seg_all.stan",
            data=data.list,
            iter=2000,
            chains=3,
            core=3)
sum.fit=summary(fit.2sm)
save(fit.2sm,file="fit_4sp_all.RData")

#plot results
load("fit_4sp_all.RData")
fit.sample=as.data.frame(fit.2sm)
plot.fit <-fit.sample %>% 
  select(-matches("_pred"),-matches("_app")) %>% 
  select(-"lp__") %>% 
  sample_n(800) %>% 
  # slice(1001:2000) %>% 
  mutate(it=row_number()) %>% 
  crossing(fsm=seq(min(db.pres$fsmwinter),max(db.pres$fsmwinter),by=0.01)) %>% 
  mutate(hsm=median(db.pres$fsmwinter)) %>% 
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>% 
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>% 
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  mutate(fsm_app=case_when(fsm<threshold_fsm~0,TRUE~1),
         hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
         pred= inv.logit((1-fsm_app)*(b_fsm * (fsm - threshold_fsm) + plateau_fsm) 
                + fsm_app * plateau_fsm
                + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm) 
                + hsm_app * plateau_hsm))
  
  
plot.fit %>%
  left_join(species,by=c("species1"="species.nb")) %>% 
  group_by(species,fsm) %>% 
  summarise(mean=mean(pred),
            q95=quantile(pred,probs=0.95),
            q05=quantile(pred,probs=0.05)) %>% 
  ggplot()+
  geom_line(aes(fsm,mean))+
  geom_ribbon(aes(fsm,ymin=q05,ymax=q95),alpha=0.2)+
  # geom_point(data=db.pres,aes(fsmwinter,presence))+
  facet_wrap(~species)



# make predictions on a subdataset
set.seed(45844)
sample_length=9500
slice_post=1500
auc=matrix(rep(0,slice_post*10), 
           nrow = 10 )
fit <- fit.sample[,!grepl("presence_pred",colnames(fit.sample))] %>%
  select(-"lp__") %>% 
  slice(1:slice_post) %>% 
  mutate(it=row_number()) %>% 
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>% 
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>% 
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2)
for (i in 1:10){
  db.pres.test <- db.clim %>%
    filter(species.binomial%in%c("Fagus sylvatica","Pinus sylvestris","Fagus sylvatica","Fraxinus excelsior","Quercus ilex")) %>% #"
    #filter(hsm>hsm.q01) %>% 
    sample_n(sample_length) %>% 
    select(presence,species.binomial,fsmwinter,hsm) %>% 
    left_join(species,by=c("species.binomial"="species")) %>% 
    mutate(pred_n=row_number(),
           hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
           fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1]) 
  pred.test <-  fit %>% 
    crossing(db.pres.test) %>% 
    filter(species1==species.nb) %>% 
    mutate(fsm_app=case_when(fsmwinter<threshold_fsm~0,TRUE~1),
           hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
           prob=inv.logit(
                           (1-fsm_app)*(b_fsm * (fsmwinter - threshold_fsm) + plateau_fsm) 
                           + fsm_app * plateau_fsm
                           + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm) 
                           + hsm_app * plateau_hsm
                           ),
           pred=as.numeric(purrr::rbernoulli(dim(.)[1],prob))) %>% 
    select("presence","pred_n","pred") %>%
    group_by(pred_n) %>%
    mutate(n=row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = n,values_from=pred) %>%
    select(-pred_n)

  pred.test=as.data.frame(pred.test)
  # auc <- numeric(dim(pred.test)[2])
  for (n in 2:dim(pred.test)[2]){
    auc[i,n-1] <- pROC::auc(pred.test[,1],pred.test[,n]
    )
  }
}
as.data.frame(t(auc)) %>% pivot_longer(cols=colnames(.)) %>% ggplot(aes(value,color=name))+geom_density()

# 2 SM, 4SP, varying slopes
data.list <- list(N=dim(db.pres)[1],
                  S=nlevels(factor(db.pres$species.binomial)),
                  presence=db.pres$presence,
                  #hsm=as.numeric(db.pres.sb$hsm),
                  species=as.numeric(factor(db.pres$species.binomial)),
                  fsm=db.pres$fsmwinter,
                  hsm=db.pres$hsm,
                  prior_threshold_fsm=-fsm.mu/fsm.sd,
                  prior_threshold_hsm=-hsm.mu/hsm.sd,
                  NULL)
fit.2sm.slope <- stan(file = "glm_seg_allslope.stan",
                data=data.list,
                iter=2000,
                chains=3)
save(fit.2sm.slope,file="fit_4sp_allslope.RData")

#plot results
load("fit_4sp_allslope.RData")

fit.sample=as.data.frame(fit.2sm.slope)
plot.fit <- fit.sample%>% 
  select(-matches("_pred")) %>% 
  select(-"lp__") %>% 
  sample_n(400) %>% 
  # slice(1001:2000) %>% 
  mutate(it=row_number()) %>% 
  crossing(fsm=seq(min(db.pres$fsmwinter),max(db.pres$fsmwinter),by=0.01)) %>% 
  mutate(hsm=median(db.pres$fsmwinter)) %>% 
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>% 
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>% 
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  pivot_longer(cols=matches("b_fsm"),names_to = "species2",values_to="b_fsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  pivot_longer(cols=matches("b_hsm"),names_to = "species2",values_to = "b_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  mutate(fsm_app=case_when(fsm<threshold_fsm~0,TRUE~1),
         hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
         pred= inv.logit((1-fsm_app)*(b_fsm * (fsm - threshold_fsm) + plateau_fsm) 
                         + fsm_app * plateau_fsm
                         + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm) 
                         + hsm_app * plateau_hsm))


plot.fit %>%
  left_join(species,by=c("species1"="species.nb")) %>% 
  group_by(species,fsm) %>% 
  summarise(mean=mean(pred),
            q95=quantile(pred,probs=0.95),
            q05=quantile(pred,probs=0.05)) %>% 
  ggplot()+
  geom_line(aes(fsm,mean))+
  geom_ribbon(aes(fsm,ymin=q05,ymax=q95),alpha=0.2)+
  # geom_point(data=db.pres,aes(fsmwinter,presence))+
  facet_wrap(~species)

# compute auc
fit.sample=as.data.frame(fit.2sm.slope) %>% 
  select(-matches("_app")) %>% 
  select(-"lp__") %>% 
  # slice(1001:2000) %>% 
  mutate(it=row_number()) %>% 
  crossing(fsm=seq(min(db.pres$fsmwinter),max(db.pres$fsmwinter),by=0.01)) %>% 
  mutate(hsm=median(db.pres$fsmwinter)) %>% 
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>% 
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>% 
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  pivot_longer(cols=matches("b_fsm"),names_to = "species2",values_to="b_fsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  pivot_longer(cols=matches("b_hsm"),names_to = "species2",values_to = "b_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2)

# make predictions on a subdataset
set.seed(45844)
sample_length=5000
slice_post=1000
auc=matrix(rep(0,slice_post*10), 
           nrow = 10 )
fit <- fit.sample[,!grepl("_app",colnames(fit.sample))] %>%
  select(-"lp__") %>% 
  slice(1:slice_post) %>% 
  mutate(it=row_number()) %>% 
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>% 
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>% 
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  pivot_longer(cols=matches("b_fsm"),names_to = "species2",values_to="b_fsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2) %>% 
  pivot_longer(cols=matches("b_hsm"),names_to = "species2",values_to = "b_hsm") %>% 
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>% 
  filter(species1==species2) %>% 
  select(-species2)
for (i in 1:10){
  db.pres.test <- db.clim %>%
    filter(species.binomial%in%c("Fagus sylvatica","Pinus sylvestris","Fraxinus excelsior","Quercus ilex")) %>% #"
    #filter(hsm>hsm.q01) %>% 
    sample_n(sample_length) %>% 
    select(presence,species.binomial,fsmwinter,hsm) %>% 
    left_join(species,by=c("species.binomial"="species")) %>% 
    mutate(pred_n=row_number(),
           hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
           fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1]) 
  pred.test <-  fit %>% 
    crossing(db.pres.test) %>% 
    filter(species1==species.nb) %>% 
    mutate(fsm_app=case_when(fsm<threshold_fsm~0,TRUE~1),
           hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
           prob=inv.logit(
             (1-fsm_app)*(b_fsm * (fsm - threshold_fsm) + plateau_fsm) 
             + fsm_app * plateau_fsm
             + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm) 
             + hsm_app * plateau_hsm
           ),
           pred=as.numeric(purrr::rbernoulli(dim(.)[1],prob))) %>% 
    select("presence","pred_n","pred") %>%
    group_by(pred_n) %>%
    mutate(n=row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = n,values_from=pred) %>%
    select(-pred_n)
  
  pred.test=as.data.frame(pred.test)
  # auc <- numeric(dim(pred.test)[2])
  for (n in 2:dim(pred.test)[2]){
    auc[i,n-1] <- pROC::auc(pred.test[,1],pred.test[,n]
    )
  }
}
as.data.frame(t(auc)) %>% pivot_longer(cols=colnames(.)) %>% ggplot(aes(value,color=name))+geom_density()





# Model with marginalized changing point
data.list <- list(N=dim(db.pres)[1],
                  presence=db.pres$presence,
                  fsm=db.pres$fsmwinter,
                  breakpoint_n=15,
                  breakpoints=seq(from=min(db.pres$fsmwinter),to=max(db.pres$fsmwinter),length.out=15))
fit.fsm.marg <- stan(file = "glm_seg_hsm_marginalized.stan",
                      data=data.list,
                      iter=1000,
                      chains=2)


#%%%%%%%%%%%%%%%%%%
#### Section 8 ####
#'@description fit log regression with "stan" package 

# Model with logistic curves fit, 1 species, both SM
db.pres.1sp <- db.clim %>%
  filter(species.binomial=="Betula pendula") %>%
  mutate(presence=as.factor(presence))
  sample_frac(0.05)
hsm.mu=mean(db.pres.1sp$hsm,na.rm=TRUE)
hsm.sd=sd(db.pres.1sp$hsm,na.rm=TRUE)
fsm.mu=mean(db.pres.1sp$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres.1sp$fsmwinter,na.rm=TRUE)
db.pres.1sp <- db.pres.1sp %>% 
  mutate(hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
         fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1]) 

data.list <- list(N=dim(db.pres.1sp)[1],
                  presence=db.pres.1sp$presence,
                  fsm=as.numeric(db.pres.1sp$fsmwinter),
                  hsm=as.numeric(db.pres.1sp$hsm),
                  prior_t_fsm=-fsm.mu/fsm.sd,
                  prior_t_hsm=-hsm.mu/hsm.sd,
                  NULL)
fit.1.sp <- stan(file = "glm_log_1sp.stan",
            data=data.list,
            iter=1000,
            chains=2,
            include=FALSE,
            pars=c("proba","K_vect"))
sum.fit1sp=summary(fit)$summary
data.frame(sm,t(sum.fit1sp[,1])) %>%
  crossing(hsm=sm) %>% 
  mutate(proba=K_int/((1 + exp(-r_fsm * (sm - t_fsm)))*(1 + exp(-r_hsm * (hsm - t_hsm))))) %>% 
  filter(hsm %in% c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) %>% 
  ggplot(aes(sm,proba,color=hsm))+
  geom_point()



# Model with logistic curves fit, 2 species, both SM
db.pres.2sp <- db.clim %>%
  filter(species.binomial%in%c("Fagus sylvatica","Quercus ilex")) %>%
  sample_frac(0.05)
hsm.mu=mean(db.pres.2sp$hsm,na.rm=TRUE)
hsm.sd=sd(db.pres.2sp$hsm,na.rm=TRUE)
fsm.mu=mean(db.pres.2sp$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres.2sp$fsmwinter,na.rm=TRUE)
db.pres.2sp <- db.pres.2sp %>% 
  mutate(hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
         fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1]) 

data.list <- list(N=dim(db.pres.2sp)[1],
                  S=nlevels(as.factor(db.pres.2sp$species.binomial)),
                  presence=db.pres.2sp$presence,
                  species=as.numeric(as.factor(db.pres.2sp$species.binomial)),
                  fsm=as.numeric(db.pres.2sp$fsmwinter),
                  hsm=as.numeric(db.pres.2sp$hsm),
                  prior_t_fsm=-fsm.mu/fsm.sd,
                  prior_t_hsm=-hsm.mu/hsm.sd,
                  NULL)
fit.2sp <- stan(file = "glm_log_all.stan",
                 data=data.list,
                 iter=1000,
                 chains=2,
                 include=FALSE,
                 pars=c("proba","K_vect"))


# Model with logistic curves fit, 4 species, both SM
db.pres.4sp <- db.clim %>%
  filter(species.binomial%in%c("Fagus sylvatica","Pinus sylvestris","Fraxinus excelsior","Quercus ilex")) %>% #"
  sample_frac(0.1)
hsm.mu=mean(db.pres.4sp$hsm,na.rm=TRUE)
hsm.sd=sd(db.pres.4sp$hsm,na.rm=TRUE)
fsm.mu=mean(db.pres.4sp$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres.4sp$fsmwinter,na.rm=TRUE)
db.pres.4sp <- db.pres.4sp %>% 
  mutate(hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
         fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1]) 

data.list <- list(N=dim(db.pres.4sp)[1],
                  S=nlevels(as.factor(db.pres.4sp$species.binomial)),
                  presence=db.pres.4sp$presence,
                  species=as.numeric(as.factor(db.pres.4sp$species.binomial)),
                  fsm=as.numeric(db.pres.4sp$fsmwinter),
                  hsm=as.numeric(db.pres.4sp$hsm),
                  prior_t_fsm=-fsm.mu/fsm.sd,
                  prior_t_hsm=-hsm.mu/hsm.sd,
                  NULL)
fit.4sp <- stan(file = "glm_log_all.stan",
                data=data.list,
                iter=2000,
                chains=2,
                core=2,
                include=FALSE,
                pars=c("proba","K_vect"))

#%%%%%%%%%%%%%%%%%%
#### Section 9 ####
#'@description Comparison of time of time to fit a model several species using
#'logistic model and different types of scaling

# scale all species together; no unskewing
sp=c("Quercus robur","Fagus sylvatica",
  "Betula pendula","Picea abies")
set.seed(4554)
db.pres <- db.clim %>%
  filter(psi>(-10000)) %>% 
  filter(species.binomial%in%sp) %>% 
  sample_frac(0.03) 


# db.pres$hsm <- -log(max(db.pres$hsm+1)-db.pres$hsm)
hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)

# scale and center al sp together 
# db.pres <- db.pres %>% 
#   mutate(hsm=(hsm-hsm.mu)/hsm.sd,
#          fsm=(fsmwinter-fsm.mu)/fsm.sd)

# scale, center, unskew by sp
db.pres <- db.pres %>%
  group_by(species.binomial) %>%
  mutate(hsm=scale(hsm,center=TRUE,scale=TRUE),
         fsm=scale(fsmwinter,center=TRUE,scale=TRUE))

data.list <- list(N=dim(db.pres)[1],
                  S=nlevels(as.factor(db.pres$species.binomial)),
                  presence=db.pres$presence,
                  species=as.numeric(as.factor(db.pres$species.binomial)),
                  fsm=as.numeric(db.pres$fsm),
                  hsm=as.numeric(db.pres$hsm),
                  prior_t_fsm=-hsm.mu/hsm.sd,
                  prior_t_hsm=-fsm.mu/fsm.sd,
                  NULL)

tic()
fit <- stan(file = "glm_log_all.stan",
            data=data.list,
            iter=1000,
            chains=2,
            core=2,
            include=FALSE,
            pars=c("proba","K_vect"))
toc()



#%%%%%%%%%%%%%%%%%%
#### Section 10 ####
#'@description Comparison of time of time to fit a model on one species using
#'logistic model and threshold model
# scale all species together; no unskewing
set.seed(4554)
db.pres <- db.clim %>%
  filter(psi>(-10000)) %>% 
  filter(species.binomial=="Fagus sylvatica") %>% 
  sample_frac(0.03) 

hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)

db.pres <- db.pres %>%
  mutate(hsm=scale(hsm,center=TRUE,scale=TRUE),
         fsm=scale(fsmwinter,center=TRUE,scale=TRUE))

data.list <- list(N=dim(db.pres)[1],
                  presence=db.pres$presence,
                  fsm=as.numeric(db.pres$fsm),
                  hsm=as.numeric(db.pres$hsm),
                  prior_t_fsm=-hsm.mu/hsm.sd,
                  prior_t_hsm=-fsm.mu/fsm.sd,
                  NULL)

tic()
fit <- stan(file = "glm_log_1sp.stan",
            data=data.list,
            iter=1000,
            chains=2,
            core=2,
            include=FALSE,
            pars=c("proba","K_vect"))
toc()


tic()
fit.1 <- stan(file = "glm_seg_1sp.stan",
            data=data.list,
            iter=1000,
            chains=2,
            core=2,
            include=FALSE,
            pars=c("proba","K_vect"))
toc()

#%%%%%%%%%%%%%%%%%%
#### Section 11 ####
#'@description Workflow to fit log regression with "stan" package per species

set.seed(4554)
# species selection for tests
species.select = c("Betula pendula",
       "Pinus sylvestris",
       "Fagus sylvatica",
       "Fraxinus excelsior")
for (sp in species.select){ #unique(db.clim$species.binomial)
  print(sp)
  db.pres <- db.clim %>%
    filter(psi>(-10000)) %>% #remove very low value of psi
    filter(species.binomial==sp) %>% 
    mutate(hsm=hsm/1000,
           fsm=fsmwinter) %>% 
    sample_frac(0.05)  
    # add_row(species.binomial="zero_like",
    #         fsmwinter=0,
    #         hsm=0)

  # unskew
  # hsm.max=max(db.pres$hsm,na.rm = TRUE)
  # db.pres$hsmsave <- db.pres$hsm
  # db.pres$hsm <- -log(max(db.pres$hsm+1)-db.pres$hsm)
  
  # scale and center  
  hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
  hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
  fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
  fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)
  db.pres.sc <- db.pres %>% 
    mutate(hsm=(hsm-hsm.mu)/hsm.sd,
           fsm=(fsmwinter-fsm.mu)/fsm.sd)
  # hsm.zero=db.pres[db.pres$species.binomial=="zero_like","hsm"][[1]]
  # fsm.zero=db.pres[db.pres$species.binomial=="zero_like","fsm"][[1]]
  # db.pres=db.pres[db.pres$species.binomial!="zero_like",]
  # db.pres %>% ggplot(aes(hsm))+geom_density()+geom_vline(xintercept = hsm.zero)
  
  
  data.list <- list(N=dim(db.pres)[1],
                    presence=db.pres$presence,
                    fsm=as.numeric(db.pres$fsm),
                    hsm=as.numeric(db.pres$hsm),
                    prior_K=df.traits[df.traits$species.binomial==sp,
                                     "prevalence"],
                    NULL)
  data.list.sc <- list(N=dim(db.pres.sc)[1],
                    presence=db.pres.sc$presence,
                    fsm=as.numeric(db.pres.sc$fsm),
                    hsm=as.numeric(db.pres.sc$hsm),
                    prior_K=df.traits[df.traits$species.binomial==sp,
                                      "prevalence"],
                    NULL)
  # model I 
  fitI <- stan(file = "glm_log_1sp_I.stan",
              data=data.list,
              iter=1000,
              chains=2,
              core=2,
              include=FALSE,
              pars=c("proba","K_vect"))
  save(fitI,file=paste0("fit_mod3/",sp,"_fitI.RData"))
  #model I scaled
  fitI.sc <- stan(file = "glm_log_1sp_I.stan",
               data=data.list.sc,
               iter=1000,
               chains=2,
               core=2,
               include=FALSE,
               pars=c("proba","K_vect"))
  save(fitI.sc,file=paste0("fit_mod3/",sp,"_fitIsc.RData"))
  # model II
  fitI <- stan(file = "glm_log_1sp_II.stan",
               data=data.list,
               iter=1000,
               chains=2,
               core=2,
               include=FALSE,
               pars=c("proba","K_vect"))
  save(fitI,file=paste0("fit_mod3/",sp,"_fitII.RData"))
  # model III
  fitI <- stan(file = "glm_log_1sp_III.stan",
               data=data.list,
               iter=1000,
               chains=2,
               core=2,
               include=FALSE,
               pars=c("proba","K_vect"))
  save(fitI,file=paste0("fit_mod3/",sp,"_fitIII.RData"))
  # model IV
  fitI <- stan(file = "glm_log_1sp_IV.stan",
               data=data.list,
               iter=1000,
               chains=2,
               core=2,
               include=FALSE,
               pars=c("proba","K_vect"))
  save(fitI,file=paste0("fit_mod3/",sp,"_fitIV.RData"))
  
  list.param=list(fsm.mu=fsm.mu,
                  fsm.sd=fsm.sd,
                  hsm.mu=hsm.mu,
                  hsm.sd=hsm.sd,
                  NULL)
  save(list.param,file=paste0("fit_mod3/",sp,"_param.RData"))
}

#%%%%%%%%%%%%%%%%%%
#### Section 11.bis ####
#'@description Workflow to fit log regression with "stan" package per species

set.seed(4554)
for (sp in unique(db.clim$species.binomial)){
  print(sp)
  db.pres <- db.clim %>%
    filter(psi>(-10000)) %>% #remove very low value of psi
    filter(species.binomial==sp) %>% 
    mutate(hsm=hsm/1000,
           fsm=fsmwinter) %>% 
    sample_frac(0.6)  

  data.list <- list(N=dim(db.pres)[1],
                    presence=db.pres$presence,
                    fsm=as.numeric(db.pres$fsm),
                    hsm=as.numeric(db.pres$hsm),
                    prior_K=df.traits[df.traits$species.binomial==sp,
                                      "prevalence"],
                    NULL)
  # model III
  fit <- stan(file = "glm_log_1sp_III.stan",
               data=data.list,
               iter=1000,
               chains=2,
               core=2,
               include=FALSE,
               pars=c("proba","K_vect"))
  save(fit,file=paste0("fit_mod4/",sp,"_fitIII.RData"))

}




#%%%%%%%%%%%%%%%
#### Section 12 ####
#'@description format outputs from models fits
#'

mod.files=list.files("fit_mod/")
species.list=unique(str_split_i(mod.files,"_",1))

df.output=df.traits %>% 
  filter(species.binomial %in% species.list) %>% 
  select(species.binomial,Group,p.trait,PX.mu,PX.sd,LT50.mean,LT50.sd) %>% 
  mutate(K_int=NA,
         r_fsm=NA,
         r_hsm=NA,
         t_hsm=NA,
         t_fsm=NA,
         hsm.mu=NA,
         hsm.sd=NA,
         fsm.mu=NA,
         fsm.sd=NA)
for (sp in species.list){
  print(sp)
  load(paste0("fit_mod/",sp,"_logfit.RData"))
  load(paste0("fit_mod/",sp,"_param.RData"))
  summary.fit=as.data.frame(summary(fit)$summary)
  if( summary.fit["lp__","Rhat"]<1.2){
    for (par in c("K_int","r_fsm","r_hsm","t_hsm","t_fsm")){
      print(par)
      df.output[df.output$species.binomial==sp,par]=summary.fit[par,"mean"]
    }
    for (par in c("hsm.mu","hsm.sd","fsm.mu","fsm.sd")){
      print(par)
      df.output[df.output$species.binomial==sp,par]=list.param[par][[1]]
    }
  }
}

df.output <- df.output %>% 
  mutate(r_hsm_app=r_hsm/hsm.sd,
         r_fsm_app=r_fsm/fsm.sd,
         t_hsm_app=hsm.mu+t_hsm*hsm.sd,
         t_fsm_app=fsm.mu+t_fsm*fsm.sd)

df.output %>% filter(is.na(K_int))

df.output %>% 
  filter(!is.na(K_int)) %>% 
  ggplot(aes(r_hsm_app,r_fsm_app,color=species.binomial))+
  geom_point(size=5)

df.output %>% 
  filter(!is.na(K_int)) %>% 
  # select(species.binomial,PX.mu,LT50.mean,t_hsm_app)
  ggplot(aes(t_fsm_app,t_hsm_app,color=species.binomial))+
  geom_point(size=5)

df.output %>% 
  filter(!is.na(K_int)) %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm_app,t_fsm_app) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm_app","t_hsm_app"),names_to="threshold",values_to="threshold_val") %>% 
  filter((traits=="PX.mu"&threshold=="t_hsm_app")|
           (traits=="LT50.mean"&threshold=="t_fsm_app")) %>% 
  ggplot(aes(traits_val,threshold_val,color=species.binomial,label=species.binomial))+
  geom_point(size=3)+
  geom_text(hjust=0.4, vjust=-0.5)+
  facet_wrap(traits~threshold,scales="free")

hsm.95=quantile(db.clim$hsm,prob=0.95)[[1]]
hsm.05=quantile(db.clim$hsm,prob=0.05)[[1]]
fsm.95=quantile(db.clim$fsmwinter,prob=0.95)[[1]]
fsm.05=quantile(db.clim$fsmwinter,prob=0.05)[[1]]
df.output %>% 
  filter(!is.na(K_int)) %>% 
  crossing(hsm=seq(-100,
                   100,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("low","indif","high"),
                      fsm=c(fsm.05,0,fsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm_app*(fsm-t_fsm_app)))*
                       (1+exp(-r_hsm_app*(hsm-t_hsm_app))))) %>% 
  ggplot(aes(hsm,pred,color=species.binomial))+
  geom_line(size=1)+
  facet_wrap(~fsm.type,scales="free")

df.output %>% 
  filter(!is.na(K_int)) %>% 
  crossing(fsm=seq(fsm.05,
                   fsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(hsm.type=c("low","indif","high"),
                      hsm=c(hsm.05,0,hsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm_app*(fsm-t_fsm_app)))*
                       (1+exp(-r_hsm_app*(hsm-t_hsm_app))))) %>% 
  ggplot(aes(fsm,pred,color=species.binomial))+
  geom_line(size=1)+
  facet_wrap(~hsm.type,scales="free")


#%%%%%%%%%%%%%%%
#### Section 13 ####
#'@description format outputs from models fits coomparison for several species
#'

for (sp in species.list){
  print(sp)
  load(paste0(folder,sp,"_","fitIII.RData"))
  print(as.data.frame(summary(fit)$summary))
  
}

folder="fit_mod4/"
mod.files=list.files(folder)[grepl(pattern = "fitI", x = list.files(folder))]
species.list=unique(str_split_i(mod.files,"_",1))
model.list=unique(str_split_i(mod.files,"_",2))
df.output=df.traits %>% 
  filter(species.binomial %in% species.list) %>% 
  select(species.binomial,Group,p.trait,PX.mu,PX.sd,LT50.mean,LT50.sd,prevalence) %>% 
  mutate(K_int=NA,
         r_fsm=NA,
         r_hsm=NA,
         t_hsm=NA,
         t_fsm=NA,
         hsm.mu=NA,
         hsm.sd=NA,
         fsm.mu=NA,
         fsm.sd=NA,
         Rhat=NA,
         lp__=NA) %>% 
  crossing(mod=model.list)


for (sp in species.list){
  print(sp)
  # load(paste0("fit_mod3/",sp,"_param.RData"))
  # for (par in c("hsm.mu","hsm.sd","fsm.mu","fsm.sd")){
  #   print(par)
  #   df.output[df.output$species.binomial==sp,par]=list.param[par][[1]]
  # }
  for(mod in model.list){
    if (mod=="fitIsc.RData"){
      load(paste0(folder,sp,"_",mod))
      summary.fit=as.data.frame(summary(fitI.sc)$summary)
      for (par in c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","lp__")){
        print(par)
        df.output[df.output$species.binomial==sp
                  & df.output$mod==mod,par]=summary.fit[par,"mean"]
      }
      df.output[df.output$species.binomial==sp
                &df.output$mod==mod,"Rhat"] = max(summary.fit[,"Rhat"])
    } else {
      load(paste0(folder,sp,"_",mod))
      summary.fit=as.data.frame(summary(fit)$summary)#fit
      for (par in c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","lp__")){
        print(par)
        df.output[df.output$species.binomial==sp
                  & df.output$mod==mod,par]=summary.fit[par,"mean"]
      }
      df.output[df.output$species.binomial==sp
                &df.output$mod==mod,"Rhat"] = max(summary.fit[,"Rhat"])
    }
  }
}

df.output.sub <- df.output %>% 
  filter(mod=="fitIsc.RData") %>% 
  mutate(r_hsm_app=r_hsm/hsm.sd,
         r_fsm_app=r_fsm/fsm.sd,
         t_hsm_app=hsm.mu+t_hsm*hsm.sd,
         t_fsm_app=fsm.mu+t_fsm*fsm.sd) %>% 
  mutate(r_hsm=r_hsm_app,
         r_fsm=r_fsm_app,
         t_hsm=t_hsm_app,
         t_fsm=t_fsm_app) %>% 
  select(-c("r_hsm_app","r_fsm_app","t_hsm_app","t_fsm_app"))

df.output <- df.output %>% 
  filter(!mod=="fitIsc.RData") %>% 
  bind_rows(df.output.sub) %>% 
  mutate(bic=2*5-2*lp__)


df.output %>% 
  filter(Rhat<1.1) %>% 
  ggplot(aes(r_hsm,r_fsm,size=-log(Rhat),color=species.binomial,shape=mod))+
  geom_point()+
  scale_y_log10()



df.output %>% 
  filter(Rhat<1.1) %>% 
  ggplot(aes(t_fsm,t_hsm,color=species.binomial,shape=mod))+
  geom_point(size=5)

df.output %>% 
  filter(Rhat<1.1) %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm,t_fsm,mod) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to="threshold",values_to="threshold_val") %>% 
  filter((traits=="PX.mu"&threshold=="t_hsm")|
           (traits=="LT50.mean"&threshold=="t_fsm")) %>% 
  ggplot(aes(traits_val,threshold_val,color=species.binomial,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  geom_abline(slope=-1)+
  facet_wrap(traits~threshold,scales="free")

df.output %>% 
  select(species.binomial,PX.mu,LT50.mean,r_hsm,r_fsm,mod) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("r_fsm","r_hsm"),names_to="slope",values_to="slope_val") %>% 
  filter((traits=="PX.mu"&slope=="r_hsm")|
           (traits=="LT50.mean"&slope=="r_fsm")) %>% 
  filter(slope_val<20) %>% 
  ggplot(aes(traits_val,slope_val,color=species.binomial,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  facet_wrap(traits~slope,scales="free")


hsm.95=quantile(db.clim$hsm,prob=0.95)[[1]]/1000
hsm.05=quantile(db.clim$hsm,prob=0.05)[[1]]/1000
fsm.95=quantile(db.clim$fsmwinter,prob=0.95)[[1]]
fsm.05=quantile(db.clim$fsmwinter,prob=0.05)[[1]]
df.output %>% 
  crossing(hsm=seq(hsm.05,
                   hsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("low","indif","high"),
                      fsm=c(fsm.05,0,fsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(hsm,pred,color=mod))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=prevalence))+
  facet_wrap(fsm.type~species.binomial,scales="free", nrow = 3)

df.output %>% 
  crossing(fsm=seq(fsm.05,
                   fsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(hsm.type=c("low","indif","high"),
                      hsm=c(hsm.05,0,hsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(fsm,pred,color=mod))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=prevalence))+
  facet_wrap(hsm.type~species.binomial,scales="free", nrow = 3)

df.output %>% 
  crossing(hsm=seq(hsm.05,
                   hsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("low","indif","high"),
                      fsm=c(fsm.05,0,fsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(hsm,pred,color=species.binomial))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=prevalence,color=species.binomial))+
  facet_wrap(fsm.type~mod, nrow = 3)


df.output %>% 
  filter(Rhat<1.1) %>% 
  crossing(hsm=seq(hsm.05,
                   hsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("low","indif","high"),
                      fsm=c(fsm.05,0,fsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(hsm,pred,color=fsm.type))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=prevalence))+
  facet_wrap(~species.binomial,scales="free", nrow = 3)
df.output %>% 
  filter(Rhat<1.1) %>% 
  crossing(fsm=seq(fsm.05,
                   fsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(hsm.type=c("low","indif","high"),
                      hsm=c(hsm.05,0,hsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(fsm,pred,color=hsm.type))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=prevalence))+
  facet_wrap(~species.binomial,scales="free", nrow = 3)

#%%%%%%%%%%%%%%%
#### Section 14 ####
#'@description dig into non converging species
#'
#_####
sp = "Betula pendula"
# sp = "Pinus nigra"
# sp = "Picea abies"
# sp = "Quercus ilex"
set.seed(4554)
print(sp)
db.pres <- db.clim %>%
  filter(psi>(-10000)) %>% #remove very low value of psi
  filter(species.binomial==sp) %>% 
  mutate(hsm=hsm/1000,
         fsm=fsmwinter) %>% 
  sample_frac(0.6)  

data.list <- list(N=dim(db.pres)[1],
                  presence=db.pres$presence,
                  fsm=as.numeric(db.pres$fsm),
                  hsm=as.numeric(db.pres$hsm),
                  prior_K=df.traits[df.traits$species.binomial==sp,
                                    "prevalence"],
                  NULL)
# model III
fit.bepeIII <- stan(file = "glm_log_1sp_III.stan",
            data=data.list,
            iter=1000,
            chains=2,
            core=2,
            include=FALSE,
            pars=c("proba","K_vect"))
fit.bepeIV <- stan(file = "glm_log_1sp_IV.stan",
            data=data.list,
            iter=1000,
            chains=2,
            core=2,
            include=FALSE,
            pars=c("proba","K_vect"))

#_####
sp = "Picea abies"
set.seed(4554)
print(sp)
db.pres <- db.clim %>%
  filter(psi>(-10000)) %>% #remove very low value of psi
  filter(species.binomial==sp) %>% 
  mutate(hsm=hsm/1000,
         fsm=fsmwinter) %>% 
  sample_frac(0.6)  

data.list <- list(N=dim(db.pres)[1],
                  presence=db.pres$presence,
                  fsm=as.numeric(db.pres$fsm),
                  hsm=as.numeric(db.pres$hsm),
                  prior_K=df.traits[df.traits$species.binomial==sp,
                                    "prevalence"],
                  NULL)
# model III
fit.piabIII <- stan(file = "glm_log_1sp_III.stan",
            data=data.list,
            iter=1000,
            chains=2,
            core=2,
            include=FALSE,
            pars=c("proba","K_vect"))
fit.piabIV <- stan(file = "glm_log_1sp_IV.stan",
               data=data.list,
               iter=1000,
               chains=2,
               core=2,
               include=FALSE,
               pars=c("proba","K_vect"))

#_#### 
sp = "Pinus nigra"
set.seed(4554)
print(sp)
db.pres <- db.clim %>%
  filter(psi>(-10000)) %>% #remove very low value of psi
  filter(species.binomial==sp) %>% 
  mutate(hsm=hsm/1000,
         fsm=fsmwinter) %>% 
  sample_frac(0.6)  

data.list <- list(N=dim(db.pres)[1],
                  presence=db.pres$presence,
                  fsm=as.numeric(db.pres$fsm),
                  hsm=as.numeric(db.pres$hsm),
                  prior_K=df.traits[df.traits$species.binomial==sp,
                                    "prevalence"],
                  NULL)
# model III
fit.piniIII <- stan(file = "glm_log_1sp_III.stan",
                    data=data.list,
                    iter=1000,
                    chains=2,
                    core=2,
                    include=FALSE,
                    pars=c("proba","K_vect"))
fit.piniIV <- stan(file = "glm_log_1sp_IV.stan",
                   data=data.list,
                   iter=1000,
                   chains=2,
                   core=2,
                   include=FALSE,
                   pars=c("proba","K_vect"))

#_####
sp = "Quercus ilex"
set.seed(4554)
print(sp)
db.pres <- db.clim %>%
  filter(psi>(-10000)) %>% #remove very low value of psi
  filter(species.binomial==sp) %>% 
  mutate(hsm=hsm/1000,
         fsm=fsmwinter) %>% 
  sample_frac(0.6)  

data.list <- list(N=dim(db.pres)[1],
                  presence=db.pres$presence,
                  fsm=as.numeric(db.pres$fsm),
                  hsm=as.numeric(db.pres$hsm),
                  prior_K=df.traits[df.traits$species.binomial==sp,
                                    "prevalence"],
                  NULL)
# model III
fit.quilIII <- stan(file = "glm_log_1sp_III.stan",
                    data=data.list,
                    iter=1000,
                    chains=2,
                    core=2,
                    include=FALSE,
                    pars=c("proba","K_vect"))
fit.quilIV <- stan(file = "glm_log_1sp_IV.stan",
                   data=data.list,
                   iter=1000,
                   chains=2,
                   core=2,
                   include=FALSE,
                   pars=c("proba","K_vect"))



#%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%
#### Section 15 ####
#'@description fit 4 models for each species and select the best
#'

#_ Fit model ----

species.list=unique(db.clim$species.binomial)
# c("Abies alba",
#   "Betula pendula",
#   "Pinus sylvestris",
#   "Picea abies",
#   "Fagus sylvatica",
#   "Fraxinus excelsior",
#   "Pinus nigra",
#   "Quercus ilex")
df.output=df.traits %>% 
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
set.seed(4554)
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
                    prior_K=df.traits[df.traits$species.binomial==sp,
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
                    prior_K=df.traits[df.traits$species.binomial==sp,
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
                    prior_K=df.traits[df.traits$species.binomial==sp,
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
                    prior_K=df.traits[df.traits$species.binomial==sp,
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

#_ Format output ----

#' Compute divergence
df.output$divergence=NA

for (sp in species.list){
  print(sp)
  
  ### 2 var
  print("2sm")
  load(file=paste0("fit_mod6/",sp,"_fitIII.RData"))
  par_mod=get_sampler_params(fit.2var)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="2sm",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
  # HSM
  print("hsm")
  load(file=paste0("fit_mod6/",sp,"_fitIII_hsm.RData"))
  par_mod=get_sampler_params(fit.hsm)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="hsm",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
  
  # FSM
  print("fsm")
  load(file=paste0("fit_mod6/",sp,"_fitIII_fsm.RData"))
  par_mod=get_sampler_params(fit.fsm)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="fsm",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
  # None
  print("none")
  load(file=paste0("fit_mod6/",sp,"_fitIII_none.RData"))
  par_mod=get_sampler_params(fit.none)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="none",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
}

#' Compute BIC
df.output <- df.output %>% 
  mutate(nb.par=case_when(mod=="2sm"~5,
                          mod=="hsm"|mod=="fsm"~3,
                          mod=="none"~1),
         bic=2*nb.par-2*lp__)
save(df.output,file="fit_mod6/df.output.csv")

#' Model selection
df.mod.select <- df.output %>% 
  filter(Rhat<1.2) %>% 
  filter(divergence <0.1) %>% 
  group_by(species.binomial) %>% 
  slice(which.min(bic)) %>% 
  ungroup()

#' Index of inflexion
hsm.95=5.5#quantile(db.clim$hsm,prob=0.95)[[1]]/1000
hsm.05=-2.5#quantile(db.clim$hsm,prob=0.05)[[1]]/1000
fsm.95=quantile(db.clim$fsmwinter,prob=0.95)[[1]]
fsm.05=quantile(db.clim$fsmwinter,prob=0.05)[[1]]
df.mod.select <- df.mod.select %>% 
  mutate(inflex_fsm=100*((16*exp(-r_fsm*(fsm.95-t_fsm)))/(2+2*exp(-r_fsm*(fsm.95-t_fsm)))^2),#+
           #(K_int*r_fsm*exp(-r_fsm*(fsm.05-t_fsm)))/(1+exp(-r_fsm*(fsm.05-t_fsm)))^2,
         #K_int*r_fsm/4,
         inflex_hsm=100*((16*exp(-r_hsm*(hsm.95-t_hsm)))/(2+2*exp(-r_hsm*(hsm.95-t_hsm)))^2)#+
           #(K_int*r_hsm*exp(-r_hsm*(hsm.05-t_hsm)))/(1+exp(-r_hsm*(hsm.05-t_hsm)))^2
           #K_int*r_hsm/4
         )

#_ Plot parameters ----

#' Slope/slope
df.mod.select %>% 
  ggplot(aes(r_hsm,r_fsm,color=species.binomial,shape=mod,label=species.binomial))+ #size=-log(Rhat),
  geom_point()+
  geom_text(hjust=0,vjust=0)


#' Threshold/threshold
df.mod.select %>% 
  # filter(mod=="2sm") %>%
  ggplot(aes(t_fsm,t_hsm,color=species.binomial,shape=mod,label=species.binomial))+
  geom_text(hjust=0,vjust=0)+
  geom_point(size=5)

#' thresholds/traits contribution
df.mod.select %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm,t_fsm,mod) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to="threshold",values_to="threshold_val") %>% 
  filter((mod%in%c("hsm","2sm")&traits=="PX.mu"&threshold=="t_hsm")|
           (mod%in%c("fsm","2sm")&traits=="LT50.mean"&threshold=="t_fsm")) %>% 
  mutate(app_thresh=-traits_val-threshold_val) %>% 
  ggplot(aes(traits_val,threshold_val,color=species.binomial,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  geom_abline(slope=-1)+
  geom_smooth(method="gam",aes(color=NULL,label=NULL,shape=NULL))+
  facet_wrap(traits~threshold,scales="free")+
  theme(legend.position = "none")+
  theme_minimal()+
  NULL

df.mod.select %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm,t_fsm,mod) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to="threshold",values_to="threshold_val") %>% 
  filter((mod%in%c("hsm","2sm")&traits=="PX.mu"&threshold=="t_hsm")|
           (mod%in%c("fsm","2sm")&traits=="LT50.mean"&threshold=="t_fsm")) %>% 
  mutate(app_thresh=traits_val-threshold_val) %>% 
  ggplot(aes(traits_val,app_thresh,color=species.binomial,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  geom_abline(slope=1)+
  geom_smooth(method="gam",aes(color=NULL,label=NULL,shape=NULL))+
  facet_wrap(~traits,scales="free")+
  theme(legend.position = "none")+
  theme_minimal()+
  NULL


#' thresholds/traits contribution
df.mod.select %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm,t_fsm,mod,inflex_fsm,inflex_hsm) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to="threshold",values_to="threshold_val") %>% 
  pivot_longer(cols=c("inflex_fsm","inflex_hsm"),names_to="inflex",values_to="inflex_val") %>% 
  filter((traits=="PX.mu"&threshold=="t_hsm"&inflex=="inflex_hsm")|
           (traits=="LT50.mean"&threshold=="t_fsm"&inflex=="inflex_fsm")) %>% 
  filter((mod=="hsm"&threshold=="t_hsm")|
           (mod=="fsm"&threshold=="t_fsm")|
           (mod=="2sm"&threshold%in%c("t_hsm","t_fsm"))) %>% 
  mutate(contribution=100*(-traits_val-threshold_val)/(-traits_val)) %>% 
  ggplot(aes(contribution,species.binomial,color=inflex_val,shape=mod))+
  geom_point(size=3)+
  geom_vline(xintercept = 100)+
  scale_color_gradientn(colours = terrain.colors(5))+
  facet_wrap(~traits)

df.mod.select %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm,t_fsm,mod,inflex_fsm,inflex_hsm) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to="threshold",values_to="threshold_val") %>% 
  pivot_longer(cols=c("inflex_fsm","inflex_hsm"),names_to="inflex",values_to="inflex_val") %>% 
  filter((traits=="PX.mu"&threshold=="t_hsm"&inflex=="inflex_hsm")|
           (traits=="LT50.mean"&threshold=="t_fsm"&inflex=="inflex_hsm")) %>% 
  filter((mod=="hsm"&threshold=="t_hsm")|
           (mod=="fsm"&threshold=="t_fsm")|
           (mod=="2sm"&threshold%in%c("t_hsm","t_fsm"))) %>% 
  mutate(contribution=100*(threshold_val)/(-traits_val)) %>% 
  ggplot(aes(traits_val,contribution,color=inflex_val,shape=mod))+
  geom_point(size=3)+
  scale_color_gradientn(colours = viridis(5))+
  geom_smooth(method="gam",aes(color=NULL,shape=NULL))+
  facet_wrap(~traits,scales="free")+
  theme_minimal()


#' slopes/traits
df.mod.select %>% 
  select(species.binomial,PX.mu,LT50.mean,r_hsm,r_fsm,mod) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("r_fsm","r_hsm"),names_to="slope",values_to="slope_val") %>% 
  filter((traits=="PX.mu"&slope=="r_hsm")|
           (traits=="LT50.mean"&slope=="r_fsm")) %>% 
  filter(slope_val<20) %>% 
  ggplot(aes(traits_val,slope_val,color=species.binomial,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  theme(legend.position = "none")+
  facet_wrap(traits~slope,scales="free")

#' Inflexion



#_ Fit curves ----
hsm.95=quantile(db.clim$hsm,prob=0.95)[[1]]/1000
hsm.05=quantile(db.clim$hsm,prob=0.05)[[1]]/1000
fsm.95=quantile(db.clim$fsmwinter,prob=0.95)[[1]]
fsm.05=quantile(db.clim$fsmwinter,prob=0.05)[[1]]

#' 1 plot by species

# for hsm
# different type of FSM
df.mod.select %>% 
  crossing(hsm=seq(hsm.05,
                   hsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("low","indif","high"),
                      fsm=c(fsm.05,0,fsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(hsm,pred,color=fsm.type))+
  geom_line(size=1)+
  # geom_hline(aes(yintercept=prevalence))+
  facet_wrap(~species.binomial,scales="free", nrow = 3)
# colored by inflexion
df.mod.select %>% 
  crossing(hsm=seq(hsm.05,
                   hsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("indif"),
                      fsm=c(0))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(hsm,pred,color=inflex_hsm,linesize=species.binomial))+
  geom_line(size=0.6)+
  scale_color_gradientn(colours = viridis(15))+
  # geom_hline(aes(yintercept=prevalence))+
  # scale_y_log10()+
  theme_minimal()+
  NULL
df.mod.select %>% 
    crossing(fsm=seq(fsm.05,
                     fsm.95,
                     length.out=400)) %>% 
    crossing(data.frame(hsm.type=c("indif"),
                        hsm=c(0))) %>% 
    mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                         (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
    ggplot(aes(fsm,pred,color=inflex_fsm,linesize=species.binomial))+
    geom_line(size=0.6)+
    scale_color_gradientn(colours = viridis(15))+
    # geom_hline(aes(yintercept=prevalence))+
    # scale_y_log10()+
    theme_minimal()+
  NULL

#for fsm
# different type of HSM
df.mod.select %>% 
  crossing(fsm=seq(fsm.05,
                   fsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(hsm.type=c("low","indif","high"),
                      hsm=c(hsm.05,0,hsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(fsm,pred,color=hsm.type))+
  geom_line(size=1)+
  # geom_hline(aes(yintercept=prevalence))+
  facet_wrap(~species.binomial,scales="free", nrow = 3)

#' for a species selection
df.mod.select %>% 
  filter(species.binomial%in%c("Picea abies","Fagus sylvatica","Ilex aquifolium")) %>% 
  crossing(hsm=seq(hsm.05,
                   hsm.95,
                   length.out=400)) %>% 
  crossing(data.frame(fsm.type=c("low","indif","high"),
                      fsm=c(fsm.05,0,fsm.95))) %>% 
  mutate(pred=K_int/((1+exp(-r_fsm*(fsm-t_fsm)))*
                       (1+exp(-r_hsm*(hsm-t_hsm))))) %>% 
  ggplot(aes(hsm,pred))+
  geom_line(aes(color=fsm.type),size=1)+ #as.factor(PX.mu)
  facet_wrap(~species.binomial,scales="free")



#_ Parameters/predictors ----
# graph pour voir valeur des paramètre en fonction des modèles
# ID: voir a quel point le choix du modèlee donne une tendance ou non dans les paramètres
df.mod.select %>% 
  pivot_longer(cols=c("K_int","r_fsm","r_hsm","t_hsm","t_fsm"),
               names_to = "parameter",
               values_to = "val") %>% 
  filter(val!=0) %>% 
  ggplot(aes(mod,val))+
  geom_boxplot()+
  facet_wrap(~parameter,scales="free")


# compare mod selection to niche caracteristics of species
# ID : species with extreme niches are explained by only one sm.
df.niche <- db.clim %>% 
  filter(presence==1) %>% 
  filter(psi>(-20000)) %>% 
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
df.mod.select %>% 
  left_join(df.niche,by="species.binomial") %>% 
  # filter(!species.binomial%in%c("Picea abies","Pinus cembra","Quercus ilex","Acer monspessulanum")) %>% 
  pivot_longer(cols=c("K_int","r_fsm","r_hsm","t_hsm","t_fsm"),
               names_to = "parameter",
               values_to = "val") %>% 
  filter(!(mod=="hsm"&parameter%in%c("r_fsm","t_fsm"))) %>% 
  filter(!(mod=="fsm"&parameter%in%c("r_hsm","t_hsm"))) %>% 
  pivot_longer(cols=c("prevalence","lat.mean","lat.sd","long.mean","long.sd","lat.q05","lat.q95"),
               names_to="predictor",
               values_to = "val_pred") %>% 
  ggplot(aes(val_pred,val,color=mod))+
  geom_point()+
  geom_smooth(method="lm",aes(color=NULL))+
  facet_wrap(parameter~predictor,scales="free",nrow = 5)
df.mod.select %>% 
  left_join(df.niche,by="species.binomial") %>% 
  pivot_longer(cols=c("prevalence","lat.mean","lat.sd","long.mean","long.sd","lat.q05","lat.q95"),
               names_to="predictor",
               values_to = "val_pred") %>% 
  ggplot(aes(mod,val_pred)) +
  geom_boxplot() +
  facet_wrap(~predictor,scales="free")

#' Continentality / lat
df.mod.select %>% 
  left_join(df.niche,by="species.binomial") %>% 
  left_join(df.cont,by="species.binomial") %>% 
  pivot_longer(cols=c("jci","lat.mean","long.mean","long.sd","lat.q05","lat.q95"),
               names_to="predictor",
               values_to = "val_pred") %>% 
  ggplot(aes(mod,val_pred)) +
  geom_boxplot() +
  facet_wrap(~predictor,scales="free")
#' Shade tolerance
df.shadetol=read.csv2("data/Species traits/data_Niinemets&Valladares_2006.csv")


df.mod.select %>% 
  left_join(df.shadetol,by=c("species.binomial"="Species")) %>% 
  mutate(across(c(shade_tolerance.mean,drought_tolerance.mean,waterlogging_tolerance.mean),
                as.numeric)) %>% 
  pivot_longer(cols=c("K_int","r_fsm","r_hsm","t_hsm","t_fsm","auc"),
               names_to = "parameter",
               values_to = "val") %>% 
  pivot_longer(cols=c("shade_tolerance.mean","drought_tolerance.mean","waterlogging_tolerance.mean"),
               names_to="predictor",
               values_to = "val_pred") %>% 
  filter(!(mod=="hsm"&parameter%in%c("r_fsm","t_fsm"))) %>% 
  filter(!(mod=="fsm"&parameter%in%c("r_hsm","t_hsm"))) %>% 
  ggplot(aes(val_pred,val,color=Group))+
  geom_point()+
  geom_smooth(method="lm",aes(color=NULL))+
  facet_wrap(predictor~parameter,scales="free",nrow=3)


# distributional overlap
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

df.mod.select %>% 
  left_join(df.overlap,by=c("species.binomial")) %>% 
  pivot_longer(cols=c("K_int","r_fsm","r_hsm","t_hsm","t_fsm"),
               names_to = "parameter",
               values_to = "val") %>% 
  pivot_longer(cols=c("overlap","overlap.hsm","overlap.fsm"),
               names_to="overlap",
               values_to = "overlap_val") %>% 
  filter(!(mod=="hsm"&parameter%in%c("r_fsm","t_fsm"))) %>% 
  filter(!(mod=="fsm"&parameter%in%c("r_hsm","t_hsm"))) %>% 
  ggplot(aes(overlap_val,val,color=mod))+
  geom_point()+
  geom_smooth(method="gam",aes(color=NULL))+
  facet_wrap(overlap~parameter,scales="free",nrow=3)


#' inflexion
# mod type
df.mod.select %>% 
  pivot_longer(cols=c("inflex_fsm","inflex_hsm"),
               names_to = "inflex",
               values_to="inflex_val") %>% 
  filter(!(mod=="hsm"&inflex%in%c("inflex_fsm"))) %>% 
  filter(!(mod=="fsm"&inflex%in%c("inflex_hsm"))) %>% 
  ggplot(aes(mod,inflex_val,color=mod))+
  geom_boxplot()+
  scale_y_log10()+
  facet_wrap(~inflex,scales="free")
# niche caract
df.mod.select %>% 
  filter(mod %in% c('fsm','2sm')) %>% 
  left_join(df.niche,by="species.binomial") %>% 
  pivot_longer(cols=c("prevalence","lat.mean","lat.sd","long.mean","long.sd","lat.q05","lat.q95"),
               names_to="predictor",
               values_to = "val_pred") %>% 
  ggplot(aes(log(inflex_fsm),val_pred)) +
  geom_point() +
  geom_smooth(method="lm")+
  facet_wrap(~predictor,scales="free")
df.mod.select %>% 
  filter(mod %in% c('fsm','2sm')) %>% 
  left_join(df.niche,by="species.binomial") %>% 
  ggplot(aes(lat.q95,log(inflex_fsm),label=species.binomial))+
  geom_text()+
  geom_point()+
  geom_smooth(method="gam")
# overlap
df.mod.select %>% 
  left_join(df.overlap,by=c("species.binomial")) %>% 
  pivot_longer(cols=c("overlap","overlap.hsm","overlap.fsm"),
               names_to="overlap",
               values_to = "overlap_val") %>% 
  pivot_longer(cols=c("inflex_fsm","inflex_hsm"),
               names_to = "inflex",
               values_to="inflex_val") %>% 
  filter(!(mod=="hsm"&inflex%in%c("inflex_fsm"))) %>% 
  filter(!(mod=="fsm"&inflex%in%c("inflex_hsm"))) %>% 
  ggplot(aes(overlap_val,log(inflex_val),color=mod))+
  geom_point()+
  geom_smooth(method="gam",aes(color=NULL))+
  facet_wrap(overlap~inflex,scales="free",nrow=3)
#shadetol
df.mod.select %>% 
  left_join(df.shadetol,by=c("species.binomial"="Species")) %>% 
  mutate(across(c(shade_tolerance.mean,drought_tolerance.mean,waterlogging_tolerance.mean),
                as.numeric)) %>% 
  pivot_longer(cols=c("shade_tolerance.mean","drought_tolerance.mean","waterlogging_tolerance.mean"),
               names_to="predictor",
               values_to = "val_pred") %>% 
  pivot_longer(cols=c("inflex_fsm","inflex_hsm"),
               names_to = "inflex",
               values_to="inflex_val") %>% 
  filter(!(mod=="hsm"&inflex%in%c("inflex_fsm"))) %>% 
  filter(!(mod=="fsm"&inflex%in%c("inflex_hsm"))) %>% 
  ggplot(aes(val_pred,log(inflex_val),color=Group))+
  geom_point()+
  geom_smooth(method="lm",aes(color=NULL))+
  facet_wrap(predictor~inflex,scales="free",nrow=3)

#_ analysis of inflexion index ----
# inflexion hsm/fsm
df.mod.select %>% 
  rename(HSM=inflex_hsm,
         FSM=inflex_fsm) %>% 
  pivot_longer(col=c("HSM","FSM")) %>% 
  filter(mod!="none") %>% 
  filter(!(mod=="fsm"&name=="HSM")&
           !(mod=="hsm"&name=="FSM")) %>% 
  ggplot(aes(name,value))+
  geom_boxplot()+
  theme_minimal()+
  theme(axis.title.x = element_blank())+
  ylab("Inflexion index")

# inflexion / trait
df.mod.select %>% 
  filter(mod!="none") %>% 
  select(species.binomial,PX.mu,LT50.mean,inflex_hsm,inflex_fsm,mod) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("inflex_hsm","inflex_fsm"),names_to="inflex",values_to="inflex_val") %>% 
  filter(!(mod=="fsm"&((inflex=="inflex_hsm")|(traits=="PX.mu")))&
           !(mod=="hsm"&((inflex=="inflex_fsm")|(traits=="LT50.mean")))) %>% 
  ggplot(aes(traits_val,inflex_val,color=mod,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  geom_smooth(method="lm",aes(color=NULL,shape=NULL,label=NULL))+
  # theme(legend.position = "none")+
  facet_wrap(~traits,scales="free")+
  theme_minimal()

# inflexion / preval
df.mod.select %>% 
  filter(mod!="none") %>% 
  select(species.binomial,prevalence,inflex_hsm,inflex_fsm,mod) %>% 
  pivot_longer(cols=c("inflex_hsm","inflex_fsm"),names_to="inflex",values_to="inflex_val") %>% 
  filter(!(mod=="fsm"&((inflex=="inflex_hsm")))&
           !(mod=="hsm"&((inflex=="inflex_fsm")))) %>% 
  ggplot(aes(prevalence,inflex_val,color=mod,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  geom_smooth(method="lm",aes(color=NULL,shape=NULL,label=NULL))+
  # theme(legend.position = "none")+
  facet_wrap(~inflex,scales="free")+
  theme_minimal()

# inflexion / seuil
df.mod.select %>% 
  filter(mod!="none") %>% 
  select(species.binomial,t_fsm,t_hsm,inflex_hsm,inflex_fsm,mod) %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to = "threshold",values_to = "threshold_val") %>% 
  filter(!(mod=="hsm"&threshold=="t_fsm")&
           !(mod=="fsm"&threshold=="t_hsm")) %>% 
  pivot_longer(cols=c("inflex_hsm","inflex_fsm"),names_to="inflex",values_to="inflex_val") %>% 
  filter(!(mod=="hsm"&inflex=="inflex_fsm")&
           !(mod=="fsm"&inflex=="inflex_hsm")) %>% 
  ggplot(aes(inflex_val,threshold_val,color=mod,label=species.binomial,shape=mod))+
  geom_point(size=3)+
  geom_text(hjust=0,vjust=0)+
  geom_smooth(method="lm",aes(color=NULL,shape=NULL,label=NULL))+
  # theme(legend.position = "none")+
  facet_wrap(~inflex,scales="free")+
  theme_minimal()

# contribution / inflexion 
df.mod.select %>% 
  select(species.binomial,PX.mu,LT50.mean,t_hsm,t_fsm,mod,inflex_fsm,inflex_hsm) %>% 
  pivot_longer(cols=c("PX.mu","LT50.mean"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("t_fsm","t_hsm"),names_to="threshold",values_to="threshold_val") %>% 
  pivot_longer(cols=c("inflex_fsm","inflex_hsm"),names_to="inflex",values_to="inflex_val") %>% 
  filter((traits=="PX.mu"&threshold=="t_hsm"&inflex=="inflex_hsm")|
           (traits=="LT50.mean"&threshold=="t_fsm"&inflex=="inflex_fsm")) %>% 
  filter((mod=="hsm"&threshold=="t_hsm")|
           (mod=="fsm"&threshold=="t_fsm")|
           (mod=="2sm"&threshold%in%c("t_hsm","t_fsm"))) %>% 
  mutate(contribution=100*(-traits_val-threshold_val)/(-traits_val)) %>% 
  ggplot(aes(contribution,label=species.binomial,y=inflex_val,shape=mod))+
  geom_point(size=3)+
  scale_color_gradientn(colours = terrain.colors(5))+
  geom_smooth(method="gam",aes(shape=NULL,label=NULL))+
  facet_wrap(~inflex)+
  theme_minimal()

#_ Statistical relation between vars ----

for (par in c("K_int","r_fsm","r_hsm","t_hsm","t_fsm")){
  print(par)
  df.lm <- df.mod.select %>% 
    # filter(!species.binomial%in%c("Picea abies","Pinus cembra","Quercus ilex","Acer monspessulanum")) %>% 
    left_join(df.shadetol,by=c("species.binomial"="Species")) %>% 
    mutate(across(c(shade_tolerance.mean,drought_tolerance.mean,waterlogging_tolerance.mean),
                  as.numeric)) %>% 
    filter(!is.na(shade_tolerance.mean))
  par.list=df.lm[,par][[1]]
  tol.lm=lm(par.list~shade_tolerance.mean+drought_tolerance.mean+waterlogging_tolerance.mean+mod,
            data=df.lm)
  print(summary(tol.lm))
  df.lm <- df.mod.select %>% 
    # filter(!species.binomial%in%c("Picea abies","Pinus cembra","Quercus ilex","Acer monspessulanum")) %>% 
    left_join(df.niche,by="species.binomial") 
  par.list=df.lm[,par][[1]]
  niche.lm=lm(par.list~prevalence+lat.mean+lat.sd+long.mean+long.sd+mod,
              data=df.lm)
  print(summary(niche.lm))
}


#_ Compute prevalence ----
df.mod.select$preval.theo=NA
for (sp in df.mod.select$species.binomial){
  print(sp)
  if (df.mod.select[df.mod.select$species.binomial==sp,"mod"]=="2sm"){
    print("2sm")
    preval=db.clim %>% 
      filter(species.binomial==sp) %>% 
      mutate(pres.theo=case_when(hsm>0&fsmwinter>0~1,
                                 TRUE~0)) %>% 
      group_by(pres.theo) %>% 
      summarise(prop=sum(presence==1)/n())
    df.mod.select[df.mod.select$species.binomial==sp,"preval.theo"]=
      preval[preval$pres.theo==1,"prop"]/(preval[preval$pres.theo==1,"prop"]+preval[preval$pres.theo==0,"prop"])
    
  }else if (df.mod.select[df.mod.select$species.binomial==sp,"mod"]=="hsm"){
    print("hsm")
    preval=db.clim %>% 
      filter(species.binomial==sp) %>% 
      mutate(pres.theo=case_when(hsm>0~1,
                                 TRUE~0)) %>% 
      group_by(pres.theo) %>% 
      summarise(prop=sum(presence==1)/n())
    df.mod.select[df.mod.select$species.binomial==sp,"preval.theo"]=
      preval[preval$pres.theo==1,"prop"]/(preval[preval$pres.theo==1,"prop"]+preval[preval$pres.theo==0,"prop"])
    
  }else if (df.mod.select[df.mod.select$species.binomial==sp,"mod"]=="fsm"){
    print("fsm")
    preval=db.clim %>% 
      filter(species.binomial==sp) %>% 
      mutate(pres.theo=case_when(fsmwinter>0~1,
                                 TRUE~0)) %>% 
      group_by(pres.theo) %>% 
      summarise(prop=sum(presence==1)/n())
    df.mod.select[df.mod.select$species.binomial==sp,"preval.theo"]=
      preval[preval$pres.theo==1,"prop"]/(preval[preval$pres.theo==1,"prop"]+preval[preval$pres.theo==0,"prop"])
    
  }
}

df.mod.select %>% 
  pivot_longer(cols=c("K_int","r_fsm","r_hsm","t_hsm","t_fsm"),
               names_to = "parameter",
               values_to = "val") %>% 
  filter(!(mod=="hsm"&parameter%in%c("r_fsm","t_fsm"))) %>% 
  filter(!(mod=="fsm"&parameter%in%c("r_hsm","t_hsm"))) %>% 
  ggplot(aes(preval.theo,val,color=mod))+
  geom_point()+
  geom_smooth(method="gam",aes(color=NULL))+
  facet_wrap(~parameter,scales="free",nrow=3)


#_ Compute AUC ----
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



#_ continentality index ----
library(raster)
rast.cont=as.data.frame(mean(rast("jci_year.nc")),xy=TRUE) %>% 
  mutate(z=x,
         x=y,
         y=z) %>% 
  select(-z) %>% 
  rast(crs="epsg:4326")
df.cont=data.frame(species.binomial=df.mod.select$species.binomial,
                   jci=NA)
for (sp in df.mod.select$species.binomial){
  print(sp)
  db.pres <- db.clim %>% 
    filter(species.binomial==sp) %>% 
    filter(presence==1)
  db.pres <- cbind(db.pres,
                   jci=extract(rast.cont,db.pres[,c("x","y")])[["mean"]])
  jci=mean(db.pres$jci,na.rm=TRUE)
  df.cont[df.cont$species.binomial==sp,"jci"]=jci
}
rm(rast.cont)


#### Section 16 ####
#'@description fit linear models
#'


#_ Fit model ----

species.list=unique(db.clim$species.binomial)
# c("Abies alba",
#   "Betula pendula",
#   "Pinus sylvestris",
#   "Picea abies",
#   "Fagus sylvatica",
#   "Fraxinus excelsior",
#   "Pinus nigra",
#   "Quercus ilex")
df.output=df.traits %>% 
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
set.seed(4554)
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
                    hsm=as.numeric(db.pres$hsm)
                    )
  fit.2var <- stan(file = "glm_lin_1sp_2var.stan",
                   data=data.list,
                   iter=1000,
                   chains=3,
                   core=3,
                   include=FALSE,
                   pars=c("proba"))
  save(fit.2var,file=paste0("fit_mod7/",sp,"_fit2var.RData"))
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
                    hsm=as.numeric(db.pres$hsm)
                    )
  fit.hsm <- stan(file = "glm_lin_1sp_1varh.stan",
                  data=data.list,
                  iter=1000,
                  chains=3,
                  core=3,
                  include=FALSE,
                  pars=c("proba"))
  save(fit.hsm,file=paste0("fit_mod7/",sp,"_fit_hsm.RData"))
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
                    fsm=as.numeric(db.pres$fsm)
                    )
  fit.fsm <- stan(file = "glm_lin_1sp_1var.stan",
                  data=data.list,
                  iter=1000,
                  chains=3,
                  core=3,
                  include=FALSE,
                  pars=c("proba"))
  save(fit.fsm,file=paste0("fit_mod7/",sp,"_fit_fsm.RData"))
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
  
  
}

#_ Format output ----

#' Compute divergence

for (sp in species.list){
  print(sp)
  
  ### 2 var
  print("2sm")
  load(file=paste0("fit_mod6/",sp,"_fitIII.RData"))
  par_mod=get_sampler_params(fit.2var)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="2sm",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
  # HSM
  print("hsm")
  load(file=paste0("fit_mod6/",sp,"_fitIII_hsm.RData"))
  par_mod=get_sampler_params(fit.hsm)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="hsm",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
  
  # FSM
  print("fsm")
  load(file=paste0("fit_mod6/",sp,"_fitIII_fsm.RData"))
  par_mod=get_sampler_params(fit.fsm)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="fsm",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
  # None
  print("none")
  load(file=paste0("fit_mod6/",sp,"_fitIII_none.RData"))
  par_mod=get_sampler_params(fit.none)
  df.output[df.output$species.binomial==sp
            & df.output$mod=="none",
            "divergence"]=
    mean(sapply(par_mod, function(x) mean(x[, "divergent__"])))
  
}

#' Compute BIC
df.output <- df.output %>% 
  mutate(nb.par=case_when(mod=="2sm"~5,
                          mod=="hsm"|mod=="fsm"~3,
                          mod=="none"~1),
         bic=2*nb.par-2*lp__)
write.csv(df.output,file="fit_mod6/df.output.csv")

#' Model selection
df.mod.select <- df.output %>% 
  filter(Rhat<1.2) %>% 
  filter(divergence <0.1) %>% 
  group_by(species.binomial) %>% 
  slice(which.min(bic)) %>% 
  ungroup()

#### Annexe ####
#'@description check figure
#'
#Proba per class
breaks <- seq(min(db.pres$fsmwinter),max(db.pres$fsmwinter),length.out=30)
labels <- round((breaks[2:30]+breaks[1:29])/2,digit=1)
# db.prob <- db.pres %>% 
db.prob <- db.pres %>% 
  sample_frac(0.01) %>% 
  mutate(fsmwinter=cut(fsmwinter,breaks=breaks,labels=labels),
         fsmwinter=as.numeric(as.character(fsmwinter))) %>% 
  group_by(species.binomial,fsmwinter) %>% 
  summarise(prob=sum(presence==1)/n()) %>% 
  filter(!is.na(fsmwinter)) %>% 
  ungroup()
db.prob %>% 
  filter(species.binomial %in% c("Quercus ilex")) %>% 
  ggplot(aes(fsmwinter,prob,color=species.binomial))+
  geom_point()


# proba per class but for 2 dimensions
#Proba per class
breaks.fsm <- seq(min(db.pres$fsmwinter),max(db.pres$fsmwinter),length.out=30)
labels.fsm <- round((breaks.fsm[2:30]+breaks.fsm[1:29])/2,digit=1)
breaks.hsm <- seq(min(db.pres$hsm),max(db.pres$hsm),length.out=30)
labels.hsm <- round((breaks.hsm[2:30]+breaks.hsm[1:29])/2,digit=1)
db.prob.2d<- db.pres %>% 
  mutate(fsmwinter=cut(fsmwinter,breaks=breaks.fsm,labels=labels.fsm),
         fsmwinter=as.numeric(as.character(fsmwinter)),
         hsm=cut(hsm,breaks=breaks.hsm,labels=labels.hsm),
         hsm=as.numeric(as.character(hsm))) %>% 
  group_by(species.binomial,fsmwinter,hsm) %>% 
  summarise(prob=sum(presence==1)/n()) %>% 
  filter(!is.na(fsmwinter)) %>% 
  filter(!is.na(hsm)) %>% 
  ungroup()
db.prob.2d %>% 
  # filter(species.binomial %in% c("Abies alba")) %>% 
  ggplot(aes(fsmwinter,hsm,color=prob))+
  geom_point(size=5)+
  facet_wrap(~species.binomial)


x <- seq(-1,1, len = 100)
z <- seq(-1,1, len = 100)
df<- expand.grid(x = x, z = z)
slope1 <- 1

slope2 <- 2

interc <- 0.5
df$y <- (df$x>-0.3)*(interc) + (df$x<=-0.3)*((interc +0.3*slope1)+slope1*df$x) + (df$z>-0.3)*(interc) + (df$z<=-0.3)*((interc+0.3*slope2)+slope2*df$z)

df$y.ilogit <- boot::inv.logit(df$y)

image(x=x,y= z, matrix(df$y, 100, 100))
ggplot(df, aes(x, y, color = z)) + geom_point()
ggplot(df, aes(x, y.ilogit, color = z)) + geom_point()

sm=seq(-5,5,by=0.1)
data.frame(sm=sm) %>% 
  crossing(r=c(0.01,0.1,1,10)) %>% 
  crossing(t=c(-2,-1,1,2)) %>% 
  crossing(K=c(0.01,0.1,0.9)) %>% 
  mutate(prob=K/(1+exp(-r*(sm-t)))) %>% 
  ggplot(aes(sm,prob,color=as.factor(r)))+
  geom_line()+                        
  facet_grid(K~t+.,scales="free_y")

sm=seq(-5,5,by=0.1)
data.frame(sm=sm) %>% 
  crossing(a=c(0.1,1,10)) %>% 
  crossing(b=c(-1,0.1,1,10)) %>% 
  crossing(c=c(-2,0.2,2)) %>% 
  mutate(prob=exp(a+b*sm)/((exp(a+b*sm)+1))) %>% 
  ggplot(aes(sm,prob,color=as.factor(b)))+
  geom_line()+                        
  facet_grid(a~c+.,scales="free_y")

# compare quantile of psi/temp with traits
db.clim %>% 
  filter(presence==1) %>% 
  filter(psi>(-15000)) %>% 
  mutate(psi=-log(max(psi+1)-psi)) %>% 
  group_by(species.binomial) %>% 
  summarise(psi05=quantile(psi,prob=0.05),
            psi95=quantile(psi,prob=0.95),
            t05=quantile(frost.winter,prob=0.05),
            t95=quantile(frost.winter,prob=0.95)) %>% 
  left_join(df.traits[,c("species.binomial","PX.mu","LT50.mean","Group")],by="species.binomial") %>% 
  # pivot_longer(cols=c("psi05","psi95","t05","t95"),
  #              names_to = "quantile",
  #              values_to = "quant_val") %>%
  pivot_longer(cols=c("psi05","t05"),
               names_to = "start",
               values_to = "start_val") %>%
  pivot_longer(cols=c("psi95","t95"),
               names_to = "end",
               values_to = "end_val") %>% 
  pivot_longer(cols = c("PX.mu","LT50.mean"),
               names_to="trait",
               values_to="trait_val") %>% 
  filter(!(trait=="PX.mu"&(start=="t05"|end=="t95"))) %>% 
  filter(!(trait=="LT50.mean"&(start=="psi05"|end=="psi95"))) %>% 
  # filter(!(trait=="PX.mu"&quantile%in%c("t05","t95"))) %>%
  # filter(!(trait=="LT50.mean"&quantile%in%c("psi05","psi95"))) %>%
  ggplot()+
  geom_point(aes(trait_val,start_val,color=Group))+
  geom_point(aes(trait_val,end_val,color=Group))+
  geom_segment(aes(x=trait_val,xend=trait_val,y=start_val,yend=end_val),alpha=0.3)+
  facet_wrap(~trait,scales="free")
  


# compare quantile of psi/temp with traits
db.clim %>% 
  filter(presence==1) %>% 
  filter(psi>(-15000)) %>% 
  # mutate(psi=-log(max(psi+1)-psi)) %>% 
  group_by(species.binomial) %>% 
  summarise(psi05=quantile(psi,prob=0.05),
            psi95=quantile(psi,prob=0.95),
            t05=quantile(frost.winter,prob=0.05),
            t95=quantile(frost.winter,prob=0.95)) %>% 
  left_join(df.traits[,c("species.binomial","PX.mu","LT50.mean","Group")],by="species.binomial") %>% 
  # pivot_longer(cols=c("psi05","psi95","t05","t95"),
  #              names_to = "quantile",
  #              values_to = "quant_val") %>%
  pivot_longer(cols=c("psi05","t05"),
               names_to = "start",
               values_to = "start_val") %>%
  pivot_longer(cols=c("psi95","t95"),
               names_to = "end",
               values_to = "end_val") %>% 
  pivot_longer(cols = c("PX.mu","LT50.mean"),
               names_to="trait",
               values_to="trait_val") %>% 
  filter(!(trait=="PX.mu"&(start=="t05"|end=="t95"))) %>% 
  filter(!(trait=="LT50.mean"&(start=="psi05"|end=="psi95"))) %>% 
  # filter(!(trait=="PX.mu"&quantile%in%c("t05","t95"))) %>%
  # filter(!(trait=="LT50.mean"&quantile%in%c("psi05","psi95"))) %>%
  ggplot()+
  geom_point(aes(trait_val,start_val-end_val,color=Group))+
  facet_wrap(~trait,scales="free")


trait_P50=read.csv(file="output/df_P50_filtered.csv")
trait_LT50=read.csv(file="output/df_LT50_filtered.csv")

db.cont[is.na(db.cont$presence)] <- 0
to_sample= db.cont %>% 
  group_by(species.binomial) %>% 
  summarise(prop=sum(presence==1)/n()) %>% 
  left_join(trait_P50[,c("species.binomial","P50.mu","Group","P88.mu")],
            by="species.binomial") %>% 
  left_join(trait_LT50[,c("Species","LT50.mean","data.quality")],
            by=c("species.binomial"="Species")) %>% 
  filter(is.na(LT50.mean)|
           data.quality==7|
           (is.na(P88.mu)&Group=="angiosperm")|
           is.na(P50.mu)) %>% 
  filter(prop>0.00001)
