## load required packages
lapply(c("data.table","tidyr","dplyr","rstan","boot"),require,character.only=TRUE)

## load data
db.clim=fread("output/db_EuForest.csv")

df.traits=read.csv("output/df_trait_filtered.csv") %>% 
  filter(species.binomial %in% unique(db.clim$species.binomial)) %>% 
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

db.pres <- db.clim %>%
  rename_with(.cols=everything(),
              tolower) %>%
  filter(species.binomial %in% df.traits$species.binomial) %>% 
  # filter(species.binomial %in% c("Abies alba","Quercus ilex","Olea europaea", 
  #                                "Pinus pinaster","Quercus robur","Fagus sylvatica",
  #                                "Betula pendula","Picea abies")) %>% 
  filter(!is.na(hsm)) %>% 
  filter(!is.na(fsmwinter)) %>% 
  filter(!is.na(fsmspring)) %>% 
  filter(!is.na(sgdd)) %>% 
  filter(!is.na(wai)) %>% 
  sample_frac(0.05)
rm(db.clim)
hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)

db.pres <- db.pres %>% 
  mutate(hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
         fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1],
         pet=scale(pet,scale=TRUE,center=TRUE)[, 1],
         sgdd=scale(sgdd,scale=TRUE,center=TRUE)[, 1]) 


#species equivalence with factor number
prop.sp <- db.pres %>% 
  group_by(species.binomial) %>% 
  mutate(ntot=n()) %>% 
  summarize(nprop=100*sum(presence==1)/ntot) %>% 
  unique()
 
prop.sp.occ <- db.pres %>% 
  filter(fsmwinter>-fsm.mu/fsm.sd&hsm>-hsm.mu/hsm.sd) %>% 
  group_by(species.binomial) %>% 
  mutate(ntot=n()) %>% 
  summarize(npropocc=100*sum(presence==1)/ntot) %>% 
  unique()

species <- data.frame(species=unique(factor(db.pres$species.binomial)),
                      species.nb=as.character(as.numeric(unique(factor(db.pres$species.binomial))))) %>% 
  left_join(prop.sp,by=c("species"="species.binomial")) %>% 
  left_join(prop.sp.occ,by=c("species"="species.binomial"))
rm(prop.sp,prop.sp.occ)

# stan model, S species, FSM & HSM
data.list <- list(N=dim(db.pres)[1],
                  S=nlevels(factor(db.pres$species.binomial)),
                  presence=db.pres$presence,
                  species=as.numeric(factor(db.pres$species.binomial)),
                  fsm=db.pres$fsmwinter,
                  hsm=db.pres$hsm,
                  prior_threshold_fsm=-fsm.mu/fsm.sd,
                  prior_threshold_hsm=-hsm.mu/hsm.sd)
fit.2sm <- stan(file = "glm_seg_all.stan",
                data=data.list,
                iter=2000,
                chains=4,
                cores=4)
save(fit.2sm,file="fit_allsp_all.RData")
 
# fit.2sm.slope <- stan(file = "glm_seg_allslope.stan",
#                 data=data.list,
#                 iter=2000,
#                 chains=4,
#                 cores=4)
# save(fit.2sm,file="fit_allsp_allslope.RData")

### plots ### don't run on server
# load("df_fit_allsp.RData")
load("fit_allsp_all.RData")
df.fit=as.data.frame(fit.2sm)
fit.sample= df.fit%>%
  sample_frac(0.2) %>% 
  # slice(1001:2000) %>%
  mutate(it=row_number()) %>%
  crossing(fsm=seq(min(db.pres$fsmwinter),max(db.pres$fsmwinter),by=0.1)) %>%
  mutate(hsm=median(db.pres$fsmwinter)) %>%
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>%
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>%
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>%
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>%
  filter(species1==species2) %>%
  select(-species2) 
fit.sample= fit.sample%>%
  mutate(fsm_app=case_when(fsm<threshold_fsm~0,TRUE~1),
         hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
         pred= inv.logit((1-fsm_app)*(b_fsm * (fsm - threshold_fsm) + plateau_fsm)
                         + fsm_app * plateau_fsm
                         + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm)
                         + hsm_app * plateau_hsm))


fit.sample %>%
  left_join(species,by=c("species1"="species.nb")) %>%
  filter(species1 %in% c("11","12","13","14","15","16","17","18","19")) %>% 
  group_by(species,fsm) %>%
  summarise(mean=mean(pred),
            q95=quantile(pred,probs=0.95),
            q05=quantile(pred,probs=0.05)) %>%
  ggplot()+
  geom_line(aes(fsm,mean))+
  geom_ribbon(aes(fsm,ymin=q05,ymax=q95),alpha=0.2)+
  # geom_point(data=db.pres,aes(fsmwinter,presence))+
  facet_wrap(~species)

# Study of parameters
summary.fit=as.data.frame(t(summary(fit.2sm)$summary)) %>% 
  slice(1) %>% 
  pivot_longer(cols=matches("plateau_fsm"),names_to = "species1",values_to="plateau_fsm") %>%
  mutate(species1= sub(".*\\[([^][]+)].*","\\1",species1)) %>%
  pivot_longer(cols=matches("plateau_hsm"),names_to = "species2",values_to = "plateau_hsm") %>%
  mutate(species2= sub(".*\\[([^][]+)].*","\\1",species2)) %>%
  filter(species1==species2) %>%
  select(-species2) %>% 
  left_join(species,by=c("species1"="species.nb"))

summary.fit %>% 
  ggplot(aes(inv.logit(plateau_fsm),npropocc/100))+
  geom_point()+
  xlim(0,0.5)+
  ylim(0,0.5)

# loop for each spcies

for (sp in species$species){
  print(sp)
  
  # FSM 
  print("fsm")
  #Proba per class
  db.sp <- db.pres[db.pres$species.binomial==sp,]
  breaks <- seq(min(db.sp$fsmwinter),max(db.sp$fsmwinter),length.out=26)
  labels <- round((breaks[2:26]+breaks[1:25])/2,digit=1)
  db.sp.prob <- db.sp %>% 
    mutate(fsmwinter=cut(fsmwinter,breaks=breaks,labels=labels),
           fsmwinter=as.numeric(as.character(fsmwinter))) %>% 
    group_by(species.binomial,fsmwinter) %>% 
    summarise(prob=sum(presence==1)/n()) %>% 
    filter(!is.na(fsmwinter)) %>% 
    ungroup()
  db.sp.prob %>% 
    ggplot(aes(fsmwinter,prob))+
    geom_point()
  
  # posterior 
  fit.sp <- df.fit %>%
    select(plateau_int_fsm,
           plateau_int_hsm,
           matches(paste0("\\b",species[species$species==sp,"species.nb"],"\\b")),
           b_fsm,
           b_hsm,
           threshold_fsm,
           threshold_hsm) %>% 
  rename_at(vars(matches("plateau_fsm")),
            ~"plateau_fsm") %>% 
  rename_at(vars(matches("plateau_hsm")),
              ~"plateau_hsm") %>% 
  crossing(fsm=seq(min(db.sp$fsmwinter),max(db.sp$fsmwinter),by=0.1)) %>%
  mutate(hsm=median(db.sp$hsm)) %>% 
  mutate(fsm_app=case_when(fsm<threshold_fsm~0,TRUE~1),
          hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
          pred= inv.logit((1-fsm_app)*(b_fsm * (fsm - threshold_fsm) + plateau_fsm)
                           + fsm_app * plateau_fsm
                           + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm)
                           + hsm_app * plateau_hsm))
  gg <- fit.sp %>% 
    group_by(fsm) %>%
    summarise(mean=mean(pred),
              q95=quantile(pred,probs=0.95),
              q05=quantile(pred,probs=0.05)) %>%
    ggplot()+
    geom_line(aes(fsm,mean))+
    geom_ribbon(aes(fsm,ymin=q05,ymax=q95),alpha=0.2)+
    geom_point(data=db.sp.prob,aes(fsmwinter,prob))+
    theme_bw()+
    xlab("fsm scaled")+
    ylab("Probability of presence")+
    ggtitle(paste0("Posterior and observed probability for ",sp,
                   "\nHSM fixed"))
  ggsave(paste0("figs3/",sp,"_fsmposterior.png"), plot = gg)
  
  # HSM 
  print("hsm")
  #Proba per class
  db.sp <- db.pres[db.pres$species.binomial==sp,]
  breaks <- seq(min(db.sp$hsm),max(db.sp$hsm),length.out=26)
  labels <- round((breaks[2:26]+breaks[1:25])/2,digit=1)
  db.sp.prob <- db.sp %>% 
    mutate(hsm=cut(hsm,breaks=breaks,labels=labels),
           hsm=as.numeric(as.character(hsm))) %>% 
    group_by(species.binomial,hsm) %>% 
    summarise(prob=sum(presence==1)/n()) %>% 
    filter(!is.na(hsm)) %>% 
    ungroup()
  db.sp.prob %>% 
    ggplot(aes(hsm,prob))+
    geom_point()
  
  # posterior 
  fit.sp <- df.fit %>%
    select(plateau_int_fsm,
           plateau_int_hsm,
           matches(paste0("\\b",species[species$species==sp,"species.nb"],"\\b")),
           b_fsm,
           b_hsm,
           threshold_fsm,
           threshold_hsm) %>% 
    rename_at(vars(matches("plateau_fsm")),
              ~"plateau_fsm") %>% 
    rename_at(vars(matches("plateau_hsm")),
              ~"plateau_hsm") %>% 
    crossing(hsm=seq(min(db.sp$hsm),max(db.sp$hsm),by=0.1)) %>%
    mutate(fsm=median(db.sp$fsmwinter)) %>% 
    mutate(fsm_app=case_when(fsm<threshold_fsm~0,TRUE~1),
           hsm_app=case_when(hsm<threshold_hsm~0,TRUE~1),
           pred= inv.logit((1-fsm_app)*(b_fsm * (fsm - threshold_fsm) + plateau_fsm)
                           + fsm_app * plateau_fsm
                           + (1-hsm_app) * (b_hsm * (hsm - threshold_hsm) + plateau_hsm)
                           + hsm_app * plateau_hsm))
  gg <- fit.sp %>% 
    group_by(hsm) %>%
    summarise(mean=mean(pred),
              q95=quantile(pred,probs=0.95),
              q05=quantile(pred,probs=0.05)) %>%
    ggplot()+
    geom_line(aes(hsm,mean))+
    geom_ribbon(aes(hsm,ymin=q05,ymax=q95),alpha=0.2)+
    geom_point(data=db.sp.prob,aes(hsm,prob))+
    theme_bw()+
    xlab("hsm scaled")+
    ylab("Probability of presence")+
    ggtitle(paste0("Posterior and observed probability for ",sp,
                   "\nFSM fixed"))
  ggsave(paste0("figs3/",sp,"_hsmposterior.png"), plot = gg)
  
  }