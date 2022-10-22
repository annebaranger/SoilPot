
psi_min=psihorday_100 %>% 
  mutate(across(.cols=matches(c("psi")),
                ~ cut(.,
                      breaks=c(-Inf, -50000,-25000, -10000, -5000,-3000,-1000,-500,-300,-200,-100,Inf), 
                      labels=c("<-50MPa","-50<Psi<-25MPa","-25<Psi<-10MPa","-10<Psi<-5MPa",
                               "-5<Psi<-3MPa","-3<Psi<-1MPa","-1<Psi<-0.5MPa",
                               "-0.5/-0.3","-0.3/-0.2","-0.2/-0.1","-0.1/0"))))  %>% 
  # dplyr::select(matches(c("^x$","^y$","psi"))) %>%
  dplyr::select(x,y,psi) %>%
  relocate(y,.after=x) %>% 
  #pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=psi)) +
  #geom_sf(data = occ.fagsy)+
  theme_bw() +
  theme(axis.title=element_blank(),
        legend.key.size = unit(0.5,"cm"))+
  labs(fill="Potentiel (MPa)")+
  #facet_wrap(~name)+
  scale_fill_brewer(palette="RdYlBu")+
  coord_equal()

safety.margins.day$tmin %>% 
  mutate(t.min=scale(t.min,scale=TRUE,center=TRUE)) %>% 
  ggplot(aes(x=x,y=y,fill=t.min))+
  geom_tile()+
  scale_fill_gradientn(name = "Sensitivity of performance to x ",colours = turbo(6),na.value="transparent")+
  labs(fill="Performance")+
  theme_bw() +
  theme(axis.title=element_blank())+
  coord_equal()



safety.margins.day$tmin %>% 
  mutate(t.min=cut(t.min,
                   breaks=c(-Inf,-20,-15,-10,-5,0,5,10, Inf),
                   labels=c("<-20", "-20<Tmin<-15", "-15/-10","-10/-5","-5/0","0/5","5/10",">10"))) %>% 
  ggplot(aes(x=x,y=y,fill=t.min))+
  geom_tile()+
  labs(fill="Min temperature (Celsius)")+
  theme_bw() +
  theme(axis.title=element_blank(),
        legend.key.size = unit(0.5,"cm"))+
  scale_fill_brewer(palette="RdYlBu")+
  coord_equal()
ggsave(psi_min,
       filename = paste0("psi_min2.png"),
       path="figs/",
       device="png",
       scale=2)

safety.margins.day$tmin %>% 
    mutate(across(.cols=colnames(.)[c(-1,-2,-3)],
                  ~ cut(.,
                        breaks=c(-Inf,0,3,7,10,15,20,30, Inf),
                        labels=c("<0", "0<FSM<3", "3/7","7/10","10/15","15/20","20/30",">30"))))  %>% 
    # dplyr::select(matches(c("^x$","^y$","psi"))) %>%
    dplyr::select(-t.min) %>%
    pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
    filter(name %in% c("Abal","Bepe","Fasy","Piab","Psme","Quro")) %>% 
    mutate(name=as.factor(name)) %>% 
    ggplot()+
    geom_tile(aes(x=x,y=y,fill=value)) +
    theme_bw() +
    # theme(axis.title=element_blank(),
    #       legend.key.size = unit(0.5,"cm"))+
    # labs(fill="Potentiel (MPa)")+
    facet_wrap(~name)+
    scale_fill_brewer(palette="RdYlBu")+
    coord_equal()

## load chelsa and compute wai


psi_min %>% 
  select(x,y,psi) %>% 
  mutate(psi=cut(psi,
                 breaks=c(-Inf, -50000,-25000, -10000, -5000,-3000,-1000,-500,-300,-200,-100,Inf), 
                 labels=c("<-50MPa","-50<Psi<-25MPa","-25<Psi<-10MPa","-10<Psi<-5MPa",
                          "-5<Psi<-3MPa","-3<Psi<-1MPa","-1<Psi<-0.5MPa",
                          "-0.5/-0.3","-0.3/-0.2","-0.2/-0.1","-0.1/0"))) %>% 
  #pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=psi)) +
  #geom_sf(data = occ.fagsy)+
  theme_bw() +
  theme(axis.title=element_blank(),
        legend.key.size = unit(0.5,"cm"))+
  labs(fill="Potentiel (MPa)")+
  #facet_wrap(~name)+
  scale_fill_brewer(palette="RdYlBu")+
  coord_equal()



psi_h1 %>% 
  filter(is.na(psi_w)) %>% 
  filter(x<11) %>% 
  filter(x>8) %>% 
  filter(y<48) %>% 
  filter(y>43)
  ggplot()+
  geom_tile(aes(x=x,y=y),fill="red")+
  geom_sf(data=europe,fill=NA)+
  coord_sf()

  

#test plot LT50_spring for different LT50_winter and burburst date

df.LT50.t=data.frame(LT50_winter=seq(-60,-20,1)) %>% 
  mutate(spring_80=-5+30*(5+LT50_winter)/(2*80),
         spring_90=-5+30*(5+LT50_winter)/(2*90),
         spring_100=-5+30*(5+LT50_winter)/(2*100),
         spring_110=-5+30*(5+LT50_winter)/(2*110))
df.LT50.t %>% 
  pivot_longer(cols=matches("spring")) %>% 
  ggplot(aes(LT50_winter,value,color=name))+
  geom_point()+
  theme_bw()
       

as.data.frame(rast.fdg.sd,xy=TRUE) %>% 
  mutate(std=cut(std,
                 breaks=c(0,10,20,30,60,Inf))) %>% 
  ggplot(aes(x=x,y=y,fill=std))+
  geom_tile()
