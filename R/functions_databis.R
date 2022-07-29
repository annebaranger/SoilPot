
psihor_100 %>% 
  mutate(across(.cols=matches(c("psi")),
                ~ cut(.,
                      breaks=c(-Inf, -20000, -10000, -5000, -3000, -1000, -500, -300, -200, -100, Inf),
                      labels=c("<-20MPa", "-20<Psi<-10MPa", "-10<Psi<-5MPa", "-5<Psi<-3MPa","-3<Psi<-1MPa",
                               "-1<Psi<-0.5MPa","-0.5/-0.3","-0.3/-0.2","-0.2/-0.1","-0.1/0"))))  %>% 
  # dplyr::select(matches(c("^x$","^y$","psi"))) %>%
  dplyr::select(x,y,psi) %>%
  relocate(y,.after=x) %>% 
  pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=value)) +
  geom_sf(data = occ.fagsy)+
  theme_bw() +
  theme(axis.title=element_blank(),
        legend.key.size = unit(0.5,"cm"))+
  labs(fill="Potentiel (MPa)")+
  facet_wrap(~name)+
  scale_fill_brewer(palette="RdYlBu")+
  coord_sf()

occ.fagsy = db.tree %>%  
  filter(ESPAR=="09") %>% 
  dplyr::select(IDP) %>% 
  distinct() %>% 
  left_join(db.stand,by="IDP") %>% 
  dplyr::select(IDP,XL,YL) %>% 
  distinct()
occ.fagsy=st_as_sf(occ.fagsy, coords = c("XL", "YL"), crs = 2154)
occ.fagsy=st_transform(occ.fagsy,crs="epsg:4326")

occ.fagsy=cbind(occ.fagsy,extract(rast(psihor_100,crs="epsg:4326"),vect(occ.fagsy)))


occ.fagsy %>% 
  ggplot(aes(abs(psi)))+
  geom_density()+
  scale_x_log10()+
  geom_vline(xintercept = 3150,color="red")+
  geom_vline(xintercept = 3556.795, color="green")
plot(vect(occ.fagsy))


# list of targetted species in fni data
tree.target=data.frame(species=c("quro","fagsy","abab","pial"),
                       fni.key=c("02","09","61","62"),
                       p50=c(6880,3150,3790,5810))

for (i in 1:dim(tree.target)[1]){
  occ.tree = db.tree %>%  
    filter(ESPAR==tree.target$fni.key[i]) %>% 
    dplyr::select(IDP) %>% 
    distinct() %>% 
    left_join(db.stand,by="IDP") %>% 
    dplyr::select(IDP,XL,YL) %>% 
    distinct()
  occ.tree=st_as_sf(occ.tree, coords = c("XL", "YL"), crs = 2154)
  occ.tree=st_transform(occ.tree,crs="epsg:4326")
  
  occ.tree=cbind(occ.tree,extract(rast(psihor_100,crs="epsg:4326"),vect(occ.tree)))
  
  q05=quantile(occ.tree$psi,probs=0.05,na.rm=TRUE)
  plot=occ.tree %>%
    ggplot()+
    geom_density(aes(abs(psi)))+
    geom_vline(xintercept = tree.target$p50[i],colour="red")+
    geom_vline(xintercept = -q05,colour="green")+
    scale_x_log10()
  print(plot)
  print(q05+tree.target$p50[i])
}



### safety margins
psimin[,1:6] %>% 
  pivot_longer(cols=colnames(.)[c(-1,-2)]) %>% 
  filter(value>-10000) %>% 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=value))+
  facet_wrap(~name)


## look for correlations between tmin and psimin
as.data.frame(c(rast(psimin[,1:3],crs="epsg:4326"),
  resample(rast(tmin[,1:3],crs="epsg:4326"),rast(psimin[,1:3],crs="epsg:4326"))),xy=TRUE) %>% 
  filter(psi>(-20000)) %>% 
  ggplot(aes(x=psi,y=t.min)) +
  geom_point(alpha=0.3)


## buil