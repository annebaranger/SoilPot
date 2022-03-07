#### Plotting environmental/climatic data  ####
###############################################

# Plot_Raw <- function(Raw_df, Shp = NULL, Dates, Legend = "Blabla") {
#   Raw_df <- pivot_longer(data = psi_min, cols = colnames(psi_min)[c(-1, -2)]) 
#   Raw_plot <- ggplot() + # create a plot
#     geom_raster(data = Raw_df, aes(x = x, y = y, fill = value)) + # plot the raw data
#     facet_wrap(~Values) + # split raster layers up
#     theme_bw() +
#     labs(x = "Longitude", y = "Latitude") + # make plot more readable
#     scale_fill_gradientn(name = Legend, colours = inferno(100)) # add colour and legend
#   if (!is.null(Shp)) { # if a shape has been designated
#     Raw_plot <- Raw_plot + geom_polygon(data = Shp, aes(x = long, y = lat, group = group), colour = "black", fill = "NA") # add shape
#   }
#   return(Raw_plot)
# }
# 


##Comparison of Psi_min according to the type of water computation
psi=psi_min %>%
  select(x,y,Psi_w_min,Psi_m_min) %>%
  pivot_longer(cols=colnames(.)[c(-1,-2)]) %>%
  ggplot() + # create a plot
  geom_raster( aes(x = x, y = y, fill = -value)) + # plot the raw data
  facet_wrap(~name) + # split raster layers up
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") + # make plot more readable
  scale_fill_gradientn(name = "potential min",trans="log",colours = viridis(80))+#,0.8,0.85,0.87,0.9,0.92
  coord_quickmap()

##Comparison of Psi_min according to the type of water computation
ratio=psi_min %>%
  select(x,y,SWC_w_ratio,SWC_m_ratio) %>%
  pivot_longer(cols=colnames(.)[c(-1,-2)]) %>%
  ggplot() + # create a plot
  geom_raster( aes(x = x, y = y, fill = value)) + # plot the raw data
  facet_wrap(~name) + # split raster layers up
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") + # make plot more readable
  scale_fill_gradientn(name = "SWC ratio",trans="log",colours = viridis(80))+#,0.8,0.85,0.87,0.9,0.92
  coord_quickmap()

cowplot::plot_grid(psi,ratio,nrow = 2,align="v") 

##Comparison Psi_min for constent texture and different WC computation
psi_min %>%
  select(x,y,Psi_w_min,Psi_m_min) %>%
  mutate(Psi_w_min=cut(Psi_w_min, breaks=c(-Inf, -50, -25, -20,-15,-10,-5,-3,-2,-1,Inf), labels=c("<-50","-50<Psi<-25","-25<Psi<-20","-20<Psi<-15","-15<Psi<-10","-10<Psi<-5","-5/-3","-3/-2","-2/-1","-1/0")),
         Psi_m_min=cut(Psi_m_min,breaks=c(-Inf, -50, -25, -20,-15,-10,-5,-3,-2,-1,Inf), labels=c("<-50","-50<Psi<-25","-25<Psi<-20","-20<Psi<-15","-15<Psi<-10","-10<Psi<-5","-5/-3","-3/-2","-2/-1","-1/0"))) %>%
  pivot_longer(cols=colnames(.)[c(-1,-2)]) %>%
  ggplot() + # create a plot
  geom_raster( aes(x = x, y = y, fill = value)) + # plot the raw data
  facet_wrap(~name) + # split raster layers up
  theme_bw() +
  scale_fill_gradientn(colours=viridis(10))+
  labs(x = "Longitude", y = "Latitude") + # make plot more readable
  #scale_fill_gradientn(name = Legend, colours = inferno(100),values=c(0,0.5,0.5,0.8,0.9,0.95,0.975,0.98,0.99,1))+
  coord_quickmap()


### Making temporal graphs of SWC 
chronology_swc <- function(point_to_extract=data.frame(Luberon=c(5.387792,43.814892),
                                                     Landes= c(-0.438943,44.427567), #Landes
                                                     Maison= c(-1.844066,47.649440), #la maison
                                                     Chamrousse= c(5.860739,45.092439), #Chamrousse
                                                     GrandEst=c(5.569744,48.475034)), #Grand-Est
                         moy=5,
                         SWC,
                         method="Mean") {
  ## Format points to extract
  loc= colnames(point_to_extract)
  row.names(point_to_extract)=c("Lon","Lat")
  point_to_extract=t(point_to_extract)
  row.names(point_to_extract) <- NULL
  
  ## Format SWC_stack
  SWC_time=terra::rast(SWC,crs="epsg:4326")
  SWC_time=as.data.frame(extract(SWC_time,point_to_extract),
                           row.names = loc) %>% 
    tibble::rownames_to_column(var="Loc") %>% 
    pivot_longer(cols=colnames(.)[c(-1)]) %>% 
    separate(name,sep="\\.",c("drop","time")) %>% 
    select(-drop) %>% 
    mutate(time=as.numeric(time)) %>% 
    group_by(Loc) %>% 
    arrange(value) %>% 
    mutate(rnk=row_number()) %>% 
    mutate(maxima=dplyr::case_when(rnk<moy+1~"Min",
                                   rnk>(max(rnk)-moy)~"Max")) %>% 
    ungroup() %>% 
    arrange(Loc,time)
  SWC_time %>% 
    ggplot(aes(time,value,maxima))+
    geom_point(aes(color=maxima))+
    geom_line()+
    facet_wrap(~Loc)+
    labs(y=method)
}




