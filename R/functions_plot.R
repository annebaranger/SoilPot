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




