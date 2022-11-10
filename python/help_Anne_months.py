#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Script pour extraire des valeurs d'un NETCDF multilayer

@author: Sylvain DUPIRE, LESSEM, INRAE 2022
Pour disposer de l'environnement Python qui va bien il faut installer miniconda : https://docs.conda.io/en/latest/miniconda.html
Une fois l'installation terminée il faut créer un environnement  que l'on peut appelé pyclim avec la commande suivante :

conda create -n pyclim python=3.9 numpy netCDF4

"""

import os,sys
from datetime import datetime, timedelta
from calendar import monthrange
import numpy as np
from netCDF4 import Dataset,MFDataset
import matplotlib.pyplot as plt


### Nom du fichier netcf
metfile = 'C:/Users/anne.baranger/Documents/Era5-land/data/CHELSA/CHELSA_EUR11_tasmin_day_2005_V1.1.nc'


### Nom du fichier de sortie

filecsv = 'C:/Users/anne.baranger/Documents/Era5-land/data/CHELSA/CHELSA_EUR11_tasmin_day_2005_V1.1.nc'



def plotraster(raster):
    plt.imshow(raster)
    plt.colorbar()
    plt.show()



def get_start_time(ncdf_dat):
    t=ncdf_dat.variables['time']
    t_split = t.units.split(' ')
    day_split = t_split[2].split('-')
    h_split = t_split[3].split(':')
    d0 = datetime(int(day_split[0]),
                  int(day_split[1]),
                  int(day_split[2]),
                  int(h_split[0]))
    return d0   



#open file in read only mode
netMet = Dataset(metfile)

#get time vector
t=netMet.variables['time']

#get start date
dstart = get_start_time(netMet)+timedelta(days=t[0].data.item())  

#get end date
dend = get_start_time(netMet) + timedelta(days=t[-1].data.item()) 

#get timestep
deltat = t[1].data.item()-t[0].data.item()

dates = np.arange(np.datetime64(dstart), np.datetime64(dend), timedelta(days=deltat ))

#get a list of months
month_list = dates.astype('datetime64[M]').astype(int)-(2005-1970)*12
months=np.unique(month_list,return_index=True)


varmet = netMet.variables['tasmin']


#create yearly min raster
MonthMinRast = np.zeros((len(months[0]),varmet.shape[1],varmet.shape[2]))


#loop to n-1 months
prevind=0
for i in months[0][:-1]:
    MonthMinRast[i]=np.minimum.reduce(varmet[prevind:months[1][i+1]])
    prevind=months[1][i+1]

#add last year
MonthMinRast[-1]=np.minimum.reduce(varmet[prevind:])

#Write new netcdf
g = Dataset('monthlymin.nc', 'w') # w


for attname in netMet.ncattrs():
   setattr(g, attname, getattr(netMet, attname))

g.createDimension('time', 12)
g.createDimension('lon', netMet.dimensions['lon'].size)
g.createDimension('lat', netMet.dimensions['lat'].size)

#Create variables

# Time
time = g.createVariable('time', netMet.variables["time"].dtype, netMet.variables["time"].dimensions)

#Proceed to copy the variable attributes
for attname in netMet.variables["time"].ncattrs():
   setattr(time, attname, getattr(netMet.variables["time"], attname))
#Finally copy the variable data to the new created variable
time [: ] = months[0] +1


# Tasmin
tasmin = g.createVariable('tasmin', netMet.variables["tasmin"].dtype,netMet.variables["tasmin"].dimensions,fill_value=getattr(netMet.variables['tasmin'],"_FillValue"))
#Proceed to copy the variable attributes
for attname in netMet.variables["tasmin"].ncattrs():
   setattr(tasmin, attname, getattr(netMet.variables["tasmin"], attname))
setattr(tasmin, 'grid_mapping', getattr(netMet.variables["tasmin"], 'grid_mapping'))
setattr(tasmin, '_FillValue', getattr(netMet.variables["tasmin"], '_FillValue'))
setattr(tasmin, 'missing_value', getattr(netMet.variables["tasmin"], 'missing_value'))
setattr(tasmin, 'long_name', getattr(netMet.variables["tasmin"], 'long_name'))
setattr(tasmin, 'standard_name', getattr(netMet.variables["tasmin"], 'standard_name'))
setattr(tasmin, 'units', getattr(netMet.variables["tasmin"], 'units'))
#Finally copy the variable data to the new created variable
tasmin [: ] = MonthMinRast[: ]


# Crs
crs = g.createVariable('crs', netMet.variables["crs"].dtype, netMet.variables["crs"].dimensions)

#Proceed to copy the variable attributes
for attname in netMet.variables["crs"].ncattrs():
   setattr(crs, attname, getattr(netMet.variables["crs"], attname))

# Lon
lon = g.createVariable('lon', netMet.variables["lon"].dtype, netMet.variables["lon"].dimensions)

#Proceed to copy the variable attributes
for attname in netMet.variables["lon"].ncattrs():
   setattr(lon, attname, getattr(netMet.variables["lon"], attname))



# lat
lat = g.createVariable('lat', netMet.variables["lat"].dtype, netMet.variables["lat"].dimensions)

#Proceed to copy the variable attributes
for attname in netMet.variables["lat"].ncattrs():
   setattr(lat, attname, getattr(netMet.variables["lat"], attname))