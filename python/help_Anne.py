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
metfile = 'C:/Users/anne.baranger/Documents/Era5-land/data/ERA5-land/swcd-1950-2021-layer3.nc'

### Nom du fichier de sortie
filecsv = 'C:/Users/anne.baranger/Documents/Era5-land/data/ERA5-land/swcd-1950-2021-layer3_min.csv'

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
dstart = get_start_time(netMet)+timedelta(hours=t[0].data.item())  
#get end date
dend = get_start_time(netMet) + timedelta(hours=t[-1].data.item()) 
#get timestep
deltat = t[1].data.item()-t[0].data.item()

dates = np.arange(np.datetime64(dstart), np.datetime64(dend), timedelta(hours=deltat ))

#get a list of years
year_list = dates.astype('datetime64[Y]').astype(int)+1970
years=np.unique(year_list,return_index=True)

varmet = netMet.variables['swvl3']

#create yearly min raster
YearMinRast = np.zeros((len(years[0]),varmet.shape[1],varmet.shape[2]))

#loop to n-1 year
prevind=0
for i,year in enumerate(years[0][:-1]):
    YearMinRast[i]=np.minimum.reduce(varmet[prevind:years[1][i+1]])
    prevind=years[1][i+1]
#add last year
YearMinRast[-1]=np.minimum.reduce(varmet[prevind:])

#for raster in range(YearMinRast.shape[0]):
#    plotraster(YearMinRast[raster])


rastMin = np.ones((varmet.shape[1],varmet.shape[2]))*varmet.missing_value

#calculation of 5 min by point
globalmax = np.maximum.reduce(YearMinRast)
#get coordinate of good values
inds = np.argwhere(globalmax!=varmet.missing_value)
##calculation of 5 min by point
for pt in inds:
    rastMin[pt[0],pt[1]] = np.mean(np.sort(YearMinRast[:,pt[0],pt[1]])[0:5])
    

np.savetxt(filecsv,rastMin,delimiter=",")
