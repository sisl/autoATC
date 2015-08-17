#!/usr/bin/env python


import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


    

def plotCSV(map, csvfile, latlims, lonlims):
    lats = []
    lons = [] 
    with open(csvfile, 'rb') as f:
        reader = csv.reader(f)
        cnt = 0
        for row in reader:
            lat = float(row[0])
            lon = float(row[1])
            
            if lon < lonlims[0] or lon > lonlims[1] or lat < latlims[0] or lat > latlims[1]:
                continue
            
            cnt += 1
            lats.append(lat)
            lons.append(lon)
    
    
    print cnt , ' total of encounters'
    
    x,y = map(lons, lats)
    map.plot(x, y, 'ko', markersize=1)
    
    
    minlat = min(lats)
    maxlat = max(lats)
    minlon = min(lons)
    maxlon = max(lons)
    return (minlat, maxlat, minlon, maxlon)


def plotAirports(map, minlat, maxlat, minlon, maxlon):
    lats = []
    lons = []
    airports = []
    is_small = []
    with open('airports.csv', 'rb') as f:
        reader = csv.reader(f)
        headers = reader.next()
        
        
        cnt = 0
        for row in reader:
            lat = float(row[4])
            lon = float(row[5])
            
            if lon > maxlon or lon < minlon or lat > maxlat or lat < minlat:
                continue
            
            type = row[2]
            region = row[9]
            
            if type == 'heliport':
                continue
            
            ident = row[1]
            cnt += 1
            
            lats.append(lat)
            lons.append(lon)
            airports.append(ident)
            is_small.append(type=='small_airport')
            
            if ident in ['K39N', 'KBLM', 'KDED', 'KEVB', 'KCPM', 'KEMT']:
                print '\'%s\': [%.16f , %.16f], ' %(ident, lat, lon)
                
            if cnt > 5000:
                break
            
    
    
    x,y = map(lons, lats)
    
    #small airports
    x_s = [xx for (xx,t) in zip(x, is_small) if t]
    x_b = [xx for (xx,t) in zip(x, is_small) if not t]
    
    #big airports
    y_s = [xx for (xx,t) in zip(y, is_small) if t]
    y_b = [xx for (xx,t) in zip(y, is_small) if not t]
    
    if len(x_s) > 0:
        map.plot(x_s, y_s, 'mo', markersize=30, alpha=.3)
    if len(x_b) > 0:
        map.plot(x_b, y_b, 'bo', markersize=30, alpha=.3)
    
    for name, xpt, ypt in zip(airports, x, y):
        plt.text(xpt+100, ypt+100, name, fontsize=10)
 
