#!/usr/bin/env python


import csv

 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.


state = 'CA'

if state == 'CA':
    lonlims = [-124., -112.]
    latlims = [33., 41.]
elif state == 'FL':
    lonlims = [-88., -79.]
    latlims = [24., 31.] 
elif state == 'NY':
    lonlims = [-78, -70]
    latlims = [38, 43]    
    
map = Basemap(projection='merc',lat_0=37,lon_0=-119,resolution='l',
                  llcrnrlon=lonlims[0], llcrnrlat=latlims[0],
                  urcrnrlon=lonlims[1], urcrnrlat=latlims[1]) #area_thresh = 0.1?


# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
#map.drawcountries(linewidth=0.25)
map.drawstates(linewidth=0.25)
map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30))
# make up some data on a regular lat/lon grid.


lats = []
lons = [] 
with open('extractedPoints.csv', 'rb') as f:
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


print cnt , ' total of ecounters'

x,y = map(lons, lats)
map.plot(x, y, 'ko', markersize=1)


minlat = min(lats)
maxlat = max(lats)
minlon = min(lons)
maxlon = max(lons)


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

x_s = [xx for (xx,t) in zip(x, is_small) if t]
x_b = [xx for (xx,t) in zip(x, is_small) if not t]

y_s = [xx for (xx,t) in zip(y, is_small) if t]
y_b = [xx for (xx,t) in zip(y, is_small) if not t]

map.plot(x_s, y_s, 'mo', markersize=30, alpha=.3)
map.plot(x_b, y_b, 'bo', markersize=30, alpha=.3)

for name, xpt, ypt in zip(airports, x, y):
    plt.text(xpt+100, ypt+100, name, fontsize=10)
                 
# labels = ['Sitka', 'Baranof Warm Springs', 'Port Alexander']
# for label, xpt, ypt in zip(labels, x, y):
#     plt.text(xpt+10000, ypt+5000, label)
# plt.title('contour lines over filled continent background')

map.plot([x_s[0], x_s[-1]], [y_s[0], y_s[-1]])
 
#plt.show()
