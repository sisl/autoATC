#!/usr/bin/env python


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import h5py

from plotGeosHelpers import *


#Choose which state for the limits!
state = 'FL'
if state == 'CA':
    lonlims = [-124., -112.]
    latlims = [33., 41.]
elif state == 'FL':
    lonlims = [-88., -79.]
    latlims = [24., 31.] 
elif state == 'NY':
    lonlims = [-78, -70]
    latlims = [38, 43]    

plt.ion() 
plt.clf()

##########################################
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
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
##########################################


#def plotPath(map, encounter):
fin = h5py.File('/home/zouhair/Downloads/KDED.hd5')


# Each encounter scenario, whether in text or mat file, follows the same format. 
# This format is a matrix of 8 columns and a varying number of rows. The columns are as follows:
# 0 - Time (Seconds)
# 1 - Identifier (Unique to each aircraft in the scenario)
# 2 - Latitude (degrees)
# 3 - Longitude (degrees)
# 4 - Altitude (ft)
# 5 - North Velocity (meters/second)
# 6 - East Velocity  (meters/second)
# 7 - Vertical Velocity (feet/second)


encounters = fin['encounters']


minlat = minlon = minalt = np.Inf
maxlat = maxlon = maxalt = -np.Inf


cnt = 0
for enc in encounters.values():
    for acIdx in np.unique(enc[1,]):
        idx = enc[1,] == acIdx
        
        meanAlt = np.mean(enc[4,])
        if meanAlt < 1400: #TODO: Figure out how to account for airport altitude...
            #map takes lons, lats
            x, y = map(enc[3,idx], enc[2,idx])
            c = plt.cm.bone(meanAlt / 3000)
            map.plot(x, y, color=c, alpha = .5)
            cnt += 1
        
    minlat = min(minlat, min(enc[2,]))
    maxlat = max(maxlat, max(enc[2,]))
    minlon = min(minlon, min(enc[3,]))
    maxlon = max(maxlon, max(enc[3,]))
    minalt = min(minalt, min(enc[4,]))
    maxalt = max(maxalt, max(enc[4,]))

#(minlat, maxlat, minlon, maxlon) = plotCSV(map, 'extractedPoints.csv', latlims, lonlims)

plotAirports(map, minlat, maxlat, minlon, maxlon)
                 

plt.show()
