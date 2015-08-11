#!/usr/bin/env python

import h5py
import os  
import math
import time

def distance_nm(lat1, long1, lat2, long2):
 
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    dtheta = (long1-long2) * degrees_to_radians

         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(dtheta) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length. (3440 NM)
    return arc*3440



airportLocs = {
                'K39N': [40.3992004394999995 , -74.6588973999000132], 
                'KBLM': [40.1869010899999992 , -74.1249008199999935], 
                'KDED': [29.0669994400000000 , -81.2837982199999942], 
                'KEVB': [29.0557003021240234 , -80.9488983154296875],
                'KCPM': [33.8899993895999998 , -118.2440032960000025], 
                'KEMT': [34.0861015320000007 , -118.0350036619999941], 
            }


outFiles = {k : h5py.File(k + '.hd5', 'w') for k in airportLocs }

for f in outFiles.values():
    f.create_group('encounters')
    
#exit(0)    

def extractPaths(fname):
    try:
        fdata = h5py.File(fname);
    except:
        print '\tFailed to read', fname
        return 

    try:
        if 'PseudoTruth' in fdata:
            truthGrp = fdata['PseudoTruth']            
            for (enckey, enc) in truthGrp.iteritems():
                for (id, latlon) in airportLocs.iteritems():
                    dist = min(distance_nm(latlon[0], latlon[1], enc[2,0], enc[3,0]), distance_nm(latlon[0], latlon[1], enc[2,-1], enc[3,-1]))                     
                    if dist < 6: #close to an airport of interest!
                        truthGrp.copy(enckey, outFiles[id]['encounters'])  
                        break #only one airport for each trajectory
                    
        else:
            print '\t %s does not have truth!'%fname
                    
        
    except:
        print '\tFailed to process', fname                                 
        return
    


for fn in os.listdir('.'):
     if os.path.isfile(fn) and os.path.splitext(fn)[1] == '.mat':
        tstart = time.time()
        print 'processing ' + fn
        extractPaths(fn)
        
        # do stuff
        telapsed = time.time() - tstart

        print '\t (%.2fs)' % telapsed 

    
print 'Dumping files ... '
for f in outFiles.values():
    print f.filename , ' has ' , len(f['encounters']), ' encounters'
    f.close()    
