#!/usr/bin/env python



import h5py

import os  


def extractGeo(fname, fout):
    try:
        fdata = h5py.File(fname);
        if 'PseudoTruth' in fdata:
            for enc in fdata['PseudoTruth'].values():
                fout.write('%.16f, %.16f\n'%(enc[2,0] , enc[3,0]))
    except:
        print '\tFailed to process ' + fname



fout = open('extractedPoints.csv', 'w')    

for fn in os.listdir('.'):
     if os.path.isfile(fn) and os.path.splitext(fn)[1] == '.mat':
        print 'processing ' + fn
        extractGeo(fn, fout)
        
fout.close()