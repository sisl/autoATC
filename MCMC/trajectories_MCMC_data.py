import numpy as np


import scipy.io
data = scipy.io.loadmat("../4D_points_sim_v7.3.mat");
Xfull = data['XBig']

data = scipy.io.loadmat("../turningpoints_means_sim_v7.3.mat");
XYnodes_heuristic = data['Xnodes']

#xy_meas = Xfull[0:3000, 0:2]

idx = np.arange(Xfull.shape[0])
np.random.shuffle(idx)
xy_meas = Xfull[idx[0:3000], 0:2]