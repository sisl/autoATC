
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as scipyDist
import pymc 
from numpy import array, empty
from numpy.random import randint


__all__ = ['state_origin', 'state_dest', 'frac', 'x_nodes', 'y_nodes',
           'xy_points', 'dataScore']



def seeP(P_array):
    Pleft = np.row_stack([p.value for p in P_array])
    p = 1-np.sum(Pleft, axis=1)
    return np.column_stack([Pleft, p])


x_true = np.array([-5., 0., -10., 10., 4.])
y_true = np.array([0., 10., 14., 7., -5.])

P_true = np.array( [[0., .4, .2]])

Nsamples_data = 400

# s_o = np.random.random_integers(0,1, Nsamples_data)
# 
# s_d = np.array(s_o)
# for i in range(2):
#     idx = s_o == i
#     d = 1;
#     if i == 1:
#       d = 3
#     s_d[idx] =  np.random.random_integers(i+1, d, np.sum(idx))    

s_o = np.random.random_integers(0, len(x_true)-1, Nsamples_data)    
s_d = np.random.random_integers(0, 1, Nsamples_data)

x_orig = x_true[s_o]; y_orig = y_true[s_o]
x_dest = x_true[s_d]; y_dest = y_true[s_d]

fracs_true = np.random.random_sample(Nsamples_data) 
x_meas = (x_dest - x_orig) * fracs_true + x_orig + np.random.randn(Nsamples_data)*0.25
y_meas = (y_dest - y_orig) * fracs_true + y_orig + np.random.randn(Nsamples_data)*0.25
xy_meas = np.column_stack([x_meas, y_meas]);


Nnodes= 5
Nsamples = 300
#state_origin = pymc.DiscreteUniform('origin',    lower=0, upper=Nnodes-1, size=Nsamples)
state_origin = np.array(range(Nnodes)) #np.random.randint(low=0, high=Nnodes, size=Nsamples)

#state_dest = pymc.DiscreteUniform('destination', lower=0, upper=Nnodes-1, size=Nsamples)

#frac = pymc.Uniform('fraction', lower=0, upper=1, size=Nsamples)

x_nodes = np.array([pymc.Uniform('x_nodes_%i'%i, lower=-15., upper=15.) for i in range(Nnodes)])
y_nodes = np.array([pymc.Uniform('y_nodes_%i'%i, lower=-15., upper=15.) for i in range(Nnodes)])



# @pymc.deterministic(plot=False)
# def xy_points(o=state_origin, d=state_dest, f=frac, x_n = x_nodes, y_n = y_nodes):
#     x_n_s = x_n #np.sort(x_n)
#     y_n_s = y_n #np.sort(y_n)
#     x = (x_n_s[d]-x_n_s[o])*f + x_n_s[o]
#     y = (y_n_s[d]-y_n_s[o])*f + y_n_s[o]
#     out = np.column_stack([x, y])
#     return out


Prows =  np.empty(Nnodes, dtype=object)
for i in range(Nnodes):
    t = np.ones(Nnodes)*10; t[i] = 0.1
    Prows[i] = pymc.Dirichlet('Dir_%i'%i, theta=t)

@pymc.deterministic
def P_s_tp1(P=Prows, st=state_origin):
  #dd = Prows[st].value
  #p_last = np.array([1. - np.sum(dd)])
  #p = np.concatenate([dd, p_last])
  return np.array([np.concatenate([p.value, np.array([1.-sum(p.value)])]) for p in Prows[st]])
  
  
Nsamples_multi = Nsamples/Nnodes
# s_tp1 = np.array([pymc.Multinomial('multi_%i'%i, p=P_s_tp1[i], n=Nsamples_multi, plot=False) for i in range(Nnodes)])
frac = np.random.rand(Nsamples_multi)

@pymc.deterministic
def xy_points(s_o=state_origin, f=frac,
              x_n=x_nodes, y_n=y_nodes, 
              Prows=Prows): 
    #Note that we manually define the Dirichlet as a parent (Prows)
    #Somehow it doesn't recognize that!
    
    xy_m = np.empty((Nsamples_multi*Nnodes, 2))
    popidx = range(Nsamples_multi);
    
    for s_origin in s_o:
        x_o = x_n[s_origin]; 
        y_o = y_n[s_origin];
        
#         s_d = [ [s for i in range(v)] for (s,v) in enumerate(s_t1[s_origin])]
#         #turn into numpy array
#         s_d = np.array([s for ll in s_d for s in ll])

        counts = (Prows[s_origin] * Nsamples_multi).astype(int)
        
        #counts = np.round(np.round(Prows[s_origin], 3) * Nsamples_multi).astype(int)        
        counts = np.append(counts, Nsamples_multi - counts.sum())
        
        s_d = np.array([s for (s, c) in enumerate(counts) for i in range(c)]) 
        
        #grab destination nodes
        x_d = x_n[s_d]
        y_d = y_n[s_d]
        
        #compute measurements based on fractions
        x_m = (x_d-x_o)*f + x_o
        y_m = (y_d-y_o)*f + y_o
    
        #populate output
        xy_m[popidx, :] = np.column_stack([x_m,y_m])
        popidx = range(popidx[-1]+1, popidx[-1]+1+Nsamples_multi)
#         import ipdb
#         ipdb.set_trace()

    
    return xy_m
    
    
      


def score(samples, data):
    d2 = scipyDist.cdist(samples, data, 'sqeuclidean')
    #only keep the closest. Maybe consider the 3 closest points?
    d2_0 = np.min(d2, axis=0)
    d2_1 = np.min(d2, axis=1)
    #assuming a normal likelihood, loglikelihood is -val
    return -np.sum(d2_0) - np.sum(d2_1)

@pymc.potential(verbose=True)
def dataScore(samples=xy_points):
    return score(samples, xy_meas)
    




def sample(this):
    for p in this.extended_parents:
        p.random()

    if hasattr(this, 'random'):
        this.random()
    return this.value

