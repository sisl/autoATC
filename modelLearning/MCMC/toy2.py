
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pymc 
from numpy import array, empty
from numpy.random import randint

import networkx as nx

__all__ = ['state_origin', 'state_dest', 'frac', 'x_nodes', 'y_nodes',
           'xy_points', 'dataScore']



x_true = np.array([0., 0.,  10., 10., 10.,   0., -10.,  -5., -10.])
y_true = np.array([0., 10., 10., 0., -10., -10., -10.,   0., 10.])
Nnodes_true = len(x_true)

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

s_o = np.random.random_integers(0, Nnodes_true-1, Nsamples_data)    
s_d = s_o + 1;
s_d[s_d >= Nnodes_true] = 1



x_orig = x_true[s_o]; y_orig = y_true[s_o]
x_dest = x_true[s_d]; y_dest = y_true[s_d]

fracs_true = np.random.random_sample(Nsamples_data) 
x_meas = (x_dest - x_orig) * fracs_true + x_orig + np.random.randn(Nsamples_data)*0.25
y_meas = (y_dest - y_orig) * fracs_true + y_orig + np.random.randn(Nsamples_data)*0.25



xy_meas = np.column_stack([x_meas, y_meas]);


Nnodes= 7
Nsamples = 300

state_origin = np.array(range(Nnodes))

x_nodes = np.array([pymc.Uniform('x_nodes_%i'%i, lower=-15., upper=15.) for i in range(Nnodes)])
y_nodes = np.array([pymc.Uniform('y_nodes_%i'%i, lower=-15., upper=15.) for i in range(Nnodes)])




p_rewire = pymc.Uniform('p_rewire',lower=0., upper=1.)
@pymc.stochastic(dtype=int)
def adjMatrix(n=Nnodes, k=2, p=0.): #make p random?
    """connectendeness graph."""
    
    def random(n,k,p):
        G = nx.connected_watts_strogatz_graph(n, k, p)
        #returning array to avoid weirdness with stochastic shape        
        return nx.adjacency_matrix(G).toarray()
    
    
    def logp(value, n, k, p):
        if value.max() > 1 or value.min() < 0:
             return -np.inf
         
        G = nx.Graph(value)
        if not nx.is_connected(G):
            return -np.inf
        
        deg = nx.degree(G)
        return -(np.array(deg.values())**4).sum()

  
Nsamples_multi = Nsamples/Nnodes
# s_tp1 = np.array([pymc.Multinomial('multi_%i'%i, p=P_s_tp1[i], n=Nsamples_multi, plot=False) for i in range(Nnodes)])
#frac = np.linspace(0, 1, Nsamples_multi) 
frac = np.random.rand(Nsamples_multi)

tauNoise = 25 #pymc.Normal('sampleNoise', mu=25, tau=1/2.) #pymc.Uniform('sampleNoise', lower=5, upper=50);

x_noise = 0#pymc.Normal('x_noise', mu=0, tau=tauNoise, size=Nsamples_multi*Nnodes)
y_noise = 0#pymc.Normal('y_noise', mu=0, tau=tauNoise, size=Nsamples_multi*Nnodes)


@pymc.deterministic
def xy_points(s_o=state_origin, P=adjMatrix, f = frac,
              x_n=x_nodes, y_n=y_nodes, x_noise = x_noise, y_noise = y_noise): 
    #Note that we manually define the Dirichlet as a parent (Prows)
    #Somehow it doesn't recognize that!
    
    
    
    xy_m = np.empty((Nsamples_multi*Nnodes, 2))
    popidx = range(Nsamples_multi);
    
    for s_origin in s_o:
        x_o = x_n[s_origin]; 
        y_o = y_n[s_origin];
        
        s_d = np.array([d for d in range(Nnodes) if P[s_origin, d] != 0])
               
        if(len(s_d) == 0):
            s_d = np.array([s_origin])
            
        #compute measurements based on fractions
        kn = Nsamples_multi/len(s_d) * np.ones(len(s_d))
        kn[-1] += Nsamples_multi - kn.sum()
        
        #grab destination nodes
        s_d = np.concatenate([np.repeat(s,k) for (s,k) in zip(s_d, kn)])

        x_d = x_n[s_d]
        y_d = y_n[s_d]
        
        f = np.concatenate([np.random.rand(k) for k in kn])

        x_m = (x_d-x_o)*f + x_o #+ x_noise[popidx]
        y_m = (y_d-y_o)*f + y_o #+ y_noise[popidx]
    
        #populate output
        xy_m[popidx, :] = np.column_stack([x_m,y_m])
        popidx = range(popidx[-1]+1, popidx[-1]+1+Nsamples_multi)
    
    return xy_m
    
          

from scipy.misc import logsumexp
import scipy.spatial.distance as scipyDist
def scoreSimple(samples, data):
    d2 = scipyDist.cdist(samples, data, 'sqeuclidean')
    

    #only keep the closest. Maybe consider the 3 closest points?
    
    #return logsumexp(-d2) 
    #these are closest points to the data
    d2_data = np.min(d2, axis=0)
    sigma2_data = 1.
    
    #these are closest points to samples
    d2_samples = np.min(d2, axis=1)
    sigma2_samples = 20.
    #we want a smaller sigma2_samples since we don't expect samples
    #to just be out there on their own! We expect them to be close to the
    #data somehow. 
    
    
    #assuming a normal likelihood, loglikelihood is -val
    return -d2_data.sum()/sigma2_data - d2_samples.sum()/sigma2_samples







import scipy
Sigma_x = 10*np.eye(2);
Sigma_z = .1*np.eye(2);
T_xz = scipy.linalg.inv(Sigma_x + Sigma_z);
T_xx = scipy.linalg.inv(Sigma_x + Sigma_x);
T_zz = scipy.linalg.inv(Sigma_z + Sigma_z);

T_xz_sq = scipy.linalg.sqrtm(T_xz)
T_xx_sq = scipy.linalg.sqrtm(T_xx)
T_zz_sq = scipy.linalg.sqrtm(T_zz)

T_xz_sq_det = scipy.linalg.det(T_xz_sq)
T_xx_sq_det = scipy.linalg.det(T_xx_sq)
T_zz_sq_det = scipy.linalg.det(T_zz_sq)
def scoreKernel(x, z): #x are samples, z is measured data
    k = 2
    n = x.shape[0]
    m = z.shape[0]
    
    pi_f = (2*np.pi)**(k/2.)
        
          
    xt = np.dot(x, T_xz_sq)
    zt = np.dot(z, T_xz_sq)
    D_xz = scipyDist.cdist(xt,zt,'sqeuclidean'); D_xz *= -0.5;
    
    v = -np.log(np.exp(D_xz).sum()*T_xz_sq_det/n/m/pi_f)   
     
    
    xt = np.dot(x, T_xx_sq)
    D_xx = scipyDist.pdist(xt,'sqeuclidean'); D_xx *= -0.5;
     
    zt = np.dot(z, T_zz_sq)
    D_zz = scipyDist.pdist(zt,'sqeuclidean'); D_zz *= -0.5;
    v += 0.5*np.log((2*np.exp(D_xx).sum() + n) *T_xx_sq_det/n/n/pi_f)  + 0.5*np.log((2*np.exp(D_zz).sum()+m)*T_zz_sq_det/m/m/pi_f )
      
    return -v*500
    

@pymc.potential(verbose=True)
def dataScore(samples=xy_points):
    return scoreSimple(samples, xy_meas)
    #return scoreKernel(samples, xy_meas)
    
    




def sample(this):
    for p in this.extended_parents:
        p.random()

    if hasattr(this, 'random'):
        this.random()
    return this.value





def LDST(G, d=3, beta=1.):
    """ Create a PyMC Stochastic for a random lower degree Spanning Tree on
    base graph G
    Parameters
    ----------
    G : nx.Graph, base graph to span
    d : int, degree bound parameter
    beta : float, "inverse-temperature parameter" for depth bound
    """

    T = nx.minimum_spanning_tree(G)
    T.base_graph = G
    T.d = d

    @mc.stoch(dtype=nx.Graph)
    def ldst(value=T, beta=beta):
        return -beta * pl.sum(pl.array(T.degree().values()) >= d)

    return ldst

class STMetropolis(pymc.Metropolis):
    """ A PyMC Step Method that walks on spanning trees by adding a
    uniformly random edge not in the tree, removing a uniformly random
    edge from the cycle created, and keeping it with the appropriate
    Metropolis probability (no Hastings factor necessary, because the
    chain is reversible, right?)
    Parameters
    ----------
    stochastic : nx.Graph that is a tree and has a base_graph which it
                 spans
    """
    def __init__(self, stochastic):
        # Initialize superclass
        pymc.Metropolis.__init__(self, stochastic, scale=1., tally=False)

    def propose(self):
        
        if np.random.rand() < 0.2:
            self.stochastic.random()
        else:
            self.stochastic.value = nx.adjacency_matrix(nx.algorithms.double_edge_swap(nx.Graph(self.stochastic.value))).toarray()


    def reject(self):
        """ Restore the graph to its state before more recent edge swap"""
        self.rejected += 1
        for s in self.stochastics:
            s.value = s.last_value 
        
    @staticmethod
    def competence(s):
        return 0
        
        