
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pymc 
import networkx as nx
from numpy import array, empty
from numpy.random import randint


from trajectories_MCMC_data import xy_meas, XYnodes_heuristic



def seeP(P_array, k=None):
    Pleft = np.row_stack([p.value for p in P_array])
    pright = 1-np.sum(Pleft, axis=1)
    P = np.column_stack([Pleft, pright])
    
    if k != None:
        for r in range(P.shape[0]):
            p = P[r, :]
            p_sort_idx = np.argsort(p)[::-1] #decreasing order of probabilities
            #Only keep the Vk highest probabilities
            if hasattr(k,'__len__'):
                kidx = k[r]
            else:
                kidx = k
                
            p[p_sort_idx[kidx:]] = 0. 
            #renomalize so that probabilities add up to 1!
            P[r, :] = p / p.sum();
            
        
    return P 



Nnodes= 18
Nsamples = 3000
#state_origin = pymc.DiscreteUniform('origin',    lower=0, upper=Nnodes-1, size=Nsamples)
state_origin = np.array(range(Nnodes)) #np.random.randint(low=0, high=Nnodes, size=Nsamples)

#state_dest = pymc.DiscreteUniform('destination', lower=0, upper=Nnodes-1, size=Nsamples)

#frac = pymc.Uniform('fraction', lower=0, upper=1, size=Nsamples)

xmin = np.min(xy_meas[:,0])*1.1;ymin = np.min(xy_meas[:,1])*1.1;
xmax = np.max(xy_meas[:,0])*1.1;ymax = np.max(xy_meas[:,1])*1.1;  

# x_nodes = np.array([pymc.Uniform('x_nodes_%i'%i, lower=xmin, upper=xmax) for i in range(Nnodes)])
# y_nodes = np.array([pymc.Uniform('y_nodes_%i'%i, lower=ymin, upper=ymax) for i in range(Nnodes)])

xy_nodes = np.empty(Nnodes, dtype=object)
for i in range(Nnodes):
    mu = XYnodes_heuristic[i, 0:2]
    Sigma = np.empty((2,2))
    Sigma[0, :] = XYnodes_heuristic[i, 2:]
    Sigma[1, :] = Sigma[0, ::-1] #symmetric
    Tau = np.linalg.inv(Sigma)
    xy_nodes[i] = pymc.MvNormal('xy_nodes_%i'%i, mu = mu, tau=Tau)



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
    t = np.ones(Nnodes)*10; t[i] = 0.5
    Prows[i] = pymc.Dirichlet('Dir_%i'%i, theta=t)


#Cardinality / Sparsity of the transition matrices
#Vk = pymc.DiscreteUniform('Sparsity', lower=1, upper=min(3, Nnodes), size=Nnodes)
Vk = pymc.Geometric('Sparsity', p = 0.8, size=Nnodes)


  
  
Nsamples_multi = Nsamples/Nnodes
# s_tp1 = np.array([pymc.Multinomial('multi_%i'%i, p=P_s_tp1[i], n=Nsamples_multi, plot=False) for i in range(Nnodes)])
#frac = np.linspace(0, 1, Nsamples_multi) 
frac = np.random.rand(Nsamples_multi)


@pymc.deterministic
def adjMatrix(Prows = Prows, Vk = Vk):
    P = np.empty((Nnodes, Nnodes))
    
    for (row, s_o) in enumerate(np.arange(Nnodes)):
        p = np.append(Prows[s_o], 1.-Prows[s_o].sum())
        p_sort_idx = np.argsort(p)[::-1] #decreasing order of probabilities
        #Only keep the Vk highest probabilities
        k = min(Nnodes, Vk[s_o])
        p[p_sort_idx[k:]] = 0. 
        #renomalize so that probabilities add up to 1!
        p /= p.sum();
        P[row, :] = p
        
    return P
    
    
@pymc.potential(verbose=1)
def adjMatrixPotential(P=adjMatrix):
    G = nx.Graph(P)
        
    degVals = np.array(nx.degree(G).values())
    
    if degVals.max() > 5:
        return -np.inf
    
    return -((degVals**2).sum())



tauNoise = 1/3000. #pymc.Normal('sampleNoise', mu=25, tau=1/2.) #pymc.Uniform('sampleNoise', lower=5, upper=50);

x_noise = pymc.Normal('x_noise', mu=0, tau=tauNoise, size=Nsamples_multi*Nnodes)
y_noise = pymc.Normal('y_noise', mu=0, tau=tauNoise, size=Nsamples_multi*Nnodes)

    
@pymc.deterministic
def xy_points(Vk=Vk, # f = frac,
              xy_n=xy_nodes,
              x_noise = x_noise, y_noise = y_noise, 
              P=adjMatrix): 
    
    xy_m = np.empty((Nsamples_multi*Nnodes, 2))
    popidx = range(Nsamples_multi);
    
    for (row, s_o) in enumerate(np.arange(Nnodes)):
        x_o = xy_n[s_o][0]; 
        y_o = xy_n[s_o][1];
        
        #renomalize so that probabilities add up to 1!
        p = P[row, :];
        
        k = min(Nnodes, Vk[s_o])
        p_sort_idx = np.argsort(p)[::-1] #decreasing order of probabilities

        counts = np.round(p * Nsamples_multi).astype(int)
        while counts.sum() > Nsamples_multi:
            counts[p_sort_idx[0]] -= 1;
            
        while counts.sum() < Nsamples_multi:
            counts[p_sort_idx[k-1]] += 1;
            
        #counts = np.round(np.round(Prows[s_o], 3) * Nsamples_multi).astype(int)        
        #counts = np.append(counts, Nsamples_multi - counts.sum())
        
        s_d = np.array([s for (s, c) in enumerate(counts) for i in range(c)]) 
        
        #grab destination nodes
        x_d = np.array([xy[0] for xy in xy_n[s_d]])
        y_d = np.array([xy[1] for xy in xy_n[s_d]])
        
        #
        f = np.concatenate([np.linspace(0., 1., (s_d == c).sum()) for c in range(Nnodes)])

        #compute measurements based on fractions
        x_m = (x_d-x_o)*f + x_o + x_noise[popidx]
        y_m = (y_d-y_o)*f + y_o + y_noise[popidx]
    
        #populate output
        xy_m[popidx, :] = np.column_stack([x_m,y_m])
        popidx = range(popidx[-1]+1, popidx[-1]+1+Nsamples_multi)
#         import ipdb
#         ipdb.set_trace()

    
    return xy_m
    
    
      

from scipy.misc import logsumexp
import scipy.spatial.distance as scipyDist
from sklearn.neighbors import KDTree, BallTree


sigma2_data = 1.
sigma2_samples = .1
def scoreSimple(samples, data):
    d2 = scipyDist.cdist(samples, data, 'sqeuclidean')
    

    #only keep the closest. Maybe consider the 3 closest points?
    
    #return logsumexp(-d2) 
    #these are closest points to the data
    d2_data = np.min(d2, axis=0)
    
    #these are closest points to samples
    d2_samples = np.min(d2, axis=1)
    #we want a smaller sigma2_samples since we don't expect samples
    #to just be out there on their own! We expect them to be close to the
    #data somehow. 
    
    
    #assuming a normal likelihood, loglikelihood is -val
    return -d2_data.sum()/sigma2_data - d2_samples.sum()/sigma2_samples


tree_z = KDTree(xy_meas)
def scoreSimpleFast(x,z): #x are samples (change), z is measure data (unchanged!)
    (d2_samples, tmp) = tree_z.query(x)
    d2_samples **= 2
    
    tree_d = KDTree(x)
    (d2_data, tmp) = tree_d.query(z)
    d2_data **= 2
        
    return -d2_data.sum()/sigma2_data - d2_samples.sum()/sigma2_samples





import scipy
Sigma_x = .01*np.eye(2);
Sigma_z = .01*np.eye(2);
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
    
    v = -2*(np.exp(D_xz).sum()*T_xz_sq_det/n/m/pi_f)   
     
    
    xt = np.dot(x, T_xx_sq)
    D_xx = scipyDist.pdist(xt,'sqeuclidean'); D_xx *= -0.5;
     
    zt = np.dot(z, T_zz_sq)
    D_zz = scipyDist.pdist(zt,'sqeuclidean'); D_zz *= -0.5;
    v += ((2*np.exp(D_xx).sum() + n) *T_xx_sq_det/n/n/pi_f + (2*np.exp(D_zz).sum()+m)*T_zz_sq_det/m/m/pi_f )
      
    return -v*1000
    

#     
# xmin = np.min(xy_meas[:,0]);ymin = np.min(xy_meas[:,1]);
# xmax = np.max(xy_meas[:,0]);ymax = np.max(xy_meas[:,1]);  
# Ng = 100    
# x_g = np.linspace(xmin*1.2, xmax*1.2, Ng)
# y_g = np.linspace(xmin*1.2, xmax*1.2, Ng)
# xx, yy = np.meshgrid(x_g, y_g)
# xy_grid = np.column_stack([xx.flatten(), yy.flatten()])
# dz2 = (x_g[1] - x_g[0])*(y_g[1] - y_g[0])
# 
# def KLscore(x, z, returnQlog = False, qlog = None, qlogExp = None): #x are samples (change), z is measure data (unchanged!)
#     
#     tree_d = BallTree(x)
#     plog = tree_d.kernel_density(xy_grid, h=.5, kernel='gaussian',return_log=True, rtol=.1)
#     
#     computeQlog = (qlog == None)
#     if computeQlog: 
#         tree_m = BallTree(z)
#         qlog = tree_m.kernel_density(xy_grid, h=1, kernel='gaussian',return_log=True)
#         qlogExp = np.exp(qlog)
#     
#     
#     delta_log = (qlog - plog)
#     v  = (qlogExp      * delta_log).sum()
#     v -= (np.exp(plog) * delta_log).sum() #- since we want p * (log(p)-log(q))
#     
#     v *= -0.5 * dz2; #-sign to make it work as a similarity
#     
#     if returnQlog:
#         return v, qlog
#     else:
#         return v
# 
# 
# (tmp, qlog_meas) = KLscore(xy_meas, xy_meas, returnQlog = True)
# qlogExp_meas = np.exp(qlog_meas)
        
@pymc.potential(verbose=1)
def dataScore(samples=xy_points):
    return scoreSimpleFast(samples, xy_meas)
    #return scoreKernel(samples, xy_meas)

    #return KLscore(samples, xy_meas, qlog = qlog_meas, qlogExp= qlogExp_meas)
    


def sample(this):
    for p in this.extended_parents:
        p.random()

    if hasattr(this, 'random'):
        this.random()
    return this.value

