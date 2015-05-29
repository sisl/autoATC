
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pymc 
import networkx as nx
from numpy import array, empty
from numpy.random import randint


__all__ = ['state_origin', 'state_dest', 'frac', 'x_nodes', 'y_nodes',
           'xy_points', 'dataScore']



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
    t = np.ones(Nnodes)*10; t[i] = 0.5
    Prows[i] = pymc.Dirichlet('Dir_%i'%i, theta=t)


#Cardinality / Sparsity of the transition matrices
Vk = pymc.DiscreteUniform('Sparsity', lower=1, upper=min(4, Nnodes), size=Nnodes)

#Vk = pymc.Geometric('Sparsity', p = 0.7, size=Nnodes)


  
  
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
    
    if degVals.max() > 4:
        return -np.inf
    
    return -((degVals**2).sum())
    
@pymc.deterministic
def xy_points(Vk=Vk, # f = frac,
              x_n=x_nodes, y_n=y_nodes, 
              P=adjMatrix): 
    
    xy_m = np.empty((Nsamples_multi*Nnodes, 2))
    popidx = range(Nsamples_multi);
    
    for (row, s_o) in enumerate(np.arange(Nnodes)):
        x_o = x_n[s_o]; 
        y_o = y_n[s_o];
        
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
        x_d = x_n[s_d]
        y_d = y_n[s_d]
        
        #
        f = np.concatenate([np.linspace(0., 1., (s_d == c).sum()) for c in range(Nnodes)])

        #compute measurements based on fractions
        x_m = (x_d-x_o)*f + x_o
        y_m = (y_d-y_o)*f + y_o
    
        #populate output
        xy_m[popidx, :] = np.column_stack([x_m,y_m])
        popidx = range(popidx[-1]+1, popidx[-1]+1+Nsamples_multi)
#         import ipdb
#         ipdb.set_trace()

    
    return xy_m
    
    
      

from scipy.misc import logsumexp
import scipy.spatial.distance as scipyDist
from sklearn.neighbors import KDTree, BallTree

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


tree_z = KDTree(xy_meas)
def scoreSimpleFast(x,z): #x are samples (change), z is measure data (unchanged!)
    (d2_samples, tmp) = tree_z.query(x)
    d2_samples **= 2
    sigma2_data = 1.
    
    tree_d = KDTree(x)
    (d2_data, tmp) = tree_d.query(z)
    d2_data **= 2
    sigma2_samples = 20.
        
    return -d2_data.sum()/sigma2_data - d2_samples.sum()/sigma2_samples




import scipy

Sigma_x = .0001*np.eye(2);
Sigma_z = .0001*np.eye(2);
T_xz = scipy.linalg.inv(Sigma_x + Sigma_z);
T_xx = scipy.linalg.inv(Sigma_x + Sigma_x);
T_zz = scipy.linalg.inv(Sigma_z + Sigma_z);

T_xz_sq = scipy.linalg.sqrtm(T_xz)
T_xx_sq = scipy.linalg.sqrtm(T_xx)
T_zz_sq = scipy.linalg.sqrtm(T_zz)

T_xz_sq_det = scipy.linalg.det(T_xz_sq)
T_xx_sq_det = scipy.linalg.det(T_xx_sq)
T_zz_sq_det = scipy.linalg.det(T_zz_sq)


kD = 2
pi_f = (2*np.pi)**(kD/2.)


xt = np.dot(xy_meas, T_xx_sq)
n = xt.shape[0]
D_xx = scipyDist.pdist(xt,'sqeuclidean'); D_xx *= -0.5;
vxx = 0.5*np.log((2*np.exp(D_xx).sum() + n) *T_xx_sq_det/n/n/pi_f) 

def scoreKernel(x, z): #x are samples (change), z is measure data (unchanged!)
    n = x.shape[0]
    m = z.shape[0]
        
          
    xt = np.dot(x, T_xz_sq)
    zt = np.dot(z, T_xz_sq)
    D_xz = scipyDist.cdist(xt,zt,'sqeuclidean'); D_xz *= -0.5;
    
    v = -np.log(np.exp(D_xz).sum()*T_xz_sq_det/n/m/pi_f)   
     
        
    zt = np.dot(z, T_zz_sq)
    D_zz = scipyDist.pdist(zt,'sqeuclidean'); D_zz *= -0.5;
    v += vxx + 0.5*np.log((2*np.exp(D_zz).sum()+m)*T_zz_sq_det/m/m/pi_f )
      
    return -v
    
    
    

    
xmin = np.min(xy_meas[:,0]);ymin = np.min(xy_meas[:,1]);
xmax = np.max(xy_meas[:,0]);ymax = np.max(xy_meas[:,1]);  
Ng = 100    
x_g = np.linspace(xmin*1.2, xmax*1.2, Ng)
y_g = np.linspace(xmin*1.2, xmax*1.2, Ng)
xx, yy = np.meshgrid(x_g, y_g)
xy_grid = np.column_stack([xx.flatten(), yy.flatten()])
dz2 = (x_g[1] - x_g[0])*(y_g[1] - y_g[0])

def KLscore(x, z, returnQlog = False, qlog = None, qlogExp = None): #x are samples (change), z is measure data (unchanged!)
    
    tree_d = BallTree(x)
    plog = tree_d.kernel_density(xy_grid, h=.5, kernel='gaussian',return_log=True, rtol=.1)
    
    computeQlog = (qlog == None)
    if computeQlog: 
        tree_m = BallTree(z)
        qlog = tree_m.kernel_density(xy_grid, h=1, kernel='gaussian',return_log=True)
        qlogExp = np.exp(qlog)
    
    
    delta_log = (qlog - plog)
    v  = (qlogExp      * delta_log).sum()
    v -= (np.exp(plog) * delta_log).sum() #- since we want p * (log(p)-log(q))
    
    v *= -0.5 * dz2; #-sign to make it work as a similarity
    
    if returnQlog:
        return v, qlog
    else:
        return v


(tmp, qlog_meas) = KLscore(xy_meas, xy_meas, returnQlog = True)
qlogExp_meas = np.exp(qlog_meas)
        
@pymc.potential(verbose=1)
def dataScore(samples=xy_points):
    return scoreSimple(samples, xy_meas)
    #return scoreKernel(samples, xy_meas)

    #return KLscore(samples, xy_meas, qlog = qlog_meas, qlogExp= qlogExp_meas)
    
    




def sample(this):
    for p in this.extended_parents:
        p.random()

    if hasattr(this, 'random'):
        this.random()
    return this.value

