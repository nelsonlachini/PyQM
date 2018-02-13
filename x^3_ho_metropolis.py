import math
import numpy as np
import matplotlib.pyplot as plt

'''
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\METHODS\\\\\\\\\\\\\\\\\\\\\\\\\\\
'''

def acc():
    return (float(sum(acceptance))/float(N*len(acceptance)))

def reset(x):
    for j in range(0,N):                
        x[j] = 0.0

def thermalize(x):
    for j in range(0,7*N_cor):                 #thermalize     
        update(x)

def update(x):                                  # refresh the path x
    flips = 0    
    for j in range(0 , N):
        #print '\nj = %2d || x[%d] = %3g' %(j,j,x[j])
        old_x = x[j]                            # save original value
        old_Sj = S(j,x)
        x[j] += np.random.uniform(-eps,eps)     # update x[j]
        dS = S(j,x) - old_Sj                    # change in action
        
        if ((dS>0) and (math.exp(-dS)<np.random.uniform(0,1))):       
            x[j] = old_x                        # restore old value
        else:
            flips += 1          
    acceptance.append(flips)                           
                                    
def S(j,x):                             # returns the relevant part of action
    jp = (j+1)%N                        # next site
    jm = (j-1)%N                        # previous site
    return (a*(x[j]**2)/2. + 
    x[j]*(x[j]-x[jp]-x[jm])/a)


def compute_G(x,n):                     # returns the average over the paths
    g = 0.0                             # in a specific time G(t) --> Gn    
    for j in range(0,N):
        g += (x[j]**3*x[(j+n)%N]**3)
    return g/N
    
def avg(G):
    # MC avg of G
    return np.sum(G,axis=0)/len(G)
    
def sdev(G):
    # std dev of G
    g = np.asarray(G)
    return (np.abs(avg(g**2)-avg(g)**2))**0.5  

def deltaE(G):
    # Delta E(t)
    avgG = avg(G)   
    adE = avgG[:-1]/avgG[1:]
    adE = np.abs(adE)
    adE = np.log(adE)
    return adE/a
    
def bin(G):
    G_binned = []
    for i in range(0, len(G), binsize):
        G_avg = 0.0     
        for j in range(0,binsize):            
            G_avg += G[i+j]
        G_binned.append(G_avg/binsize)
    return G_binned



def plot(x,y,sdev_x,sdev_y):
    plt.errorbar(x , y , xerr=sdev_x , yerr=sdev_y , marker='s',
    mfc='red', mec='black', ms=4, mew=1 , linestyle='-' )
    plt.plot(np.arange(N-1),np.ones(N-1) ,'r')
    plt.axis([-1, 6, 0, 2.5])
    plt.show()

def bootstrap(G):
    N_cf = len(G)
    G_bootstrap = []
    for i in range(0,N_cf):
        alpha = int(np.random.uniform(0,N_cf))
        G_bootstrap.append(G[alpha])
    return G_bootstrap

def bootstrap_deltaE(G,nbstrap=100): # Delta E + errors       
    bsE = []    
    for i in range(nbstrap):        
        g = bootstrap(G)
        bsE.append(deltaE(g))
    bsE = np.array(bsE)
    sdevE = sdev(bsE) 
    return sdevE

def MCaverage(x,G):
    reset(x)
    thermalize(x)
    for alpha in range(0,N_cf):        
        for j in range(0,N_cor):
            update(x)
        for n in range(0,N):
            G[alpha][n] = compute_G(x,n)            
     

'''
\\\\\\\\\\\\\\\\\\\\\\\\PARAMETERS\\\\\\\\\\\\\\\\\\\\\\\
'''

# set parameters:
N = 20                                  # path length
N_cf = 1000                             # number of path configurations
N_cor = 20                              # number of updates to uncorrelate
a = 0.5                                 # time lattice spacing
eps = 1.4                               # random parameter
acceptance = []
binsize = 10

# create arrays:
x = np.zeros((N,),np.float)
G = np.zeros((N_cf,N),np.float)

# do the simulation:
MCaverage(x,G)

dE = deltaE(G).tolist()
sdevE = bootstrap_deltaE(bin(G)).tolist()

print '\nAcceptance average ratio: %10g' % (acc())

plot(np.arange(N-1), dE , None, sdevE )

print "\n%2s %10s %10s" % ("t","Delta E(t)","error")
print 26*"-"
for i in range(len(dE)/2):
    print "%2d %10g %10g" % (i,dE[i],sdevE[i])

print '\nAverage deltaE: %10g' % (sum(dE)/float(len(dE)))

