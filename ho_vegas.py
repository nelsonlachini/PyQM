import vegas
import math
import matplotlib.pyplot as plt


#The objective is to compute the propagator for the harmonic oscillator using 
#a path integral estimative and compare it with the theoretical values. The path integral will be
#descretized and evaluated as a N-1 multidimensional integral through the VEGAS algorithm. 


T = 4.0                             #evaluation time : T = tf-ti
xi = 0.0                            #fixed path initial position
xf = xi                             #fixed path final position
xrang = 2.0

#Lattice

def gslt(nbin):
    N = 8                           #number of nodes: xi=x(ti) ,x0=x(t0) ,... ,xN-2=x(tN-2), xf=x(tN-1)
    a = float(T/N)                  #time lattice spacing    
    m = 1                           #mass in natural units
    norm = (m/(2*math.pi*a))**(N/2) #normalization factor for the harmonic oscillator
    limit = 5                       #limits of the N-1 dimensional integration: 
                                    #   WARNING: if too high can be a problem for VEGAS    
    nitT    = 30                     #number of iterations to 'thermalize' the lattice
    nevalT  = 1000                  #number of Monte Carlo integral iterations to 'thermalize'    
    nitC    = 10                    #number of iterations for computation
    nevalC  = 1000                  #number of Monte Carlo integral iterations for computation
    
    def hoV(y):                     #evaluate the harmonic potential at x
        return (y**2)/2.0           #  omega**2 in natural units
    
    def slat(X):                    #builds and returns the harmonic oscillator action
        dx = 0.0                    #   in discretized form using xi, xf and x[i]={0,x1,...,xN-2}
        Vx = 0.0        
        Xi = x[k]
        Xf = Xi
       
        for j in range(N-1):
            Vx += hoV(X[j])
            if j == N-2:         
                dx += (Xf-X[0])**2 + (Xf - X[j])**2 
                Vx += hoV(Xi)
            else:    
                dx += (X[j+1] - X[j])**2
      
        return math.exp(-((m/(2*a))*dx + a*Vx))*norm    
    
    
    delta = xrang/nbin    
    x = []
    prop = []
    err = []
    for k in range(nbin+1):
        x.append(k*delta)
      
        integ = vegas.Integrator((N-1)*[[-limit,limit]])    #makes the object that will run the Monte Carlo integration       
        integ(slat , nitn=nitT , neval=nevalT)              #first ntherm runs just to heat the integrator    
        result = integ(slat , nitn=nitC , neval=nevalC)     #final computation        
        prop.append(result.mean)
        err.append(result.sdev)
        
    return (x,prop,err)



#Theoretical
  
def gsth(nbin):
    E0 = 0.5       
    delta = xrang/nbin
    
    x = []
    prop = []    
    for i in range(nbin+1): 
        x.append(i*delta)
        ovlap_xE0 = math.exp(-(x[i]**2)/2)/((math.pi)**(0.25))
        prop.append((ovlap_xE0**2)*math.exp(-E0*T))           
    return x,prop

laresult = gslt(15)
thresult = gsth(15)

plt.errorbar(laresult[0] , laresult[1] , yerr=laresult[2] , xerr=None , marker='s',
         mfc='red', mec='black', ms=4, mew=1 , linestyle=' ' )      #lattice result
plt.errorbar(thresult[0], thresult[1] ,  yerr=None , xerr=None)     #theoretical result
plt.axis([-0.1, 2.1, -0.01, 0.08])
plt.show()
