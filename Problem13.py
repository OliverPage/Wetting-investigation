import scipy as sp
import matplotlib
import matplotlib.pyplot as py
import math
import numpy as np

#Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'

# [Starter code for Problem 8]

# ------------------------------------------------------------
# Solving the bulk with an inhomogeneous iterative method
# ------------------------------------------------------------
def potential(y, epsilon_w=1):
    '''Calculates the potentail as according to equation 36. '''
    if y >= 1:
        return -1*epsilon_w * (y**-3)
    else:
        return 10**10

#Define some useful objects
class Latticesite:
    def __init__(self,siteNumber, epsilon_w=1):
        self.index = siteNumber
        self.coordinate = []
        self.NNs = []
        self.NNNs = []
        self.potential = potential(np.floor(siteNumber/Ly), epsilon_w)
        self.density_current = 0.2
        self.density_previous = 0.2
    
    def update(self):
        '''Update the density. Remember to account for divergences under iteration'''
        if self.density_current<1.0:
            self.density_previous = self.density_current
        elif self.density_current>1.0:
            self.density_previous = 1.0
        else:
            self.density_previous = 0.0

def iterate(sites,k,mu,beta,epsilon=1):
    '''Perform a single iteration of eq. 46 for the particle at site k, given the sitelist sites'''
    nDens = 0.0
    nnDens = 0.0
    for neighbor in sites[k].NNs:
        nDens = nDens + sites[neighbor].density_previous
    for neighbor in sites[k].NNNs:
        nnDens = nnDens + sites[neighbor].density_previous
    return 1 * (1-sites[k].density_previous)*sp.exp(beta*(mu+epsilon*nDens+0.25*epsilon*nnDens-sites[k].potential)) #Here we assume T,mu,V are all in units of the interaction strength 

#Initialize the lattice
def initialize(Lx,Ly, epsilon_w=1):

    def getIndex(Lx,nSites):
        '''Converts from coordinates (x,y) to lattice index'''
        return lambda x,y:(Lx*x + y)

    nSites = Lx*Ly
    sites = [Latticesite(k, epsilon_w) for k in range(nSites)]  
    pos = getIndex(Lx,Ly)   
    
    for site in sites:
        x,y = [site.index//Lx,site.index%Lx]                                #Get x and y coordinates and from those the coordinates of the neighbors
        nns = set([((x+1)%Lx,y),((x-1)%Lx,y),(x,(y+1)%Ly),(x,(y-1)%Ly)])    
        nnns = set([( (x+1)%Lx,(y+1)%Ly ),( (x-1)%Lx, (y-1)%Ly),((x-1)%Lx,(y+1)%Ly),((x+1)%Lx,(y-1)%Ly)])
        site.NNs = [pos(x[0],x[1]) for x in nns]           #Store the neighbor indices as instance variables
        site.NNNs = [pos(x[0],x[1]) for x in nnns]
        site.density_previous = 0.2                        #Initialize the system in the low density limit
        if (y>Ly/3 and y<(2/3)*Ly):
            if (x<Lx/3 and x>0):
                site.density_previous = 0.8
    return sites
    

#Now we iterate the solver until the density is converged
def run(mu,T,Lx,Ly,cTol=10**-8,mixing=0.1,iterMax=1000,show=True,epsilon=1,epsilon_w=1, savename=None):
    'Calculates the density profile at a given mu,T'
    sites = initialize(Lx,Ly,epsilon_w)
    convergence = 0.1
    iteration = 0.0
    while (convergence>cTol) and (iteration<iterMax):
        for k in range(len(sites)):
            sites[k].density_current = (sum([site.density_previous for site in sites])/sum([site.density_current for site in sites])) * (sites[k].density_previous*(1-mixing) + mixing*iterate(sites,k,mu,1/T,epsilon)) #Calculate new state of the system from the old state of the system
        two_norm = sum([(site.density_current-site.density_previous)**2 for site in sites])
        convergence = math.sqrt(two_norm)
        iteration = iteration +1
        for site in sites:
            site.update()
    
    'Can then return an image of the density profile'
    z = []
    for site in sites:
        z.append(site.density_previous)
    Z = sp.array(z)
    Z = sp.reshape(Z,(Lx,Ly))
    x,y = sp.meshgrid(range(Lx),range(Ly))
    contourlevels = py.MaxNLocator(nbins=10).tick_values(Z.min(), Z.max())
    cmap = py.get_cmap('PuBu')
    norm = matplotlib.colors.BoundaryNorm(contourlevels, ncolors=cmap.N, clip=True)
    if show==True:
        print(sites[5].density_previous)
        fig,ax = py.subplots()
        ax.set_xlabel("$x (\sigma )$")
        ax.set_ylabel("$y (\sigma )$")
        ax.plot()
        raw = ax.pcolormesh(x,y,Z,cmap=cmap,norm=norm)
        fig.colorbar(raw)
        if savename!=None: 
            py.savefig(savename+'.pdf')
            py.title(savename)
        py.show()
    else:
        #print(sites[5].density_previous)
        return py.pcolormesh(x,y,Z,cmap=cmap,norm=norm)


#The lattice is a (Lx x Ly square lattice)
Lx = 30
Ly = 30

Beta, Mu =1, -2.53*1.2
Iter=10


run(Mu,1/1,Lx,Ly,epsilon=1.2,epsilon_w=0.5/Beta,iterMax=100, savename="Prob13Drop")
run(Mu,1/1,Lx,Ly,epsilon=1.2,epsilon_w=2/Beta,iterMax=400, savename="Prob13Wetting")
"""
for i in range(21):
    run(-2.53,1/Beta,Lx,Ly,epsilon=1.2,epsilon_w=0.4*i/Beta,iterMax=Iter,savename="$Prob13: ε_w=$"+str(np.round(0.4*i/Beta,1)))

for i in range(21):
     run(-2.53,1/Beta,Lx,Ly,epsilon=1.2,epsilon_w=2/Beta,iterMax=i*2,savename="Prob13: Iterations="+str(2*i))
"""