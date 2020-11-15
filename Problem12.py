import scipy as sp
import matplotlib
import matplotlib.pyplot as py
import math
import numpy as np
from tqdm import tqdm 


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
    def __init__(self,siteNumber, epsilon_w=1, NoPot=False):
        self.index = siteNumber
        self.coordinate = []
        self.NNs = []
        self.NNNs = []
        self.potential = potential(np.floor(siteNumber/Ly), epsilon_w)
        if NoPot==True: self.potential = 0.0
        self.density_current = 0.0
        self.density_previous = 0.0
    
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
    return (1-sites[k].density_previous)*sp.exp(beta*(mu+epsilon*nDens+0.25*epsilon*nnDens-sites[k].potential)) #Here we assume T,mu,V are all in units of the interaction strength 

#Initialize the lattice
def initialize(Lx,Ly, epsilon_w=1, No_Pot=False):

    def getIndex(Lx,nSites):
        '''Converts from coordinates (x,y) to lattice index'''
        return lambda x,y:(Lx*x + y)

    nSites = Lx*Ly
    sites = [Latticesite(k, epsilon_w, NoPot=No_Pot) for k in range(nSites)]  
    pos = getIndex(Lx,Ly)      
    for site in sites:
        x,y = [site.index//Lx,site.index%Lx]                                #Get x and y coordinates and from those the coordinates of the neighbors
        nns = set([((x+1)%Lx,y),((x-1)%Lx,y),(x,(y+1)%Ly),(x,(y-1)%Ly)])    
        nnns = set([( (x+1)%Lx,(y+1)%Ly ),( (x-1)%Lx, (y-1)%Ly),((x-1)%Lx,(y+1)%Ly),((x+1)%Lx,(y-1)%Ly)])
        site.NNs = [pos(x[0],x[1]) for x in nns]           #Store the neighbor indices as instance variables
        site.NNNs = [pos(x[0],x[1]) for x in nnns]
        site.density_previous = 0.2                        #Initialize the system in the low density limit
    return sites

#Now we iterate the solver until the density is converged
def run(mu,T,Lx,Ly,cTol=10**-8,mixing=0.1,iterMax=500,show=True,epsilon=1,epsilon_w=1, savename=None, prob_10=False, NoPot=False):
    'Calculates the density profile at a given mu,T'
    sites = initialize(Lx,Ly,epsilon_w, No_Pot=NoPot)
    convergence = 0.1
    iteration = 0.0
    while (convergence>cTol) and (iteration<iterMax):
        for k in range(len(sites)):
            sites[k].density_current = sites[k].density_previous*(1-mixing) + mixing*iterate(sites,k,mu,1/T,epsilon) #Calculate new state of the system from the old state of the system
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
    cmap = py.get_cmap('hot')
    norm = matplotlib.colors.BoundaryNorm(contourlevels, ncolors=cmap.N, clip=True)
    if prob_10!=True:
        if show==True:
            print(sites[5].density_previous)
            fig,ax = py.subplots()
            ax.set_xlabel('x')
            ax.set_ylabel('y')
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
    
    if prob_10==True:
        x = np.linspace(0, Ly, Ly)
        Z = np.transpose(Z)     
        return Lx*Z[0]


#The lattice is a (Lx x Ly square lattice)
Lx = 40
Ly = 40


#Run a few examples
Prob_8 = None
Prob_9 = None
Prob_10 = None
Prob_11 = None
Prob_12 = None

if Prob_8 != None:
    for mu in [-3,-2.5,-2]:
        run(mu,1.0,Lx,Ly,NoPot=True)
        run(mu,1.5,Lx,Ly,NoPot=True)
#figs, ax = py.subplots()
#ax = [run(k,0.5,4,4,show=False) for k in [-1,-4]]

if Prob_9 != None:
    Beta, epsi, epsi_w = 1, 1.2, 1.6
    for mu in [-2.67,-2.53]:
        run(mu*epsi,1.0,Lx,Ly,epsilon=epsi/Beta,epsilon_w=epsi_w/Beta, savename="Problem 9 Mu="+str(mu))

if Prob_10 != None:
    Lx, Ly = 15, 15
    Beta, epsi, epsi_w = 1, 1.2, 1.6
    mu267 = run(-2.67*epsi,1/Beta,Lx,Ly,epsilon=epsi,epsilon_w=epsi_w,prob_10=True)
    mu253 = run(-2.53*epsi,1/Beta,Lx,Ly,epsilon=epsi,epsilon_w=epsi_w,prob_10=True)
    
    py.plot(np.linspace(0,Lx,Lx), mu267, color='b', label="μ/ε=-2.67")
    py.plot(np.linspace(0,Lx,Lx), mu253, color='g', label="μ/ε=-2.53") 
    py.xlabel("Perpendicular Distance (σ)")
    py.ylabel("Effective Number Density")
    py.legend()
    py.tight_layout()
    py.savefig("Problem_10.pdf")
    py.title("1D Projection of Density against Perpnedicular Distance")
    py.show()


if Prob_11 != None:
    Beta, Epislon, Epsilon_w1, Epsilon_w2 = 1, 1, 2, 0.5
    Mu_Co = -2.5*Epislon
    Mu = np.linspace(Mu_Co-0.5, Mu_Co, 10)
    
    Gamma_array_Ew1, Gamma_array_Ew2 = [], []
    for i in Mu:
        mu1 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w1,prob_10=True)
        mu2 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w1,prob_10=True, NoPot=True)
        Gamma1 = mu1-mu2
        Gamma_array_Ew1.append(np.sum(Gamma1))
        
        mu3 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w2,prob_10=True)
        mu4 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w2,prob_10=True, NoPot=True)
        Gamma2 = mu3-mu4
        Gamma_array_Ew2.append(np.sum(Gamma2))
        
    Rel_Mu = [i - Mu_Co for i in Mu]
    py.plot(Rel_Mu, Gamma_array_Ew1, label="$βε_w = 2.0$", color='r')
    py.xlabel("$β(μ-μ_{coex})$")
    py.ylabel("Adsorption, Γ")
    py.legend()
    py.savefig("Prob11_βε_w=2.pdf")
    py.title("Adsorption plot")
    py.show()
    
    py.plot(Rel_Mu, Gamma_array_Ew2, label="$βε_w = 0.5$", color='r')
    py.xlabel("$β(μ-μ_{coex})$")
    py.ylabel("Adsorption, Γ")
    py.legend()
    py.savefig("Prob11_βε_w=0_5.pdf")
    py.title("Adsorption plot")
    py.show()
    
    py.plot(Rel_Mu, Gamma_array_Ew2, label="$βε_w = 0.5$", color='r')
    py.plot(Rel_Mu, Gamma_array_Ew1, label="$βε_w = 2.0$", color='c')
    py.xlabel("$β(μ-μ_{coex})$")
    py.ylabel("Adsorption, Γ")
    py.legend()
    py.savefig("Prob11_Both.pdf")
    py.title("Adsorption plot")
    py.show()
    
# ----------------------------- Week 5 ---------------------------------

if Prob_12 != None:
    Beta = 1.2
    Epislon = 1
    Epsilon_w1, Epsilon_w2, Epsilon_w3, Epsilon_w4, Epsilon_w5 = 1.8/Beta, 1.7/Beta, 1.6/Beta, 1.5/Beta, 1.4/Beta
    Mu_Co = -2.5*Epislon
    Mu = Mu_Co - np.linspace(0.3,0,20)/Beta
    rel_mu = [Beta*(i-Mu_Co) for i in Mu]
    Gamma_array_Ew1, Gamma_array_Ew2, Gamma_array_Ew3, Gamma_array_Ew4, Gamma_array_Ew5 = [], [], [], [], []
    for i in tqdm (Mu, desc="Loading…"):
        mu1 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w1,prob_10=True)
        mu2 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w1,prob_10=True, NoPot=True)
        Gamma_array_Ew1.append(np.sum(abs(mu1-mu2)))
        
        mu3 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w2,prob_10=True)
        mu4 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w2,prob_10=True, NoPot=True)
        Gamma_array_Ew2.append(np.sum(mu3-mu4))
        
        mu5 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w3,prob_10=True)
        mu6 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w3,prob_10=True, NoPot=True)
        Gamma_array_Ew3.append(np.sum(mu5-mu6))
        
        mu7 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w4,prob_10=True)
        mu8 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w4,prob_10=True, NoPot=True)
        Gamma_array_Ew4.append(np.sum(mu7-mu8))
        
        mu9 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w5,prob_10=True)
        mu10 = run(i,1/Beta,Lx,Ly,epsilon=Epislon,epsilon_w=Epsilon_w5,prob_10=True, NoPot=True)
        Gamma_array_Ew5.append(np.sum(mu9-mu10))
        

    py.plot(rel_mu, Gamma_array_Ew1, label="$βε_w = 1.8$")
    py.plot(rel_mu, Gamma_array_Ew2, label="$βε_w = 1.7$")
    py.plot(rel_mu, Gamma_array_Ew3, label="$βε_w = 1.6$")
    py.plot(rel_mu, Gamma_array_Ew4, label="$βε_w = 1.5$")
    py.plot(rel_mu, Gamma_array_Ew5, label="$βε_w = 1.4$")
    py.legend()
    py.xlabel("$β(μ-μ_{coex})$")
    py.ylabel("$Adsorption, Γ$")
    py.title("Adsorption plot")
    py.savefig("Prob12.pdf")
    py.show()


