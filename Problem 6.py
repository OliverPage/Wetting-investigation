import numpy as np
import matplotlib.pyplot as plt

def plotter(x_axis, y_axis, Title = None, x_label="Average Number Density of Lattice", y_label="Function Value", y2_axis = None, SaveName=None, Zeros="Yes"):
    """
    Plotting function, can plot multiple functions on the same axis. 
    """       
    zeros = np.zeros(len(x_axis))
    if Zeros=="Yes": plt.plot(x_axis, zeros, color='r', label = "Zero") 
    if y2_axis != None:plt.plot(x_axis, y2_axis, label="βε = 2/3", color='g')
    
    plt.plot(x_axis, y_axis, Label = "βε = 1", color='b')
    #plt.title(Title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    if SaveName != None: plt.savefig(SaveName)
    plt.show()
    
    index1 = np.argmin(np.absolute(y_axis))
    print("The root for βε = 1 is {} ".format(x_axis[index1]))
    
    index1b = np.argmin(np.absolute(y2_axis))
    print("The root for βε = 2/3 is {} ".format(x_axis[index1b]))
    
    

def func_run(Gamma, point_no = 100, Title = None, x_label = None, y_label = None, Save_name=None):
    """
    Calculates the function values for different densities for multiple values of β, then plots it.
    """ 
    density = np.linspace(0, 1, 100) # Average number density of lattice can take any vaue from zero to one
    func_return_alpha1 = []
    func_return_alpha2 = []
    
    func = lambda rho, alpha, gamma: rho - (1-rho)*np.exp(alpha*gamma + 5*alpha*rho) # LHS of eq. 45

    for i in range(len(density)):
        func_return_alpha1.append(func(density[i], 1, Gamma))
        func_return_alpha2.append(func(density[i], 2/3, Gamma))
        
    plotter(density, func_return_alpha1, Title, y2_axis=func_return_alpha2, SaveName=Save_name)

func_run(-3, point_no = 100, Title = "μ/ε=−3.0") #, Save_name="6i.pdf")
func_run(-2.5, point_no = 100, Title = "μ/ε=−2.5" , Save_name="6ii.pdf")
func_run(-2, point_no = 100, Title = "μ/ε=−2.0") #, Save_name="6iii.pdf")

































