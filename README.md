# Wetting-investigation
Third year computational lab, code used to model the wetting transition.
For the Python code, you will need to uncommment the relevant parts to obtain the desired outputs. 

The file titled "Problem12.py" contains Python 3 code for problems 8-12 inclusive. 
This means it can be used to numerically calulate equilibrium densities, plot 2D density maps, 1D density vs distance plots, and adsorption plots.

The file titled "Problem13.py" contains Pyhton 3 code for problem 13. 
It plots the evolution of a system of fluids next to a surface. 

Droplet_Iterations.gif shows show the final state of the droplet changes after N iterations, each frame corresponds to the final state where the number of iterations has been increased from 0 in the first frame to 500 in tha final frame. This was calculated for a 30x30 lattice, with βε_w = 2, μ/ε = -2.53.

Wall Interaction.gif shows the final states as the dimensionless wall interaction strength, βε_w, is increased in increments of 0.4 from 0.0 to 6.0. This uses 100 iterations (to save on compuatational time), on a 30x30 lattice, using μ/ε = -2.53.
Wall Interaction 2.gif does the same thing under the same parameters, but in steps of 0.2 from 0.0 to 4.0.
Wall Interaction 2.gif (in my opinion) is the better of the two animations.

Problem 6.py calculates the equilibrium density of a fluid via graphical methods.  
