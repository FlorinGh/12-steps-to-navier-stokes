# Step2: Nonlinear Convection
# in this step the convection term of the NS equation is solved in 1D
# this time the wave velocity is nonlinear and this is the case in NS equation

import numpy as np
import pylab as pl
pl.ion() # all functions will be ploted in the same graph (similar to Matlab hold on)

D = 2.0 # length of the 1D domain
T = 1.0 # total amount of time for the analysis

nx = 31 # number of grid points
dx = D/(nx-1) # distance between any pair of adjacent grid points
gridx = np.linspace(0,2,nx) # creating the space grid

nt = 50 # number of time iterations
dt = T/nt # duration of each timestep

u = np.ones(nx) # initializing the matrix for velocities
u[0.5/dx:1/dx+1] = 2.0 # input of initial conditions, same as step1

un = np.ones(nx) # used only to initialize a matrix with the same dimension
# in here will be kept the velocity for current time step

for n in range(nt): # loop for time iteration
    un = u.copy() # copy the current time step velocity
    for i in range(1,nx): # loop over the entire space domain
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1]) # compute the velocity for
        # this step for the entire grid based on the velocities in previuos
        # time step and grid points
    pl.plot(np.linspace(0,2,nx),u) # the position of this statements
    # makes the plot to get the curves for each time step
    
# Discussion:
# the wave has moved to the right, with a change in shape
# in the same time, the hight is decreased
# increasing the total time will make the scheme unstable if there not enough timesteps
# that is if the time step is incrreased the problem loose resolution

# as a good recomendation, the time step dt must be corelated with the element size, dx
# if the size of element, dx is decreased, so should be the timestep
# the time should be smaller than the time needed for the wave to pass one grid cell
# for this case, the max velocity is 2, therefore the time step should be less than
