# Step3: Diffusion only
# in this step only the diffusion term from the Navier-Stokes equations
# is used to determine the velocity;
# the diffusion term is a second order derivative multiplied by
# the diffusive coefficient, which in this case is the viscosity;
# to develop the liniar algebraic aproximation of the second order derivative,
# the Taylor expansion is used;

import numpy as np
import pylab as pl
pl.ion() # all function will be ploted in the same graph 

D = 2.0 # length of the domain (irrelevant in this problem)

nx = 41 # number of gridpoints
dx = D/(nx-1) # distance between any pair of adjacents grid points
grid = np.linspace(0,D,nx) # creating the space grid

nt = 100 # the number of timesteps we want to calculate (iterations)
nu = 0.3 # the value of viscosity
sigma = 0.2 # parameter sigma, related to Courant number
dt = sigma*dx**2/nu # the duration of each timestep

u = np.ones(nx) # initializing the a matrix for velocities
u[.5/dx : 1/dx+1] = 2 # input of initial conditions, same as step1
pl.figure(figsize = (11,7), dpi = 100)
pl.plot(grid, u)

un = np.ones(nx) # used only to initialize a matrix with the same dimension as u
# in here will be kept the velocity for current time step


for n in range(nt): # loop for time iteration
    un = u.copy() # copy the current values of velocity for each time step
    for i in range(1,nx-1): # loop over the entire space domain
        u[i] = un[i] + nu*dt/dx**2*(un[i+1] - 2*un[i] + un[i-1])
        # compute the velocities; basically the diffusion coefficient is 0.2
    pl.plot(grid, u) # plots all profiles on the same graph
    pl.ylim([1.,2.2])
    pl.xlabel('X')
    pl.ylabel('Velocity')    
    pl.title('1D Diffusion only')


# Discussion:
# the time step is calculated using a relaxation factor called Courant number;
# after running the program it can be seen that the wave does not move along
# the x axis;
# but the profile changes in a way that the velocities between 0.5 and 1 are
# decreasing and the velocities in the vecinity are increasing
# given more time for diffusion, the velocities are reducing even more