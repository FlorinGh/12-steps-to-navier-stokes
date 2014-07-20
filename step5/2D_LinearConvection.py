# Step5: 2D Linear Convection
# the ecuation to be solved: du/dt + c*du/dx +c*du/dy = 0
# in this equation we have the acumulation and convection terms
# this time the problem is 2D and we want to solve for the velocities

# to prepare the scheme, finite differences are used, a forward diff for time
# and backward diff for space derivatives
# special functions are used for 2d visualisation

# the domain will be a 2 X 2 square
# Initial Conditions: velocity is 2 for x, y in [0.5, 1] and 1 elsewhere
# Boundary Conditions: velocity is 1 at the edges of the domain (x, y = 0, 2)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
pl.ion()

# Variable declaration
nx = 81 # initial value 81
ny = 81 # initial value 81
nt = 100 # number of time steps; initial value 100
c = 1

dx = 2.0 / (nx-1)
dy = 2.0 / (ny-1)
sigma = 0.2
dt = sigma*dx

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

u = np.ones((ny,nx)) # Create a 1xn vector of 1's
un = np.ones((ny,nx))

# Assign initial conditions
u[.5/dy:1/dy+1, .5/dx:1/dx+1] = 2 # This is a 2D hat function

# Plot initial condition
fig = pl.figure(figsize = (11,7), dpi = 100)
ax = Axes3D(fig)
X, Y = np.meshgrid(x,y)
surf1 = ax.plot_surface(X,Y,u[:])
pl.xlabel('X')
pl.ylabel('Y')
pl.title('2D Linear Convection: Initial condition')

# Iterating in two dimensions
# First using nested for loops - use this method delete the quotes
"""
for n in range(nt+1): # loop across number of time steps
    un = u.copy()
    for i in range(1,len(u)):
        for j in range(1, len(u)):
            u[i,j] = un[i,j] - c*dt/dx*(un[i,j]-un[i-1,j]) - c*dt/dy*(un[i,j]-un[i,j-1])
            u[0,:] = 1
            u[-1,:] = 1
            u[:,0] = 1
            u[:,-1] = 1
"""
# Second using array operations
for n in range(nt+1):
    un[:] = u[:]
    u[1:,1:] = un[1:,1:]-c*dt/dx*(un[1:,1:]-un[0:-1,1:])-c*dt/dy*(un[1:,1:]-un[1:,0:-1])
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

fig = pl.figure(figsize=(11,7), dpi=100)
ax = Axes3D(fig)
surf2 = ax.plot_surface(X, Y, u[:])
pl.xlabel('X')
pl.ylabel('Y')
pl.title('2D Linear Convection: Solution')
