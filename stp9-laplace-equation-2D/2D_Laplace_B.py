# Step 9: 2D Laplace Equation
# this is the homogenous Poisson equation
# the equation to be solved: d2p/dx2 + d2p/dy2 = 0
# as it can be seen it is steady state equations and it appears in may applications
# this code will solve the equation using functions

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from pylab import cm
pl.ion()


# This is the function that will be called for ploting
def plot2D(x, y, p):
    fig = pl.figure(figsize = (11,7), dpi = 100)
    ax = Axes3D(fig)
    X, Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm, 
            linewidth=0, antialiased=True)
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)    
    ax.view_init(30,225)


# This is the function that will be called for solving the Laplace eq
# this contains the scheme used in the A version
def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    while l1norm > l1norm_target:
        pn = p.copy()
        p[1:-1,1:-1] = (dx**2*(pn[1:-1,0:-2] + pn[1:-1,2:]) + \
        dy**2*(pn[0:-2,1:-1] + pn[2:,1:-1]))/(2*dx**2+2*dy**2)
        
        p[0,0] = (dx**2*(pn[0,-1] + pn[0,1]) + \
        dy**2*(pn[-1,0] + pn[1,0]))/(2*dx**2+2*dy**2)
        
        p[-1,-1] = (dx**2*(pn[-1,-2] + pn[-1,0]) + \
        dy**2*(pn[-2,-1] + pn[0,-1]))/(2*dx**2+2*dy**2)
        
        p[:,0] = 0
        p[:,-1] = y
        p[0,:] = p[1,:]
        p[-1,:] = p[-2,:]
        
        l1norm = (np.sum(np.abs(p[:])-np.abs(pn[:])))/np.sum(np.abs(pn[:]))
    
    return p

# Solving a specific case
# Variable declaration
nx = 51  # number of x nodes
ny = 51  # number of y nodes

dx = 2.0/(nx-1)
dy = 1.0/(ny-1)

x = np.linspace(0,2,nx)
y = np.linspace(0,1,ny)

# Initialization
p = np.zeros((ny,nx))
pn = np.zeros((ny,nx))

# Assign initial conditions
p[:,0] = 0
p[:,-1] = y
p[0,:] = p[1,:]
p[-1,:] = p[-2,:]

# To solve this case, type the following:
# To plot the initial conditions:
# plot2D(x, y, p)

# To solve the case:
# laplace2d(p, y, dx, dy, 0.01)

# To plot the results:
# plot2D(x, y, p)

# It can also be created a function that only asks for the number of nodes,
# the IC and BC and plots directly the IC and results
# equation, the number of nodes 