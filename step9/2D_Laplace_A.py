# Step 9: 2D Laplace Equation
# this is the homogenous Poisson equation
# the equation to be solved: d2p/dx2 + d2p/dy2 = 0
# as it can be seen it is steady state equations and it appears in may applications
# this code will sole the equation using loops in the style of other steps

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from pylab import cm
pl.ion()

# Variable declaration
nx = 51  # number of x nodes
ny = 51  # number of y nodes
nt = 1000 # maximum number of iterations

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

# Plot the initial conditions
fig = pl.figure(figsize = (11,7), dpi = 100)
ax = Axes3D(fig)
X, Y = np.meshgrid(x,y)
surf1 = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm, 
linewidth=0, antialiased=True)
ax.view_init(30,225)

# Applying the scheme to determine the pressure distributions
l1norm_target = 0.01
for n in range(nt+1):
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
    if l1norm < l1norm_target:
        break

# Plot the results
fig = pl.figure(figsize = (11,7), dpi = 100)
ax = Axes3D(fig)
surf2 = ax.plot_surface(X,Y, p[:], rstride=1, cstride=1, cmap=cm.coolwarm, 
linewidth=0, antialiased=True)
ax.view_init(30,225)
