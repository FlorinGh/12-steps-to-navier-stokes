# Step 10: 2D Poisson equation
# this is very similar to the Laplace eq, but this time we a have a non-homgenous
# eq, which has a free term that is a source term
# this will be expressed as two pressure spikes (one positive and one negative)
# the iterations will advance in pseudo-time to relax the spikes

import numpy as np
import pylab as pl
from pylab import cm
from mpl_toolkits.mplot3d import Axes3D
pl.ion()

# Variable declaration
nx = 100
ny = 100
nt = 1000

xmin = 0.0
xmax = 2.0
ymin = 0.0
ymax = 1.0

dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

# Initializing all matrices
p = np.zeros((nx,ny))
pn = np.zeros((nx,ny))
b = np.zeros((nx,ny))

x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)

# Defining the source
b[nx/4,ny/4] = -100
b[3*nx/4,3*ny/4] = -100
b[nx/4,3*ny/4] = -80
b[3*nx/4,ny/4] = -99


# The scheme is copied from step 9:
for n in range(nt):
    pn = p.copy()
    p[1:-1,1:-1] = (dx**2*(pn[1:-1,0:-2] + pn[1:-1,2:]) + \
    dy**2*(pn[0:-2,1:-1] + pn[2:,1:-1]) - 
    b[1:-1,1:-1]*dx**2*dy**2)/(2*dx**2+2*dy**2)
    
    p[0,:] = p[-1,:] = p[:,0] = p[:,-1] = 0.0
    
# For plotting we will use the function created in step9, B
def plot2D(x, y, p):
    fig = pl.figure(figsize = (11,7), dpi = 100)
    ax = Axes3D(fig)
    X, Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm, 
            linewidth=0, antialiased=True)
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.set_zlim(-0.015,0.015)    
    #ax.view_init(30,225)

plot2D(x, y, p)

