# Step7: 2D Diffusion only
# the ecuations to be solved are: du/dt = nu*d2u/dx2 +nu*d2u/dy2

# in this equation we have the acumulation and diffusion terms
# the problem is 2D and we want to solve for the velocities
# the expected phenomena is that the in the domain will tend to spread
# from the high values to the low values; thsi can be best seen in the plots

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from pylab import cm
pl.ion()

# Variable declaration
nx = 31
ny = 31
nt = 17
nu = 0.5

dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.25 # under relaxation factor
dt = sigma*dx*dy/nu

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

# Initializing the velocity
u = np.ones((ny,nx)) # Create a ny by nx matrix of 1's
# Change the nx and ny values to see how u looks and ...
# print u
un = np.ones((ny,nx))

# Assign initial conditions
u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2 # This is the hat function

# Plot the initial conditions
fig =pl.figure(figsize = (11,7), dpi = 100)
ax = Axes3D(fig)
X, Y = np.meshgrid(x,y)
surf1 = ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)

# Applying the scheme in a function of time
def diffuse(nt):
    u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2 # reinitialize for every call
    
    for n in range(nt+1):
        un = u.copy()
        u[1:-1,1:-1] = un[1:-1,1:-1] + nu*dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]) + nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
        
        u[0,:] = 1
        u[-1,:] = 1
        u[:,0] = 1
        u[:,-1] = 1
    
    fig = pl.figure(figsize = (11,7), dpi = 100)
    ax = Axes3D(fig)
    surf2 = ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax.set_zlim(1,2.5)
    pl.show()
    
diffuse(10)
