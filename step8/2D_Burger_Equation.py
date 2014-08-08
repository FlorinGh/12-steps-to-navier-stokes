# Step 8: Burgers' equations
# it's an equation that has both the convective and diffusive terms
# the ecuations to be solved are:
# du/dt + u*du/dx + v*du/dy = nu*d2u/dx2 +nu*d2u/dy2
# this a very complex equation but has analitycal solutions
# these equations can have discontinuous solutions, that is shocks

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
pl.ion()

# Variable declaration
nx = 51
ny = 51
nt = 1200

dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.0009
nu = 0.005
dt = sigma*dx*dy/nu

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

# Initialization
u = np.ones((ny,nx))
v = np.ones((ny,nx))
un = np.ones((ny,nx))
vn = np.ones((ny,nx))

# Assign initial conditions
u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2
v[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2

# Plot the initial condition
fig  = pl.figure(figsize = (11,7), dpi = 100)
ax = Axes3D(fig)
X, Y = np.meshgrid(x,y)
ax.plot_wireframe(X, Y, u[:])
pl.xlabel('X')
pl.ylabel('Y')
pl.title('Burgers Equation: Initial condition')
ax.set_zlim(1,2.0)

# The scheme
for n in range(nt+1):
    un = u.copy()
    vn = v.copy()
    
    u[1:-1,1:-1] = un[1:-1,1:-1] - \
    dt/dx*un[1:-1,1:-1]*(un[1:-1,1:-1]-un[0:-2,1:-1]) - \
    dt/dy*vn[1:-1,1:-1]*(un[1:-1,1:-1]-un[1:-1,0:-2]) + \
    nu*dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]) + \
    nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
    
    v[1:-1,1:-1] = vn[1:-1,1:-1] - \
    dt/dx*un[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[0:-2,1:-1]) - \
    dt/dy*vn[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[1:-1,0:-2]) + \
    nu*dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1]) + \
    nu*dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])
    
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1
        
    v[0,:] = 1
    v[-1,:] = 1
    v[:,0] = 1
    v[:,-1] = 1
    
# Plot the results
fig = pl.figure(figsize = (11,7), dpi = 100)
ax = Axes3D(fig)
X, Y = np.meshgrid(x,y)
ax.plot_wireframe(X,Y,u[:])
pl.xlabel('X')
pl.ylabel('Y')
pl.title('Burgers Equation: Solution after 1200 steps')
ax.set_zlim(1,2.0)