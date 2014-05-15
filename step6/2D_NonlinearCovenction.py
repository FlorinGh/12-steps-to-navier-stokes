# Step6: 2D Nonlinear Convection
# the ecuations to be solved are: du/dt + u*du/dx +v*du/dy = 0
#                                 dv/dt + u*dv/dx +v*dv/dy = 0
# in this equation we have the acumulation and convection terms
# but this time the convection is nonlinear and the equations are coupled
# this time the problem is 2D and we want to solve for the velocities

# First we use the equations given in the IPython notebook
import numpy as np
import pylab as pl
pl.ion()

# Variable declaration
nx = 101
ny = 101
nt = 80
c = 1

dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.2
dt = sigma*dx

x = np.linspace(0.0,2.0,nx)
y = np.linspace(0.0,2.0,ny)
 
u = np.ones((ny,nx))
v = np.ones((ny,nx))
un = np.ones((ny,nx))
vn = np.ones((ny,nx))

# Initialize the solution
u[.5/dy:1/dy+1, .5/dx:1/dx+1] = 2
v[.5/dy:1/dy+1, .5/dx:1/dx+1] = 2

# Apply the scheme
for n in range(nt+1):
    un = u.copy()
    vn = v.copy()
    
    u[1:,1:] = un[1:,1:] - un[1:,1:] * dt/dx * (un[1:,1:]-un[0:-1,1:])\
    -vn[1:,1:] * dt/dy * (un[1:,1:]-un[1:,0:-1])
    v[1:,1:] = vn[1:,1:] - un[1:,1:] * dt/dx * (vn[1:,1:]-vn[0:-1,1:])\
    -vn[1:,1:] * dt/dy * (vn[1:,1:]-vn[1:,0:-1])
    
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1
    
    v[0,:] = 1
    v[-1,:] = 1
    v[:,0] = 1
    v[:,-1] = 1

from pylab import cm
fig = pl.figure(figsize=(11,7), dpi=200)
ax = fig.gca(projection = '3d')
X,Y = np.meshgrid(x,y)
ax.plot_surface(X,Y,v,cmap=cm.coolwarm)







