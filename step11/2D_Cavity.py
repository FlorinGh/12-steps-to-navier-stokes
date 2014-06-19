# Step 11: 2D Cavity flow
# this code will solve the Navier Stokes equations in 2D
# with different boundary conditions different problems can be solved
# in this problem a cavity flow will be solve
# specific BC impose a moving wall of 1 m/s on the upper side of the square
# all other sides are no slip walls with zero velocity and constant pressure
# details of the problem description can be seen here:
# http://nbviewer.ipython.org/github/barbagroup/CFDPython/blob/master/lessons/15_Step_11.ipynb

import numpy as np
import pylab as pl
from pylab import cm
pl.ion()

# Variable declaration
xmin = 0.0
xmax = 2.0
ymin = 0.0
ymax = 2.0
nt = 10000
c = 1


# Defining the call function for the problem
def cavityFlow(xmin, xmax, ymin, ymax, nt, c):
    # Number of nodes
    nx = 51
    ny = 51
    
    # Cell size
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)
    
    # Step size
    dt = dx/100
    
    # The grid
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)
    Y, X = np.meshgrid(y,x)
    
    # Fluid properties
    rho = 1
    nu = 0.1
    
    # Matrix initialization
    u = np.zeros((nx, ny))
    v = np.zeros((nx, ny))
    p = np.zeros((nx, ny))
    
    u, v, p = cavitySolver(nt, u, v, dt, dx, dy, p, rho, nu)
    ContourPlot2D(u, v, p, Y, X)


# Defining a function for the pressure build-up
def buildUpB(b, rho, dt, u, v, dx, dy):
    b[1:-1,1:-1] = rho*(1/dt*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx) + \
    (v[1:-1,2:]-v[1:-1,0:-2])/(2*dy)) - \
    ((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx))**2 - \
    2*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dy))*((v[2:,1:-1]-v[0:-2,1:-1])/(2*dx)) - \
    ((v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))**2)
    
    return b


# Defining the pressure-Poisson solution:
def pressPoisson(p, dx, dy, b):
    pn = np.empty_like(p)
    pn = p.copy()
    nit = 50
    
    for r in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = 0.5*((pn[2:,1:-1]+pn[0:-2,1:-1])*dy**2 + \
        (pn[1:-1,2:]+pn[1:-1,0:-2])*dx**2)/(dx**2+dy**2) - \
        0.5*(dx*dy)**2/(dx**2+dy**2)*b[1:-1,1:-1]
        
        p[-1,:] = p[-2,:]        # dp/dx = 0 at x = 2
        p[0,:] = p[1,:]          # dp/dx = 0 at x = 0
        p[:,0] = p[:,1]          # dp/dy = 0 at y = 0
        p[:,-1] = 0              # p = 0 at y = 2
    
    return p
    
    
# Defining the solver function:
def cavitySolver(nt, u, v, dt, dx, dy, p, rho, nu):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.empty_like(p)
    
    for s in range(nt):
        un = u.copy()
        vn = v.copy()
        
        b = buildUpB(b, rho, dt, u, v, dx, dy)
        p = pressPoisson(p, dx, dy, b)
        
        u[1:-1,1:-1] = un[1:-1,1:-1] - \
        un[1:-1,1:-1]*(dt/dx)*(un[1:-1,1:-1]-un[0:-2,1:-1]) - \
        vn[1:-1,1:-1]*(dt/dy)*(un[1:-1,1:-1]-un[1:-1,0:-2]) - \
        (dt/(2*rho*dx))*(p[2:,1:-1]-p[0:-2,1:-1]) + \
        nu*(dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]) + \
        dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2]))
        
        v[1:-1,1:-1] = vn[1:-1,1:-1] - \
        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[0:-2,1:-1]) - \
        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[1:-1,0:-2]) - \
        dt/(2*rho*dy)*(p[1:-1,2:]-p[1:-1,0:-2]) + \
        nu*(dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1]) + \
        dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2]))
        
        u[0,:] = 0
        u[:,0] = 0
        u[:,-1] = c
        v[0,:] = 0
        v[:,0] = 0
        v[:,-1] = 0
        v[-1,:] = 0
        
    return u, v, p


# Defining the plotting function
def ContourPlot2D(u, v, p, Y, X):
    pl.figure(figsize = (11,7), dpi = 100)
    pl.contourf(X,Y,p,alpha=0.5,cmap=cm.gist_heat)# plotting the pressure field contours
    pl.colorbar()
    pl.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) # plotting velocity vectors
    pl.xlabel('X')
    pl.ylabel('Y')
    pl.title('Pressure contours and velocity vectors')

# Run the function
cavityFlow(xmin, xmax, ymin, ymax, nt, c)