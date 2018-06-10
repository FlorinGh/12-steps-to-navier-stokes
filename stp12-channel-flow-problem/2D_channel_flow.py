# Step 12: 2D Channel Flow
# this another problem solved with the Navier Stokes Eq
# this time the domain is a channel
# and in the U momentum equation we have a source term which makes the fluid
# to move form left to right
# introducing a periodic boundary condition

import numpy as np
import pylab as pl
from pylab import cm
pl.ion()

# The problem will be solved by modules, using functions


# Defininig the function for pressure build-up
def buildUpB(rho, dt, dx, dy, u, v):
    b = np.empty_like(u)
    
    # The same body as in step 11
    b[1:-1,1:-1] = rho*(1/dt*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx) + \
    (v[1:-1,2:]-v[1:-1,0:-2])/(2*dy)) - \
    ((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx))**2 - \
    2*(u[1:-1,2:]-u[1:-1,0:-2])/(2*dy)*(v[2:,1:-1]-v[0:-2,1:-1])/(2*dx) - \
    ((v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))**2)
    
    # Introducing the periodic boundary condition for x = 0
    # the idea is that the first value in b (b[0,1:-1]) must be computed 
    # with the values one place before of all arrays in the upper formula
    # [2:,] -> [1,]; [1:-1,] -> [0,]; [0:-2,] -> [-1,]
    b[0,1:-1] = rho*(1/dt*((u[1,1:-1]-u[-1,1:-1])/(2*dx) + \
    (v[0,2:]-v[0,0:-2])/(2*dy)) - \
    ((u[1,1:-1]-u[-1,1:-1])/(2*dx))**2 - \
    2*(u[0,2:]-u[0,0:-2])/(2*dy)*(v[1,1:-1]-v[-1,1:-1])/(2*dx) - \
    ((v[0,2:]-v[0,0:-2])/(2*dy))**2)    
    
    # Introducing the periodic boundary condition for x = 2
    # the idea is that the last value in b (b[-1,1:-1]) must be computed 
    # with the next values of all arrays in the upper formula
    # [2:,] -> [0,]; [1:-1,] -> [-1,]; [0:-2,] -> [-2,]
    b[-1,1:-1] = rho*(1/dt*((u[0,1:-1]-u[-2,1:-1])/(2*dx) + \
    (v[-1,2:]-v[-1,0:-2])/(2*dy)) - \
    ((u[0,1:-1]-u[-2,1:-1])/(2*dx))**2 - \
    2*(u[-1,2:]-u[-1,0:-2])/(2*dy)*(v[0,1:-1]-v[-2,1:-1])/(2*dx) - \
    ((v[-1,2:]-v[-1,0:-2])/(2*dy))**2)
    
    return b


# Defining the pressure-Poisson function:
def pressPoisson(p, dx, dy, b):
    pn = np.empty_like(p)
    pn = p.copy()
    nit = 50
    
    # The same body as in step 11
    for r in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = 0.5*((pn[2:,1:-1]+pn[0:-2,1:-1])*dy**2 + \
        (pn[1:-1,2:]+pn[1:-1,0:-2])*dx**2)/(dx**2+dy**2) - \
        0.5*(dx*dy)**2/(dx**2+dy**2)*b[1:-1,1:-1]
            
        # Introducing the periodic boundary condition for x = 0
        p[0,1:-1] = 0.5*((pn[1,1:-1]+pn[-1,1:-1])*dy**2 + \
        (pn[0,2:]+pn[0,0:-2])*dx**2)/(dx**2+dy**2) - \
        0.5*(dx*dy)**2/(dx**2+dy**2)*b[0,1:-1]

        # Introducing the periodic boundary condition for x = 2
        p[-1,1:-1] = 0.5*((pn[0,1:-1]+pn[-2,1:-1])*dy**2 + \
        (pn[-1,2:]+pn[-1,0:-2])*dx**2)/(dx**2+dy**2) - \
        0.5*(dx*dy)**2/(dx**2+dy**2)*b[-1,1:-1]
        
        # Introducing the wall boundary conditions dp/dy = 0 at y =0,2            
        p[0,:] = p[1,:]       # dp/dy = 0 at y = 0
        p[-1,:] = p[-2,:]     # dp/dy = 0 at y = 2
    
    return p


# Defining the solver function:
def channelSolver(nt, u, v, dt, dx, dy, p, rho, nu, F):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    
    udiff = 1
    stepcount = 0
    
    while udiff > 0.001:
        un = u.copy()
        vn = v.copy()
        
        b = buildUpB(rho, dt, dx, dy, u, v)
        p = pressPoisson(p, dx, dy, b)
        
        # Calculating U component
        # The same scheme as in step 11
        u[1:-1,1:-1] = un[1:-1,1:-1] - \
        un[1:-1,1:-1]*(dt/dx)*(un[1:-1,1:-1]-un[0:-2,1:-1]) - \
        vn[1:-1,1:-1]*(dt/dy)*(un[1:-1,1:-1]-un[1:-1,0:-2]) - \
        (dt/(2*rho*dx))*(p[2:,1:-1]-p[0:-2,1:-1]) + \
        nu*(dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]) + \
        dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])) + F*dt
        
        # Introducing the periodic boundary condition for x = 0
        # [2:,] -> [1,]; [1:-1,] -> [0,]; [0:-2,] -> [-1,]
        u[0,1:-1] = un[0,1:-1] - \
        un[0,1:-1]*(dt/dx)*(un[0,1:-1]-un[-1,1:-1]) - \
        vn[0,1:-1]*(dt/dy)*(un[0,1:-1]-un[0,0:-2]) - \
        (dt/(2*rho*dx))*(p[1,1:-1]-p[-1,1:-1]) + \
        nu*(dt/dx**2*(un[1,1:-1]-2*un[0,1:-1]+un[-1,1:-1]) + \
        dt/dy**2*(un[0,2:]-2*un[0,1:-1]+un[0,0:-2])) + F*dt        
        
        # Introducing the periodic boundary condition for x = 2
        # [2:,] -> [0,]; [1:-1,] -> [-1,]; [0:-2,] -> [-2,]
        u[-1,1:-1] = un[-1,1:-1] - \
        un[-1,1:-1]*(dt/dx)*(un[-1,1:-1]-un[-2,1:-1]) - \
        vn[-1,1:-1]*(dt/dy)*(un[-1,1:-1]-un[-1,0:-2]) - \
        (dt/(2*rho*dx))*(p[0,1:-1]-p[-2,1:-1]) + \
        nu*(dt/dx**2*(un[0,1:-1]-2*un[-1,1:-1]+un[-2,1:-1]) + \
        dt/dy**2*(un[-1,2:]-2*un[-1,1:-1]+un[-1,0:-2])) + F*dt        
        
        
        # Calculating V component 
        # The same scheme as in step 11
        v[1:-1,1:-1] = vn[1:-1,1:-1] - \
        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[0:-2,1:-1]) - \
        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[1:-1,0:-2]) - \
        dt/(2*rho*dy)*(p[1:-1,2:]-p[1:-1,0:-2]) + \
        nu*(dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1]) + \
        dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2]))
        
        # Introducing the periodic boundary condition for x = 0
        # [2:,] -> [1,]; [1:-1,] -> [0,]; [0:-2,] -> [-1,]        
        v[0,1:-1] = vn[0,1:-1] - \
        un[0,1:-1]*dt/dx*(vn[0,1:-1]-vn[-1,1:-1]) - \
        vn[0,1:-1]*dt/dy*(vn[0,1:-1]-vn[0,0:-2]) - \
        dt/(2*rho*dy)*(p[0,2:]-p[0,0:-2]) + \
        nu*(dt/dx**2*(vn[1,1:-1]-2*vn[0,1:-1]+vn[-1,1:-1]) + \
        dt/dy**2*(vn[0,2:]-2*vn[0,1:-1]+vn[0,0:-2]))        
              
        # Introducing the periodic boundary condition for x = 2
        # [2:,] -> [0,]; [1:-1,] -> [-1,]; [0:-2,] -> [-2,]        
        v[-1,1:-1] = vn[-1,1:-1] - \
        un[-1,1:-1]*dt/dx*(vn[-1,1:-1]-vn[-2,1:-1]) - \
        vn[-1,1:-1]*dt/dy*(vn[-1,1:-1]-vn[-1,0:-2]) - \
        dt/(2*rho*dy)*(p[-1,2:]-p[-1,0:-2]) + \
        nu*(dt/dx**2*(vn[0,1:-1]-2*vn[-1,1:-1]+vn[-2,1:-1]) + \
        dt/dy**2*(vn[-1,2:]-2*vn[-1,1:-1]+vn[-1,0:-2]))        
        
        # Introducing the wall boundary conditions: u,v=0 @ y=0,2
        u[:,0] = 0
        u[:,-1] = 0
        v[:,0] = 0
        v[:,-1] = 0
        
        # Checking for convergence
        udiff = (np.sum(np.abs(u[:])-np.abs(un[:])))/np.sum(np.abs(u[:]))
        
        # Counting the steps
        stepcount += 1        
        
    return u, v, p, stepcount


# Defining the plotting function
def VectorPlot2D(u, v, Y, X):
    pl.figure(figsize = (11,7), dpi = 100)
    pl.quiver(X[::3,::3],Y[::3,::3],u[::3,::3],v[::3,::3]) # plotting velocity vectors
    pl.xlabel('X')
    pl.ylabel('Y')    
    pl.title('Channel Flow - Velocity profiles')
    #pl.quiver(X,Y,u,v)

# Defining the call function for the problem
def channelFlow(xmin, xmax, ymin, ymax, nt):
    # Number of nodes
    nx = 41
    ny = 41
    
    # Cell size
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)
    
    # Step size
    dt = dx/5.0
    
    # The grid
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)
    Y, X = np.meshgrid(y,x)
    
    # Fluid properties
    rho = 1
    nu = 0.1
    F = 1
    
    # Matrix initialization
    u = np.zeros((nx, ny))
    v = np.zeros((nx, ny))
    p = np.zeros((nx, ny))
    
    u, v, p, steps = channelSolver(nt, u, v, dt, dx, dy, p, rho, nu, F)
    VectorPlot2D(u, v, Y, X)
    #print steps

# Variable declaration
xmin = 0.0
xmax = 2.0
ymin = 0.0
ymax = 2.0
nt = 10

# Run the function
channelFlow(xmin, xmax, ymin, ymax, nt)

