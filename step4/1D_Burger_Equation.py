# Step4: Burger's equation
# this is a combination of the last two equations solved
# it's a 1D equation of NS
# it contains the accumulation, convection and diffusion terms

import numpy as np
import sympy as sp
import pylab as pl
pl.ion()

# Setting up the symbolic variables
x, nu, t = sp.symbols('x nu t')
phi = sp.exp(-(x-4*t)**2/(4*nu*(t+1))) + sp.exp(-(x-4*t-2*np.pi)**2/(4*nu*(t+1)))

# Evaluate the partial derivative dphi/dx using sympy
phiprime = phi.diff(x)
#print phiprime

# Create the initial conditions function
u = -2*nu*(phiprime/phi)+4
#print u

# Transform the symbolic equation into a function using lambdify
from sympy.utilities.lambdify import lambdify
ufunc = lambdify ((t, x, nu), u)
#print ufunc(1,4,3)

# Variables declaration
nx = 101
nt = 100
dx = 2*np.pi/(nx-1)
nu = 0.07
dt = dx*nu
T = nt*dt

grid = np.linspace(0, 2*np.pi, nx)
un = np.empty(nx)
t = 0

# Initializing the velocity function
u = np.asarray([ufunc(t, x, nu) for x in grid])
#print u

# Make a plot of the initial conditions
pl.figure(figsize=(11,7), dpi=100)
pl.plot(grid,u, marker='o', lw=2)
pl.xlim([0,2*np.pi])
pl.ylim([0,10])

# Apply the scheme with the periodic boundary conditions in mind
for n in range(nt):
    un = u.copy()
    for i in range(nx-1):
        u[i] = un[i] - un[i] * dt/dx * (un[i]-un[i-1]) + \
            nu * dt/(dx**2) * (un[i+1] - 2*un[i] + un[i-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1]-un[-2]) + \
            nu * dt/(dx**2) * (un[0] - 2*un[-1] + un[-2])

u_analytical = np.asarray([ufunc(T, xi, nu) for xi in grid])

pl.figure(figsize=(11,7), dpi=100)
pl.plot(grid, u, marker='o', lw=2, label='Computational')
pl.plot(grid, u_analytical, label='Analytical')
pl.xlim([0, 2*np.pi])
pl.ylim([0,10])
pl.legend()