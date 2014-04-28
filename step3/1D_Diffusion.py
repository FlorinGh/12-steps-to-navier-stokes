# Step3: Diffusion only
# in this step only the diffusion term from the NS eq is used to determine the velocity
# the diffusion term is a second order derivative multiplied by the diffusive coefficient
# to develop the liniar algebraic aproximation of the second order derivative,
# the Taylor expansion is used

import numpy as np
import pylab as pl
pl.ion()

nx = 41
dx = 2./(nx-1)

nt = 20
nu = 0.3

sigma = 0.2
dt = sigma*dx**2/nu

u = np.ones(nx)
u[.5/dx : 1/dx+1] = 2

un = np.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1,nx-1):
        u[i] = un[i] + nu*dt/dx**2*(un[i+1] - 2*un[i] + un[i-1])
    pl.plot(np.linspace(0,2,nx), u)

