# Step2: Nonlinear Convection
# in this step the convection term of the NS equation is solved in 1D
# this time the wave velocity is nonlinear and this is the case in NS equation

import numpy as np
import pylab as pl
pl.ion()

nx = 41
dx = 2./(nx-1)

nt = 20
dt = 0.025

u = np.ones(nx)
u[0.5/dx:1/dx+1] = 2.0

un = np.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1,nx):
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1])

pl.plot(np.linspace(0,2,nx),u)