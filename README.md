# 12 Steps to Navier-Stokes Equations

## **Challenge**

Use computing capabilities of Python to solve the nonlinear coupled partial derivative equations that govern the dynamics of fluids, the Navier-Stokes equations:

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/navier-stokes-equations.jpg)

## **Actions**

* creating implicit numerical schemes to solve ever increasing difficult components of the NS equations
* linear convection:

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp5-linear-convection-2D/2D_linear_conv_initial_conditions.png)


![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp5-linear-convection-2D/2D_linear_conv_solution.png)

* nonlinear convection:

```python
# Step6: 2D Nonlinear Convection
# the ecuations to be solved are: du/dt + u*du/dx +v*du/dy = 0
#                                 dv/dt + u*dv/dx +v*dv/dy = 0
# in this equation we have the acumulation and convection terms
# but this time the convection is nonlinear and the equations are coupled
# this time the problem is 2D and we want to solve for the velocities

# First we use the equations given in the IPython notebook
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
pl.ion()

# Variable declaration
nx = 501 # initial value 101
ny = 501 # initial value 101
nt = 400 # initial value 80

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

from pylab import cm # this will be used to create a colormap 3d plot
fig = pl.figure(figsize=(11,7), dpi=100)
ax = Axes3D(fig)
X,Y = np.meshgrid(x,y)
ax.plot_surface(X,Y,v,cmap=cm.coolwarm) # plot the u component of the velocity
pl.xlabel('X')
pl.ylabel('Y')
pl.title('2D Nonlinear Convection: Solution')
# to plot the v component just change u with v
# because the equations are symetrical, the plots are identical
```

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp6-nonlinear-convection-2D/nonlin_conv_2D_solution_2.png)

* diffusion:

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp7-diffusion-2D/diffusion_initial_conditions.png)

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp7-diffusion-2D/sol_diffusion_10.png)

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp7-diffusion-2D/sol_diffusion_30.png)

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/blob/master/stp7-diffusion-2D/sol_diffusion_270.png)

* Burgers' equation

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/tree/master/stp8-burgers-equation-2D/burgers_ic.png)

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/tree/master/stp8-burgers-equation-2D/burgers_sol_120.png)

![](.https://github.com/FlorinGh/12-steps-to-navier-stokes/tree/master/stp8-burgers-equation-2D/burgers_sol_1200.png)

* cavity flow

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/tree/master/stp11-cavity-problem/cav_sol_10.png)

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/tree/master/stp11-cavity-problem/solution_5000.png)

* channel flow

![](https://github.com/FlorinGh/12-steps-to-navier-stokes/tree/master/stp12-channel-flow-problem/solution.png)

## **Results**

The result of this exercise was package of numerical solutions to the difficult equations of fluid dynamics; the implementation is only in 2D and can solve any problem that can be formulated in a structured 2D mesh; the main equations take also into account turbulence and as seen in the results of teh cavity problem, turbulence is modelled implicitly in the solutions of this project.

For a complete overview of this project please visit its dedicated repository on github:     [https://github.com/FlorinGh/12-steps-to-navier-stokes](https://github.com/FlorinGh/12-steps-to-navier-stokes)​.

