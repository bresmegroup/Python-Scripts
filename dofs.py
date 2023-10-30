
#%% 
# The Bresme Group 2023
# Calculate degrees of freedom of each in a water molecule modelled with rigid bonds. 
# The degrees of freedom can be used to calculate the atomic temperatures.
# Reference: Olarte-Plata and Bresme, J. Chem. Phys. 2022, 156, 204701  doi: 10.1063/5.0090983
#
from math import degrees
import numpy as np

# SPC/E
# r0 = 1.0 # Angstrom
# theta0 = 109.47 # degrees

# TIP3P, TIP4P, TIP4P/2005
r0 = 0.9572 # Angstrom
theta0 = 104.52 # degrees

theta0_rad = theta0 * np.pi / 180. # radians

# Atomic weight
mh = 1.0078 # g/mol
mo = 15.999 # g/mol

# H coordinates
xh =  r0*np.sin(theta0_rad/2)
yh = -r0*np.cos(theta0_rad/2)
print ("xh=", xh, " yh=", yh)

# Center of mass
xcom = 0
ycom = 2*mh*yh/(mo + 2*mh)
print ("xcom=", xcom, " ycom=", ycom)

# Moment of inertia
Ix = 2*mh*xh**2
Iy = 2*mh*(yh-ycom)**2 + mo*ycom**2
Iz = 2*mh*((xh**2 + (yh - ycom)**2)) + mo*ycom**2

print ("Ix=", Ix)
print ("Iy=", Iy)
print ("Iz=", Iz)

# Angular velocity
eps=1/2
wx= np.sqrt(2*eps/Ix)
wy= np.sqrt(2*eps/Iy)
wz= np.sqrt(2*eps/Iz)

print ("wx=", wx)
print ("wy=", wy)
print ("wz=", wz)

# Linear velocities
vx = vy = vz = np.sqrt(2*eps/(mo+2*mh))
print ("vx=", vx)
print ("vy=", vy)
print ("vz=", vz)

# Kinetic energy hydrogen
Ek_h = mh/2*(vx**2 + vy**2 + vz**2) + mh/2*((wx*xh)**2 + (wy*(yh-ycom))**2 + (wz*np.sqrt(xh**2 + (yh-ycom)**2))**2)

print ("Ek_h=", Ek_h) 

# Kinetic energy oxygen
Ek_o = 1/2*mo*(vx**2 + vy**2 + vz**2) + 1/2*mo*((wy*ycom)**2 + (wz*ycom)**2)
print ("Ek_o=", Ek_o)

# Number of degrees of freedom, oxygen and hydrogen
Ndof_o = 2*Ek_o
Ndof_h = 2*Ek_h

print ("====================")
print ("Ndof_o=", Ndof_o)
print ("Ndof_h=", Ndof_h)
print ("====================")


# %%
