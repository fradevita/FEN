import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

Nx = 100 
Ny = 100
Nz = 1
Sf = Nx*Ny*Nz

Lx = math.pi 
Ly = math.pi
dx = Lx/Nx
dy = Ly/Ny
X = np.arange(0, Lx, dx)
Y = np.arange(0, Ly, dy)

# Read data from file 
data = np.fromfile('data/vof_0002000')
data_raw = data[0:Sf]
fmid = np.reshape(data_raw,(Nx,Ny,Nz))
data = np.fromfile('data/vof_0004000')
data_raw = data[0:Sf]
fend = np.reshape(data_raw,(Nx,Ny,Nz))

# Plot
fig, ax = plt.subplots()
X, Y = np.meshgrid(X, Y)
plt.contour(X+0.5*dx, Y+0.5*dy, fmid[:,:,0],[0.05,0.5,0.95])
plt.contour(X+0.5*dx, Y+0.5*dy, fend[:,:,0],[0.05,0.5,0.95])
ax.set_aspect('equal', adjustable='box')
plt.savefig("fig.png")
if (display): plt.show()
plt.close(fig)
