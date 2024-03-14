import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

####################################################################################################
# SETUP
####################################################################################################
grid_file = open('grid.json')
grid = json.load(grid_file)

# Setup Grid
Nx = grid["Grid"]["Nx"]
Ny = grid["Grid"]["Ny"]
Nz = grid["Grid"]["Nz"]
Lx = grid["Grid"]["Lx"]
Ly = grid["Grid"]["Ly"]
Lz = grid["Grid"]["Lz"]
x0 = grid["Grid"]["origin"][0]
y0 = grid["Grid"]["origin"][1]
z0 = grid["Grid"]["origin"][2]
dx = Lx/Nx
dy = Ly/Ny
dz = Lz/Nz
assert(dx == dy)
assert(dx == dz)
x = np.linspace(dx/2, Lx - dx/2, Nx)
y = np.linspace(dy/2, Ly - dy/2, Ny)
z = np.linspace(dz/2, Lz - dz/2, Nz)

####################################################################################################
# Load data
####################################################################################################
# Load simulation data
u = np.fromfile('u.raw')
u = np.reshape(u, (Nx,Ny,Nz), order = 'F')
uc = 0.5*(u[int(Nx/2),:,int(Nz/2)] + u[int(Nx/2-1),:,int(Nz/2)])

v = np.fromfile('v.raw')
v = np.reshape(v, (Nx,Ny,Nz), order = 'F')
vc = 0.5*(v[:,int(Ny/2),int(Nz/2)] + v[:,int(Ny/2-1),int(Nz/2)])

# Load reference data
Uref = pd.read_csv('Uref.csv')
Vref = pd.read_csv('Vref.csv')

####################################################################################################
# Plot
####################################################################################################
fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize = (15,8))
# u(y)
ax[0].set_xlabel(r'$u/U$')
ax[0].set_xlim([-0.6, 1.0])
ax[0].set_ylabel(r'$y/L$')
ax[0].set_ylim([ 0.0, 1.0])
ax[0].plot(Uref["u"], Uref["y"], 'sk', fillstyle = 'none', label = 'Ku et al (1987)')
ax[0].plot(       uc,         y,  '-'                    , label = 'data')
ax[0].grid()
ax[0].legend()
# v(x)
ax[1].set_xlabel(r'$x/L$')
ax[1].set_xlim([0, 1])
ax[1].set_ylabel(r'$v/U$')
ax[1].set_ylim([-0.6, 0.6])
ax[1].plot(Vref["x"], Vref["v"], 'sk', fillstyle = 'none', label = 'Ku et al (1987)')
ax[1].plot(        x,        vc,  '-'                    , label = 'data')
ax[1].grid()
ax[1].legend()
plt.tight_layout()
if (display):
    plt.show()
    print('        Test completed.')
else:
    plt.savefig('plot.png')
    print('        Test completed, check image.')
plt.close()
