import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys
import json
import glob

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

###############################################################################
# READ SETUP FILE
###############################################################################
grid = json.load(open('unset.json'))

# Setup Grid
Nx = grid["Grid"]["Nx"]
Ny = grid["Grid"]["Ny"]
Nz = grid["Grid"]["Nz"]
Lx = grid["Grid"]["Lx"]
Ly = grid["Grid"]["Ly"]
Lz = grid["Grid"]["Lz"]
x0 = grid["Grid"]["origin"][0]
y0 = grid["Grid"]["origin"][1]
dx = Lx/Nx
dy = Ly/Ny
assert(dx == dy)
X = np.linspace(dx/2, Lx - dx/2, Nx)
Y = np.linspace(dy/2, Ly - dy/2, Ny)

###############################################################################
# Load data
############################################################################### 
filelist = glob.glob('data/vof_*')
f0 = np.reshape(np.fromfile(sorted(filelist)[0]), (Nx,Ny), order = 'F' ) 
f1 = np.reshape(np.fromfile(sorted(filelist)[int(len(filelist)/2)]), (Nx,Ny), order = 'F' ) 
f2 = np.reshape(np.fromfile(sorted(filelist)[-1]), (Nx,Ny), order = 'F' ) 

# Plot
fig, ax = plt.subplots()
X, Y = np.meshgrid(X, Y)
plt.contour(X, Y, np.transpose(f0), [0.5], colors = 'black', linestyles = 'solid')
plt.contour(X, Y, np.transpose(f1), [0.5], colors = 'blue', linestyles = 'solid')
plt.contour(X, Y, np.transpose(f2), [0.5], colors = 'red', linestyles = 'dashed')
ax.set_aspect('equal', adjustable='box')
plt.savefig("fig.png")
if (display): plt.show()
plt.close(fig)
