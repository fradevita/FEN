import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import math
import json
import os 
import glob
import subprocess
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

# READ SETUT FILE
setup_file = open('setup.json')
setup = json.load(setup_file)

# Setup Grid
Nx = setup["Grid"]["Nx"]
Ny = setup["Grid"]["Ny"]
Nz = setup["Grid"]["Nz"]
Lx = setup["Grid"]["Lx"]
Ly = setup["Grid"]["Ly"]
Lz = setup["Grid"]["Lz"]
x0 = setup["Grid"]["origin"][0]
y0 = setup["Grid"]["origin"][1]
dx = Lx/Nx
dy = Ly/Ny
assert(dx == dy)
X = np.linspace(dx/2, Lx - dx/2, Nx)
Y = np.linspace(dy/2, Ly - dy/2, Ny)

# Read data
data = np.fromfile('u.dat')
u = np.reshape(data,(Nx,Ny,Nz), order = 'F')
data = np.fromfile('v.dat')
v = np.reshape(data,(Nx,Ny,Nz), order = 'F')
data = np.fromfile('vort.dat')
vort = np.reshape(data,(Nx,Ny,Nz), order = 'F')

# Extract profiles
uy = np.zeros(np.size(Y))
for j in range(0, np.size(Y)):
    uy[j] = 0.5*(u[int(Nx/2),j,0] + u[int(Nx/2-1),j,0])
vx = np.zeros(np.size(X))
for i in range(0, np.size(X)):
    vx[i] = 0.5*(v[i,int(Ny/2),0] + v[i,int(Ny/2-1),0])

# Load reference solution
uref = np.genfromtxt('uref')
vref = np.genfromtxt('vref')

plt.figure()
plt.xlabel('y')
plt.ylabel('u')
plt.title('Lid driven cavity test case')
plt.xlim([0, 1])
plt.plot(Y + dy/2, uy)
plt.plot(uref[:,0] + 0.5, uref[:,1], 'o')
plt.savefig('u_y.png')
if (display): plt.show()
plt.close()

plt.figure()
plt.xlabel('x')
plt.ylabel('v')
plt.title('Lid driven cavity test case')
plt.xlim([0, 1])
plt.plot(X + dx/2, vx)
plt.plot(vref[:,0] + 0.5, vref[:,1], 'o')
plt.savefig('v_x.png')
if (display): plt.show()
plt.close()

fig, ax = plt.subplots()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Lid driven cavity test case')
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.contourf(X, Y, np.transpose(vort[:,:,0]), 128, vmin = -10, vmax = 10, cmap = 'jet')
ax.axis('equal')
plt.savefig('vort.png')
if (display): plt.show()
plt.close()
