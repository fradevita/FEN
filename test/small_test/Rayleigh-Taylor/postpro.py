import numpy as np
import matplotlib.pyplot as plt
import subprocess
import math
import os
import glob
import pandas
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

Nx = 128
Ny = 512
Nz = 1
Sf = Nx*Ny*Nz
Lx = 1.0
Ly = 4.0
dx = Lx/Nx
dy = Ly/Ny
assert(dx==dy)

X = np.linspace(dx/2, Lx-dx/2, Nx)
Y = np.linspace(dy/2, Ly-dy/2, Ny)

n = 0
data = np.fromfile('data/vof_'+str(n).zfill(1))
vof0 = np.reshape(data, (Nx,Ny,Nz), order = 'F')
n = 1
data = np.fromfile('data/vof_'+str(n).zfill(1))
vof1 = np.reshape(data, (Nx,Ny,Nz), order = 'F')
n = 2
data = np.fromfile('data/vof_'+str(n).zfill(1))
vof2 = np.reshape(data, (Nx,Ny,Nz), order = 'F')
n = 3
data = np.fromfile('data/vof_'+str(n).zfill(1))
vof3 = np.reshape(data, (Nx,Ny,Nz), order = 'F')
n = 4
data = np.fromfile('data/vof_'+str(n).zfill(1))
vof4 = np.reshape(data, (Nx,Ny,Nz), order = 'F')

fig, ax = plt.subplots(ncols = 5, nrows = 1)

ax[0].set_xlim([0.,1.])
ax[0].set_title(r'$t = 0$')
ax[0].axes.xaxis.set_visible(False)
ax[0].axes.yaxis.set_visible(False)
ax[0].contour(X, Y, np.transpose(vof0[:,:,0]), [0.5])
ax[1].set_xlim([0.,1.])
ax[1].set_title(r'$t = 0.7$')
ax[1].axes.xaxis.set_visible(False)
ax[1].axes.yaxis.set_visible(False)
ax[1].contour(X, Y, np.transpose(vof1[:,:,0]), [0.5])
ax[2].set_xlim([0.,1.])
ax[2].set_title(r'$t = 0.8$')
ax[2].axes.xaxis.set_visible(False)
ax[2].axes.yaxis.set_visible(False)
ax[2].contour(X, Y, np.transpose(vof2[:,:,0]), [0.5])
ax[3].set_xlim([0.,1.])
ax[3].set_title(r'$t = 0.9$')
ax[3].axes.xaxis.set_visible(False)
ax[3].axes.yaxis.set_visible(False)
ax[3].contour(X, Y, np.transpose(vof3[:,:,0]), [0.5])
ax[4].set_xlim([0.,1.])
ax[4].set_title(r'$t = 1.0$')
ax[4].axes.xaxis.set_visible(False)
ax[4].axes.yaxis.set_visible(False)
ax[4].contour(X, Y, np.transpose(vof4[:,:,0]), [0.5])

plt.tight_layout()
plt.savefig('plot.png')
if (display): plt.show()
