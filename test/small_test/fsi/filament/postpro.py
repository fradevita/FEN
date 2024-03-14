import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

Lx = 8
Ly = 8
Nx = 512
Ny = 512 
dx = Lx/Nx
dy = Ly/Ny
assert(dx==dy)
X = np.linspace(dx/2, Lx - dx/2, Nx) 
Y = np.linspace(dy/2, Ly - dy/2, Ny)
dt = 0.100000E-02
iout = 100

####################################################################################################
# Trailing edge plot
####################################################################################################
data = np.genfromtxt('out.txt')
Huang = np.genfromtxt('Huang.txt')
plt.figure()
plt.xlim([0, 30])
plt.ylim([-0.6, 0.8])
plt.xlabel(r'$tU/L$')
plt.ylabel(r'$y_{t.e.}/L$')
plt.plot(data[:,0], data[:,1] - Lx/2, label = 'data')
plt.plot(Huang[:,0], Huang[:,1], 'o', label = 'Huang et al 2007')
plt.legend(loc = 'upper center')
if (display): 
    plt.show()
else:
    plt.savefig('TE.png')
plt.close()

####################################################################################################
# Filament plot
####################################################################################################
plt.figure()
nin = 20000
nen = 23100
for n in range(nin,nen+iout,iout):
    XY = np.genfromtxt('data/sb_'+str(n).zfill(7))
    plt.plot(XY[:,1], XY[:,0], color = 'black')
#plt.plot(data[nin:nen,2], data[nin:nen,1])
if (display): plt.show()
plt.close()

####################################################################################################
# Forces plot
####################################################################################################
data = np.genfromtxt('forces.txt')
Cd = np.genfromtxt('Cd.txt')
Cl = np.genfromtxt('Cl.txt')

Cd_1 = np.genfromtxt('Cd_old.txt')
Cl_1 = np.genfromtxt('Cl_old.txt')
fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize = (15,10))

ax[0].set_xlim([20, 26])
ax[0].set_ylim([0., 2.])
ax[0].plot(data[:,0], 2.*data[:,2], label = 'data')
ax[0].plot(  Cd[:,0],      Cd[:,1], label = 'de Tullio & Pascazio 2016')
ax[0].plot(  Cd_1[:,0],      Cd_1[:,1], label = 'ref')
ax[0].legend(loc = 'upper center')

ax[1].set_xlim([20, 26])
ax[1].set_ylim([-2., 3.])
ax[1].plot(data[:,0], 2.*data[:,1], label = 'data')
ax[1].plot(  Cl[:,0],      Cl[:,1], label = 'de Tullio & Pascazio 2016')

ax[1].plot(  Cl_1[:,0],      Cl_1[:,1], label = 'ref')
ax[1].legend(loc = 'upper center')

plt.tight_layout()
if (display): 
    plt.show()
else:
    plt.savefig('forces.png')
plt.close()

####################################################################################################
# Vorticity plot
####################################################################################################
c = 0
Nt = [20000, 20800, 21600, 22400]
for n in Nt:
    XY = np.genfromtxt('data/sb_'+str(n).zfill(7))
    data = np.fromfile('data/vrt_'+str(n).zfill(7))
    vrt = np.reshape(data, (Nx,Ny), order = 'F')[:,:]
    fig, ax = plt.subplots()
    ax.set_xlim([1.5, 6])
    ax.set_ylim([3, 5])
    ax.plot(XY[:,1], XY[:,0], color = 'red')
    CS = ax.contour(X, Y, vrt, np.linspace(-40,40,30), colors = 'k')
    plt.savefig('data/img_'+str(c).zfill(1))
    plt.close()
    c = c + 1
