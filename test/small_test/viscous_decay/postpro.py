import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import math
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1


l = 1.0
h = l
k = 2*math.pi/l
g = 9.80665
omega = math.sqrt(g*k)
a = 0.005
rho = 1000
mu = rho*math.sqrt(g*l)*l/1.0e+4
nu = mu/rho
P = 2*math.pi/omega
Ep0 = 4903.0289577702924

Nx = 128
Ny = 256
Nz = 1
Sf = Nx*Ny*Nz
Lx = 1
Ly = 2
dx = Lx/Nx
dy = Ly/Ny
X = np.linspace(dx/2, Lx-dx/2, Nx)
Y = np.linspace(dy/2, Ly-dy/2, Ny)

energy = np.genfromtxt('energy.dat')
t = energy[:,0]
et = energy[:,1] + energy[:,2] + Ep0
gamma = 2*nu*k**2
decay = np.zeros(len(et))
for i in range(len(et)):
    decay[i] = math.exp(-2*gamma*t[i])

plt.plot(t/P, et/et[0], label = 'data')
plt.plot(t/P, decay, label = 'theory')
plt.xlim([0, 5])
plt.title('Viscous energy decay water wave')
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
plt.close()
