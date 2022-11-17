import numpy as np
import matplotlib.pyplot as plt
import math
import os 
import glob
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

kc = 1.0
n = 2.0
dpdx = 1.0
h = 1
Nx = 4
Ny = 32
Nz = 1
Ly = 1
Sf = Nx*Ny*Nz
dy = 1.0/Ny
y = np.linspace(dy/2, Ly-dy/2, Ny)
# Read data from file
data = np.fromfile('u_032')
u_raw = data[0:Sf]
u = np.transpose(np.reshape(u_raw,(Nx,Ny,Nz), order = 'F')[:,:,0])
u_sim = u[:,2]

# Soluzione analitica
ua = np.zeros(100)
ya = np.linspace(-0.5, 0.5, 100)
u0 = ((1/kc)*abs(dpdx))**(1/n)*(h/2)**(1+1/n)*(n/(n+1))
for j in range(100):
    ua[j] = u0*(1 - (2*abs(ya[j])/h)**(1 + 1/n))

plt.figure()
plt.plot(u_sim, y, 'o')
plt.plot(ua, ya+0.5)
plt.show()
plt.close()
