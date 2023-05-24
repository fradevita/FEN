import numpy as np
import math

L = 1
Nx = 512
Lx = 8.0
theta = 0.1*math.pi
x0 = Lx/2
y0 = 2.0
dx = Lx/Nx
deltal = 0.7*dx
Nv = int(L/deltal)
deltal = L/(Nv-1)

mesh = np.zeros((Nv,2))
for i in range(Nv):
    mesh[i,0] = 4.0 + i*deltal*math.sin(theta)
    mesh[i,1] = 2.0 + i*deltal*math.cos(theta)
    
# Salvo la griglia
np.savetxt('mesh.txt', mesh, delimiter = ' ')
