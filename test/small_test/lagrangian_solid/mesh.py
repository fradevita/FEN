import numpy as np

L = 1
Lx = 4.0*L
Nx = 128

# Eulerian mesh
delta = Lx/Nx

# lagrangian resolution
rl = delta*0.7

# Number of points
nv = int(L/rl) 

# Actual resolution
rl = L/(nv - 1)

mesh = np.zeros((nv,2))
for n in range(nv):
    mesh[n,0] = 0.0 + n*rl

# Salvo la griglia
np.savetxt('mesh.txt', mesh, delimiter = ' ')
