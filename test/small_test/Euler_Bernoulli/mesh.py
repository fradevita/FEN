import numpy as np

L = 1
Nv = 256

# Resolution
rl = L/(Nv - 1)

mesh = np.zeros((Nv,2))
for n in range(Nv):
    mesh[n,0] = 0.0 + n*rl

# Salvo la griglia
np.savetxt('mesh.txt', mesh, delimiter = ' ')
