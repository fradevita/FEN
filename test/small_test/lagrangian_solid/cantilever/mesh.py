import numpy as np

L = 1        # Cantilever length
N = 64       # Number of points
rl = L/(N-1) # edge size

# Create the mesh
mesh = np.zeros((N,2))
for n in range(N):
    mesh[n,0] = 0.0 + n*rl

# Salvo la griglia
np.savetxt('mesh.txt', mesh, delimiter = ' ')
