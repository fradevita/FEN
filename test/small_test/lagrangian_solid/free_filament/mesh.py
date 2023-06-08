import numpy as np
import math

L = 1
N = 100
theta = 0.1*math.pi

# Eulerian mesh
delta = L/N

mesh = np.zeros((N,2))
for i in range(N):
    mesh[i,0] = 0.0 + i*delta*math.sin(theta)
    mesh[i,1] = 0.0 - i*delta*math.cos(theta)
    
# Salvo la griglia
np.savetxt('mesh.txt', mesh, delimiter = ' ')
