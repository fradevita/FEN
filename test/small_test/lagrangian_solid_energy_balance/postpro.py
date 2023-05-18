import numpy as np
import matplotlib.pyplot as plt
import sys

# Simulation data
data  = np.genfromtxt('energy.dat')

fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize = (18,12))

t = data[:,0]
W = data[:,1]
Ep = data[:,2]
Ek = data[:,3]

ax[0].set_xlabel('t')
ax[0].plot(t, W, label = 'Work')
ax[0].plot(t, Ep, label = 'Ep')
ax[0].plot(t, Ek, label = 'Ek')
ax[0].legend()
ax[1].set_xlabel('t')
ax[1].set_ylabel('|W - Em|')
ax[1].plot(t, W - (Ek + Ep))

plt.tight_layout()
plt.show()
plt.close()
