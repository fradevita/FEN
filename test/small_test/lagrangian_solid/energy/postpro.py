import numpy as np
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

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
ax[1].plot(t, abs(W - (Ek + Ep)))

plt.tight_layout()
if (display):
    plt.show()
    print('        Test completd.')
else:
    plt.savefig("energy.png")
    print('        Test completd, check image.')
plt.close()
