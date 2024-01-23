import numpy as np
import matplotlib.pyplot as plt
import pandas
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

# Load data
data_064 = np.genfromtxt('out_064.txt')
data_128 = np.genfromtxt('out_128.txt')
data_256 = np.genfromtxt('out_256.txt')

# Reference values
Cd = 5.57953523384
Cl = 0.010618948146

# Array for plotting reference values
Cdref = np.ones(len(data_128[:,0]))

fig, ax = plt.subplots()
ax.set_title('Drag')
ax.set_xlim([0, 3])
ax.set_xlabel(r'$t$')
ax.set_ylim([4, 6])
ax.set_ylabel(r'$C_d$')
ax.plot(data_064[:,0], -data_064[:,1], color = 'blue', label = 'Ny = 64')
ax.plot(data_064[:,0],  data_064[:,3], 'o', color = 'blue', label = 'Ny = 64, probes')
ax.plot(data_128[:,0], -data_128[:,1], color = 'green', label = 'Ny = 128')
ax.plot(data_128[:,0],  data_128[:,3], 'o', color = 'green', label = 'Ny = 128, probes')
ax.plot(data_256[:,0], -data_256[:,1], color = 'red', label = 'Ny = 256')
ax.plot(data_256[:,0],  data_256[:,3], 'o', color = 'red', label = 'Ny =256, probes')
ax.plot(data_128[:,0], Cdref*Cd, '--', color = 'black', label = r'reference $C_d$')
plt.tight_layout()
plt.legend()
plt.savefig('drag.png')
if (display): plt.show()
plt.close()

# Compute error
eCd_eul = np.zeros(3)
eCd_eul[0] = abs(-data_064[-1,1] - Cd)/Cd
eCd_eul[1] = abs(-data_128[-1,1] - Cd)/Cd
eCd_eul[2] = abs(-data_256[-1,1] - Cd)/Cd

eCd_prb = np.zeros(3)
eCd_prb[0] = abs(data_064[-1,3] - Cd)/Cd
eCd_prb[1] = abs(data_128[-1,3] - Cd)/Cd
eCd_prb[2] = abs(data_256[-1,3] - Cd)/Cd

# Resolution
N = [64, 128, 256]

# Plot
scaling1 = np.zeros(len(N))
scaling2 = np.zeros(len(N))
for i in range(0,np.size(N)):
    scaling1[i] = 2/N[i]
    scaling2[i] = 100/N[i]**2

plt.figure()
plt.xlabel(r'$N$')
plt.ylabel(r'$|C_d - C_d^{ref}|/C_d^{ref}$')
plt.loglog(N, eCd_eul, '-o', label = r'$|e|, eulerian$')
plt.loglog(N, eCd_prb, '-s', label = r'$|e|, probes$')
plt.loglog(N, scaling1, '--', color = 'black', label = '1/N')
plt.loglog(N, scaling2, '-.', color = 'black', label = r'$1/N^{2}$')
plt.legend()
plt.savefig('error.png')
if (display):
    plt.show()
    print('        Test completd.')
else:
    plt.savefig("cantilever.png")
    print('        Test completd, check image.')
plt.close()
