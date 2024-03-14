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
data_064 = pandas.read_csv('data/C_064')
data_128 = pandas.read_csv('data/C_128')
data_256 = pandas.read_csv('data/C_256')

# Flow data
Umean = 0.2
L = 0.1

# Reference values
Cd = 5.57953523384
Cl = 0.010618948146

# Array for plotting reference values
Cdref = np.ones(len(data_256["t"]))

fig, ax = plt.subplots()
ax.set_xlim([0, 3])
ax.set_xlabel(r'$t$')
ax.set_ylim([3, 6])
ax.set_ylabel(r'$C_d$')
ax.plot(data_064["t"][::10], data_064["Fd"][::10]*2/Umean**2/L, color = 'red')
ax.plot(data_064["t"][::10], data_064["Fd"][::10]*2/Umean**2/L, 'o', color = 'red', label = 'N = 64, probes')

ax.plot(data_128["t"][::50], data_128["Fd"][::50]*2/Umean**2/L, color = 'blue')
ax.plot(data_128["t"][::50], data_128["Fd"][::50]*2/Umean**2/L, '^', color = 'blue', label = 'N = 128, probes')

ax.plot(data_256["t"][::50], data_256["Fd"][::50]*2/Umean**2/L, 's', color = 'green', label = 'N = 256, probes')
ax.plot(data_256["t"][::50], data_256["Fd"][::50]*2/Umean**2/L, color = 'green')

ax.plot(data_256["t"], Cdref*Cd, '--', color = 'black', label = r'reference $C_d$')

plt.tight_layout()
plt.legend()
plt.savefig('drag.png')
if (display): plt.show()
plt.close()

# Compute error
e_P = np.zeros(3)
e_F = np.zeros(3)

e_P[0] = abs(data_064["Fd"].to_numpy()[-1]*2/Umean**2/L - Cd)/Cd
e_P[1] = abs(data_128["Fd"].to_numpy()[-1]*2/Umean**2/L - Cd)/Cd
e_P[2] = abs(data_256["Fd"].to_numpy()[-1]*2/Umean**2/L - Cd)/Cd

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
plt.loglog(N, e_P, '-o', label = 'probes')
plt.loglog(N, scaling1, '--', color = 'black', label = '1/N')
plt.loglog(N, scaling2, '-.', color = 'black', label = r'$1/N^{2}$')
plt.legend()
plt.savefig('error.png')
if (display): plt.show()
plt.close()
