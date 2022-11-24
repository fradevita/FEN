import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys
from scipy.optimize import curve_fit

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

# Reference values
SIM = np.genfromtxt('fig1a.SIM')
KL = np.genfromtxt('fig1a.KL')
KLf = np.genfromtxt('fig19.f')
KLp = np.genfromtxt('fig19.p')
DATA = np.genfromtxt('out.txt')

plt.figure()
plt.xlim([0.01, 2.5])
plt.ylim([0, 2])
plt.xlabel(r'$tU/D$')
plt.ylabel(r'$C_d$')
plt.plot(SIM[:,0]/2, SIM[:,1], label = 'SIM (Mohahegh et al 2017)')
plt.plot(KL[:,0]/2, KL[:,1], 'o', label = 'K and L. 1995')
plt.plot(KLf[:,0]/2, KLf[:,1], 'x', label = 'friction, K and L. 1995')
plt.plot(KLp[:,0]/2, KLp[:,1], 'v', label = 'pressure, K and L. 1995')
plt.plot(DATA[:,0], -DATA[:,1], label = 'FEN', color = 'black')
plt.legend()
plt.show()
plt.close()
