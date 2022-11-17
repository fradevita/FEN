import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import math
import os 
import glob
import subprocess
from scipy.optimize import curve_fit
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

data = np.genfromtxt('error')
N = data[:,0]
Linf_u = data[:,1]
Linf_v = data[:,2]
Linf_w = data[:,3]
L2_u = data[:,4]
L2_v = data[:,5]
L2_w = data[:,6]

# Fit the error
def objective(x, a, b):
    return a*x + b

pars, cov = curve_fit(f=objective, xdata=np.log(N), ydata=np.log(Linf_u))
a, b = pars
pars, cov = curve_fit(f=objective, xdata=np.log(N), ydata=np.log(Linf_v))
c, d = pars
pars, cov = curve_fit(f=objective, xdata=np.log(N), ydata=np.log(Linf_w))
e, f = pars

if (-a > 1.5 and -c > 1.5 and -e > 1.5):
  print('Test completed.')
else:
  print('The convergence rate of the Poisson solver is less than second order')

# Plot
scaling1 = np.zeros(len(N))
scaling2 = np.zeros(len(N))
for i in range(0,np.size(N)):
  scaling1[i] = 1.0/N[i]
  scaling2[i] = 1.0/N[i]**2

plt.plot(N, Linf_u, 'o', label = r'$L_{\inf}(u)$')
plt.plot(N, Linf_v, '^', label = r'$L_{\inf}(v)$')
plt.plot(N, Linf_w, 'v', label = r'$L_{\inf}(w)$')
plt.plot(N, L2_u, 'o', label = r'$L_{2}(u)$')
plt.plot(N, L2_v, '^', label = r'$L_{2}(v)$')
plt.plot(N, L2_w, 'v', label = r'$L_{2}(w)$')
plt.plot(N, scaling1,  '-', color = 'black', label = r'$1/N$')
plt.plot(N, scaling2, '--', color = 'black', label = r'$1/N^2$')
plt.xlabel("N")
plt.ylabel(r"$L_{\inf},L_{2}$")
plt.title("3D ABC flow test case")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
