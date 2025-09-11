import numpy as np
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

Nz = 1
Lx = 2*math.pi
Ly = 2*math.pi

data = np.genfromtxt('error.txt')
N = data[:,0]
Linfx = data[:,1]
Linfy = data[:,2]

# Fit the error
def objective(x, a, b):
  return a*x + b

popt, _ = curve_fit(objective, np.log(N), np.log(Linfx))
a, b = popt

popt, _ = curve_fit(objective, np.log(N), np.log(Linfy))
c, d = popt

if (-a > 1.8 and -c > 1.8):
  print('        Test completed.')
else:
  print('        The convergence rate of the advective terms is less than second order')

# Plot
scalingx = np.zeros(len(N))
scalingy = np.zeros(len(N))
for i in range(0,np.size(N)):
  scalingx[i] = np.exp(b)*N[i]**a
  scalingy[i] = np.exp(d)*N[i]**c

plt.plot(N,    Linfx,     'o', label = 'data')
plt.plot(N, scalingx, 'black', label = '$N^{%2.2f}$' % a)
plt.plot(N,    Linfy,     'x', label = 'data')
plt.plot(N, scalingy, 'black', label = '$N^{%2.2f}$' % c)
plt.xlabel("N")
plt.ylabel("error")
plt.title("advective term test case")
plt.xscale("log")
plt.yscale("log")
#plt.xticks(N, N)
plt.legend()
if (display):
    plt.show()
else:
    plt.savefig("error.png")
plt.close()
