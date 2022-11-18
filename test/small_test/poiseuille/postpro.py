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

mu = 1.0
dpdx = 1.0
Nx = 4
Nz = 1
Ly = 1

N = [8, 16, 32, 64]
emax = np.zeros(len(N))

for i in range(len(N)):

  Ny = N[i]
  Sf = Nx*Ny*Nz

  # Read data from file
  data = np.fromfile('u_'+str(N[i]).zfill(3))
  u_raw = data[0:Sf]
  u = np.reshape(u_raw,(Nx,Ny,Nz), order = 'F')

  dy = Ly/Ny
  Y = np.arange(0, Ly, dy)
  Yf = Y + 0.5*dy

  uy = np.zeros(np.size(Y))
  sol = np.zeros(np.size(Y))
  e = np.zeros(np.size(Y))
  for j in range(0, np.size(Y)):
    uy[j] = u[2,j,0]
    sol[j] = -dpdx*(Yf[j]**2 - Yf[j]*Ly)/(2*mu)
    e[j] = abs(uy[j] - sol[j])
  
  emax[i] = np.amax(e)

# Fit the error
def objective(x, a, b):
    return a*x + b

pars, cov = curve_fit(f=objective, xdata=np.log(N), ydata=np.log(emax))
a, b = pars

if (-a > 1.8):
  print('Test completed.')
else:
  print('The convergence rate of the Poisson solver is less than second order')

# Plot
scaling = np.zeros(len(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("Error")
plt.title("Poiseuille test case")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
plt.close()
