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

N = [16, 32, 64, 128, 256]
emax = np.zeros(len(N))

for r in range(len(N)):

  Nx = N[r]
  Ny = N[r]
  Sf = Nx*Ny*Nz

  # Read data from file
  data = np.fromfile('u_'+str(N[r]).zfill(3))
  u_raw = data[0:Sf]
  u = np.reshape(u_raw,(Nx,Ny,Nz), order = 'F')

  dx = Lx/Nx
  dy = Ly/Ny
  X = np.arange(0, Lx, dx)
  Y = np.arange(0, Ly, dy)
  Yf = Y + 0.5*dy

  e = np.zeros((np.size(X),np.size(Y)))
  for j in range(0, np.size(Y)):
    for i in range(0, np.size(X)):
      sol = -math.cos(X[i] + dx)*math.sin(Y[j] + 0.5*dy)*math.exp(-2*0.3)
      if (sol > 1.0e-14):
        e[i][j] = abs(sol - u[i,j,0])/sol
      else:
        e[i][j] = 0
  emax[r] = np.amax(e)

# Fit the error
def objective(x, a, b):
  return a*x + b

popt, _ = curve_fit(objective, np.log(N), np.log(emax))
a, b = popt

if (-a > 1.8):
  print('        Test completed.')
else:
  print('        The convergence rate of the Poisson solver is less than second order')

# Plot
scaling = np.zeros(len(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Taylor Green Vortex test case")
plt.xscale("log")
plt.yscale("log")
#plt.xticks(N, N)
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
