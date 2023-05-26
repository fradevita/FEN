import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import math
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

data = np.genfromtxt('error.dat')
N    = data[:,0]
eFvx = data[:,1]
eFvy = data[:,2]
eFpx = data[:,3]
eFpy = data[:,4]

# Fit the error
def objective(x, a, b):
  return a*x + b

popt, _ = curve_fit(objective, np.log(N), np.log(eFvx))
a, b = popt

if (-a > 1.8):
  print('        Convergence rate for Fvx is second order.')
else:
  print('        The convergence rate for Fvx is less than second order')

popt, _ = curve_fit(objective, np.log(N), np.log(eFvy))
c, d = popt

if (-c > 1.8):
  print('        Convergence rate for Fvy is second order.')
else:
  print('        Convergence rate for Fvy is less than second order.')

popt, _ = curve_fit(objective, np.log(N), np.log(eFpx))
e, f = popt

if (-e > 1.8):
  print('        Convergence rate for Fpx is second order.')
else:
  print('        Convergence rate for Fpx is less than second order.')

popt, _ = curve_fit(objective, np.log(N), np.log(eFpy))
g, h = popt

if (-g > 1.8):
  print('        Convergence rate for Fpy is second order.')
else:
  print('        Convergence rate for Fpy is less than second order.')

# Plot
scalingFvx = np.zeros(np.size(N))
scalingFvy = np.zeros(np.size(N))
scalingFpx = np.zeros(np.size(N))
scalingFpy = np.zeros(np.size(N))
for j in range(0,np.size(N)):
  scalingFvx[j] = np.exp(b)*N[j]**a
  scalingFvy[j] = np.exp(d)*N[j]**c
  scalingFpx[j] = np.exp(f)*N[j]**e
  scalingFpy[j] = np.exp(h)*N[j]**g

plt.plot(N, eFvx, 'bo', label = 'Fvx')
plt.plot(N, scalingFvx, 'b-', label = '$N^{%2.2f}$' % a)
plt.plot(N, eFvy, 'go', label = 'Fvy')
plt.plot(N, scalingFvy, 'g-', label = '$N^{%2.2f}$' % c)
plt.plot(N, eFpx, 'ro', label = 'Fpx')
plt.plot(N, scalingFpx, 'r-', label = '$N^{%2.2f}$' % e)
plt.plot(N, eFpy, 'yo', label = 'Fpy')
plt.plot(N, scalingFpy, 'y-', label = '$N^{%2.2f}$' % g)
plt.xlabel("N")
plt.ylabel(r"$|e_{max}|$")
plt.title("Convergence rate of the MLS hydrodynamic load computation")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
