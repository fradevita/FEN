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

emax = np.zeros(5)
data = np.genfromtxt('error')
N = data[:,0]
emax = data[:,1]

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
scaling = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the Multigrid Poisson Solver with PP BC")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
