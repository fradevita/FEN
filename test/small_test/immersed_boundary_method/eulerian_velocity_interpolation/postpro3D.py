import numpy as np
import matplotlib.pyplot as plt
import math
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


data = np.genfromtxt('out.txt')
N = data[:,0]
eumax = data[:,1]
evmax = data[:,2]
ewmax = data[:,2]

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(eumax))
a, b = pars

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(evmax))
c, d = pars

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(ewmax))
e, f = pars

if (-a > 1.8 and -c > 1.8 and -e > 1.8):
  print('3D test completed.')
else:
  print('The convergence rate of the eulerian interpolation is less than second order')
 
# Plot
scalingu = np.zeros(np.size(N))
scalingv = np.zeros(np.size(N))
scalingw = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scalingu[i] = np.exp(b)*N[i]**a
  scalingv[i] = np.exp(d)*N[i]**c
  scalingw[i] = np.exp(d)*N[i]**e

plt.plot(N, eumax, 'o', label = 'e_u')
plt.plot(N, evmax, 's', label = 'e_v')
plt.plot(N, evmax, '^', label = 'e_w')
plt.plot(N, scalingu, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.plot(N, scalingv, 'black', linestyle='solid', label = '$N^{%2.2f}$' % c)
plt.plot(N, scalingw, 'black', linestyle='solid', label = '$N^{%2.2f}$' % e)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the 3D velocity interpolation")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("plot3D.png")
if (display): plt.show()
