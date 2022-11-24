import numpy as np
import matplotlib.pyplot as plt
import math
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

# Plot for the first test
data = np.genfromtxt('error')
N = data[:,0]
e_Linf = data[:,1]
e_L2 = data[:,2]

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(e_Linf))
a, b = pars

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(e_L2))
c, d = pars

if (-a > 1.8 and -c > 1.8):
  print('Test completed for the polynomial reconstruction.')
else:
  print('The convergence rate of the polynomial reconstruction is less than second order')

# Plot
scaling_Linf = np.zeros(np.size(N))
scaling_L2 = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling_Linf[i] = np.exp(b)*N[i]**a
  scaling_L2[i] = np.exp(d)*N[i]**c

plt.figure()
plt.plot(N, e_Linf, 'o', label = r'$L_{\inf}$')
plt.plot(N, scaling_Linf, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.plot(N, e_L2, 's', label = r'$L_2$')
plt.plot(N, scaling_L2, 'black', linestyle='solid', label = '$N^{%2.2f}$' % c)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the polynomial reconstruction")
plt.xscale("log")
plt.yscale("log")
plt.xticks(N, N)
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
plt.close()