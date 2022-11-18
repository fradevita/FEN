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
data = np.genfromtxt('error_NN')
N = data[:,0]
emax = data[:,1]

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(emax))
a, b = pars

if (-a > 1.8):
  print('Test completed for the Neumann-Neumann Poisson solver.')
else:
  print('The convergence rate of the Neumann-Neumann Poisson solver is less than second order')

# Plot
scaling = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.figure()
plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the Poisson Solver with NN BC")
plt.xscale("log")
plt.yscale("log")
plt.xticks(N, N)
plt.legend()
plt.savefig("plot_NN.png")
if (display): plt.show()
plt.close()

# Plot for the second test
data = np.genfromtxt('error_PP')
N = data[:,0]
emax = data[:,1]

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(emax))
a, b = pars

if (-a > 1.8):
  print('Test completed for the Periodic-Periodic Poisson solver.')
else:
  print('The convergence rate of the Periodic-Periodic Poisson solver is less than second order')

# Plot
scaling = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.figure()
plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the Poisson Solver with PP BC")
plt.xscale("log")
plt.yscale("log")
plt.xticks(N, N)
plt.legend()
plt.savefig("plot_PP.png")
if (display): plt.show()
plt.close()

# Plot for the third test
data = np.genfromtxt('error_PN')
N = data[:,0]
emax = data[:,1]

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(emax))
a, b = pars

if (-a > 1.8):
  print('Test completed for the Periodic-Neumann Poisson solver.')
else:
  print('The convergence rate of the Periodic-Neumann Poisson solver is less than second order')

# Plot
scaling = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.figure()
plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the Poisson Solver with PN BC")
plt.xscale("log")
plt.yscale("log")
plt.xticks(N, N)
plt.legend()
plt.savefig("plot_PN.png")
if (display): plt.show()
plt.close()

# Plot for the 3D test
data = np.genfromtxt('error_PPP')
N = data[:,0]
emax = data[:,1]

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(emax))
a, b = pars

if (-a > 1.8):
  print('Test completed for the Periodic-Periodic-Periodic Poisson solver.')
else:
  print('The convergence rate of the Periodic-Periodic Periodic Poisson solver is less than second order')

# Plot
scaling = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling[i] = np.exp(b)*N[i]**a

plt.figure()
plt.plot(N, emax, 'o', label = 'data')
plt.plot(N, scaling, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.xlabel("N")
plt.ylabel("error")
plt.title("Convergence rate of the Poisson Solver with PPP BC")
plt.xscale("log")
plt.yscale("log")
plt.xticks(N, N)
plt.legend()
plt.savefig("plot_PPP.png")
if (display): plt.show()
plt.close()
