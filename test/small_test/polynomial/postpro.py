import numpy as np
import matplotlib.pyplot as plt
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

# Fit the error
def objective(x, a, b):
  return a*x + b

# Plot for the first test
data = np.genfromtxt('error_2nd')
N = data[:,0]
e2nd_Linf = data[:,1]
e2nd_L2 = data[:,2]

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(e2nd_Linf))
a, b = pars

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(e2nd_L2))
c, d = pars

# Plot for the first test
data = np.genfromtxt('error_4th')
N = data[:,0]
e4th_Linf = data[:,1]
e4th_L2 = data[:,2]

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(e4th_Linf))
e, f = pars

pars, cov = curve_fit(f = objective, xdata = np.log(N), ydata = np.log(e4th_L2))
g, h = pars

if (-a > 1.8 and -c > 1.8):
  print('Test completed for the polynomial reconstruction.')
else:
  print('The convergence rate of the polynomial reconstruction is less than second order')

# Plot
scaling2nd_Linf = np.zeros(np.size(N))
scaling2nd_L2 = np.zeros(np.size(N))
scaling4th_Linf = np.zeros(np.size(N))
scaling4th_L2 = np.zeros(np.size(N))
for i in range(0,np.size(N)):
  scaling2nd_Linf[i] = np.exp(b)*N[i]**a
  scaling2nd_L2[i] = np.exp(d)*N[i]**c
  scaling4th_Linf[i] = np.exp(f)*N[i]**e
  scaling4th_L2[i] = np.exp(h)*N[i]**g

plt.figure()
plt.plot(N, e2nd_Linf, 'o', label = r'$L_{\infty}, p = 2$')
plt.plot(N, scaling2nd_Linf, 'black', linestyle='solid', label = '$N^{%2.2f}$' % a)
plt.plot(N, e2nd_L2, 's', label = r'$L_2, p = 2$')
plt.plot(N, scaling2nd_L2, 'black', linestyle='dashed', label = '$N^{%2.2f}$' % c)
plt.plot(N, e4th_Linf, 'o', label = r'$L_{\infty}, p = 4$')
plt.plot(N, scaling4th_Linf, 'black', linestyle='dotted', label = '$N^{%2.2f}$' % e)
plt.plot(N, e4th_L2, 's', label = r'$L_2$, p = 4')
plt.plot(N, scaling4th_L2, 'black', linestyle='dashdot', label = '$N^{%2.2f}$' % g)
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