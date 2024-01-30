import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

df = pd.read_csv('error_2D.csv')
N = df["N"]

plt.figure()
plt.xlabel("N")
plt.ylabel(r"$|e_{max}|$")
plt.title("Convergence rate of the MLS interpolation")
plt.xscale("log")
plt.yscale("log")
for col in df.columns.drop(['N']):
    data = df[col]
    popt, _ = curve_fit(objective, np.log(N), np.log(data))
    a, b = popt
    if (-a > 1.8):
        print('        Convergence rate for ' + col + ' is second order.')
    else:
        print('        Convergence rate for ' + col + ' is less than second order')

    scaling = np.exp(b)*N**a
    plt.plot(N, scaling, '-o', label = col + ', $N^{%2.2f}$' % a)
plt.legend()
plt.show()
plt.close()



# popt, _ = curve_fit(objective, np.log(N), np.log(emax))
# a, b = popt

# if (-a > 1.8):
#   print('        Convergence rate for f(1) is second order.')
# else:
#   print('        The convergence rate of the mls interpolation is less than second order')

# popt, _ = curve_fit(objective, np.log(N), np.log(exmax))
# c, d = popt

# if (-c > 1.8):
#   print('        Convergence rate for derivative in x with f(2) is second order.')
# else:
#   print('        Convergence rate for derivative in x with f(2) is less than second order.')

# popt, _ = curve_fit(objective, np.log(N), np.log(eymax))
# e, f = popt

# if (-e > 1.8):
#   print('        Convergence rate for derivative in y with f(3) is second order.')
# else:
#   print('        Convergence rate for derivative in y with f(3) is less than second order.')

# popt, _ = curve_fit(objective, np.log(N), np.log(efxmax))
# g, h = popt

# if (-g > 1.8):
#   print('        Convergence rate for derivative in x with finite difference is second order.')
# else:
#   print('        Convergence rate for derivative in x with finite difference is less than second order.')

# popt, _ = curve_fit(objective, np.log(N), np.log(efymax))
# i, l = popt

# if (-i > 1.8):
#   print('        Convergence rate for derivative in y with finite difference is second order.')
# else:
#   print('        Convergence rate for derivative in x with finite difference is less than second order.')
  
# # Plot
# scaling = np.zeros(np.size(N))
# scalingx = np.zeros(np.size(N))
# scalingy = np.zeros(np.size(N))
# scalingfx = np.zeros(np.size(N))
# scalingfy = np.zeros(np.size(N))
# for j in range(0,np.size(N)):
#   scaling[j] = np.exp(b)*N[j]**a
#   scalingx[j] = np.exp(d)*N[j]**c
#   scalingy[j] = np.exp(f)*N[j]**e
#   scalingfx[j] = np.exp(h)*N[j]**g
#   scalingfy[j] = np.exp(l)*N[j]**i

# plt.plot(N, emax, 'bo', label = 'f')
# plt.plot(N, scaling, 'b-', label = '$N^{%2.2f}$' % a)
# plt.plot(N, exmax, 'go', label = 'dfdx')
# plt.plot(N, scalingx, 'g-', label = '$N^{%2.2f}$' % c)
# plt.plot(N, eymax, 'ro', label = 'dfdy')
# plt.plot(N, scalingy, 'r-', label = '$N^{%2.2f}$' % e)
# plt.plot(N, efxmax, 'yo', label = 'dfdx finite difference')
# plt.plot(N, scalingfx, 'y-', label = '$N^{%2.2f}$' % g)
# plt.plot(N, efymax, 'co', label = 'dfdy finite difference')
# plt.plot(N, scalingfy, 'c-', label = '$N^{%2.2f}$' % i)
# plt.xlabel("N")
# plt.ylabel(r"$|e_{max}|$")
# plt.title("Convergence rate of the MLS interpolation")
# plt.xscale("log")
# plt.yscale("log")
# plt.legend()
# plt.savefig("plot.png")
# if (display): plt.show()
