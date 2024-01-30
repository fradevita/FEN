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

df = pd.read_csv('error_3D.csv')
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