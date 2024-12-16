import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

data = pd.read_csv('out.csv')

def objective(x, a, b):
    return a*x + b

keys = data.keys()
nkeys = len(keys)

for i in range(1,nkeys):
    pars, conv = curve_fit(f = objective, xdata = np.log(data["N"]), ydata = np.log(data["Linfx"]))
    a, b = pars
    print("      Convergence rate for " + keys[i] + f" = {a:2.2f}")
