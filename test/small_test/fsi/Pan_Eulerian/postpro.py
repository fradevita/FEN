import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

####################################################################################################
# Parameters
####################################################################################################
H = 1
Deltap = np.array([1.763e-3, 8.167e-4, 2.377e-4])
mu = np.array([3.2498036e-3, 1.5e-3, 4.2834760e-4])
Uref = Deltap*H**2/12.0/mu
Tref = H/Uref
Re = Uref*H/mu

####################################################################################################
# Load data
####################################################################################################
# simulation data
simData = [pd.read_csv('case1.csv'), pd.read_csv('case2.csv'), pd.read_csv('case3.csv')] 

# reference data
ref = pd.read_csv('reference.csv')

####################################################################################################
# Reference steady state values
####################################################################################################
x2c = [0.2732, 0.2725, 0.2722]
omegac = [-0.05345, -0.05343, -0.05052]

print('        **** Error ***************************************************')
for j in range(3):
    e_x2 = 100*np.abs(x2c[j] - simData[j]["y"].values[-1])/np.abs(x2c[j])
    e_Om = 100*np.abs(omegac[j] - simData[j]["omega"].values[-1])/np.abs(omegac[j])
    print(f'        Re = {Re[j]:0.2f}, |x2c - x2|/x2c = {e_x2:2.2f}%, |Omc - Om2|/Om2c = {e_Om:2.2f}%')
print('        **************************************************************')

####################################################################################################
# Function to apply periodicity to x 
####################################################################################################
def addPeriodicity(x0):
    x1 = np.zeros_like(x0)
    c = 0
    x1[0] = x0[0]
    for i in range(1,len(x0)):
        if (x0[i] < x0[i-1]):
            c = c + 1
        x1[i] = x0[i] + c
    return x1 - 0.5

####################################################################################################
# Plot
####################################################################################################
styles = ['-b', '-r', '-k']
styles2 = ['xb', 'xr', 'xk']
fig, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (11,5))

ax.set_xlim([0, 50])
ax.set_ylim([0.25, 0.43])
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.plot(ref["x1"], ref["Re12.78"], 'ob', fillstyle = 'none')
ax.plot(ref["x1"], ref["Re27.73"], 'or', fillstyle = 'none')
ax.plot(ref["x1"], ref["Re96.74"], 'ok', fillstyle = 'none')
for j in range(3):
    ax.plot(addPeriodicity(simData[j]["x"]), simData[j]["y"], styles[j], label = f'Re = {Re[j]:2.2f}')
ax.legend()
ax.grid(True)

plt.tight_layout()
if (display):
    plt.show()
    print('        Test completd.')
else:
    plt.savefig("plot.png")
    print('        Test completd, check image.')
plt.close()