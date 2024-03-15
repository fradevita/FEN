import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

# Characterstic scales
D = 15.0e-3  # Particle diameter
U = 0.128    # Free fall velocity

###################################################################################################
# Load data
###################################################################################################
# Simulation data
data1 = pd.read_csv('out1.csv')
data2 = pd.read_csv('out2.csv')
data3 = pd.read_csv('out3.csv')
data4 = pd.read_csv('out4.csv')

###################################################################################################
# Plot results
###################################################################################################
fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize = (15,10))

# Position plot
ax[0].set_xlabel('t [s]')
ax[0].set_xlim([0, 5])
ax[0].set_ylim([0., 9.0])
ax[0].set_ylabel('x/d []')
ax[0].plot(data1["t"], data1["x"]/D-0.3,  'm-', label = 'Re = 1.5, simulation')
ax[0].plot(data2["t"], data2["x"]/D-0.3,  'b-', label = 'Re = 4.1, simulation')
ax[0].plot(data3["t"], data3["x"]/D-0.3,  'g-', label = 'Re = 11.6, simulation')
ax[0].plot(data4["t"], data4["x"]/D-0.3,  'r-', label = 'Re = 31.9, simulation')
ax[0].grid()
ax[0].legend()

# Velocity plot
ax[1].set_xlim([0, 5])
ax[1].set_ylim([-0.22, 0.])
ax[1].set_xlabel('t [s]')
ax[1].set_ylabel('u [m/s]')

ax[1].plot(data1["t"], data1["w"],  'm-', label = 'Re = 1.5, simulation')
ax[1].plot(data2["t"], data2["w"],  'b-', label = 'Re = 4.1, simulation')
ax[1].plot(data3["t"], data3["w"],  'g-', label = 'Re = 11.6, simulation')
ax[1].plot(data4["t"], data4["w"],  'r-', label = 'Re = 31.9, simulation')

ax[1].grid()
ax[1].legend()
plt.tight_layout()
if (display):
    plt.show()
else:
    plt.savefig('plot.png', dpi = 300)
plt.close()
