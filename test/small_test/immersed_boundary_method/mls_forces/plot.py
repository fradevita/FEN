import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('local_forces.csv')
data.sort_values(["theta"], inplace=True)
ref  = pd.read_csv('ref_forces.csv')
ref.sort_values(["theta"], inplace=True)

fig, ax = plt.subplots(ncols = 2, nrows = 2)

ax[0,0].set_xlabel('l')
ax[0,0].set_ylabel('F')
ax[0,0].grid()
ax[0,0].set_title('Fvx')
ax[0,0].plot(data["theta"], data["Fvx"], '-o', fillstyle = 'none', label = 'MLS')
ax[0,0].plot( ref["theta"],  ref["Fvx"], '-s', fillstyle = 'none', label = 'ref')
ax[0,0].legend()

ax[0,1].set_xlabel('l')
ax[0,1].set_ylabel('F')
ax[0,1].grid()
ax[0,1].set_title('Fvy')
ax[0,1].plot(data["theta"], data["Fvy"], '-o', fillstyle = 'none', label = 'MLS')
ax[0,1].plot( ref["theta"],  ref["Fvy"], '-s', fillstyle = 'none', label = 'ref')
ax[0,1].legend()

ax[1,0].set_xlabel('l')
ax[1,0].set_ylabel('F')
ax[1,0].grid()
ax[1,0].set_title('Fpx')
ax[1,0].plot(data["theta"], data["Fpx"], '-o', fillstyle = 'none', label = 'MLS')
ax[1,0].plot( ref["theta"],  ref["Fpx"], '-s', fillstyle = 'none', label = 'ref')
ax[1,0].legend()

ax[1,1].set_xlabel('l')
ax[1,1].set_ylabel('F')
ax[1,1].grid()
ax[1,1].set_title('Fpy')
ax[1,1].plot(data["theta"], data["Fpy"], '-o', fillstyle = 'none', label = 'MLS')
ax[1,1].plot( ref["theta"],  ref["Fpy"], '-s', fillstyle = 'none', label = 'ref')
ax[1,1].legend()

plt.show()
plt.close()

fig, ax = plt.subplots(ncols = 1, nrows = 1)

ax.plot(data["theta"], data["nx"], '-o', fillstyle = 'none', label = 'nx')
ax.plot(data["theta"], data["ny"], '-s', fillstyle = 'none', label = 'ny')

plt.show()
plt.close()