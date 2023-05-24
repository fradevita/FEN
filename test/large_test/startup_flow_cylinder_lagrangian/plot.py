import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('data/C')

fig, ax = plt.subplots(ncols = 2, nrows = 2)

ax[0,0].set_xlabel('l')
ax[0,0].set_ylabel('F')
ax[0,0].grid()
ax[0,0].set_title('Fvx')
ax[0,0].plot(data["t"], data["x"], '-o', fillstyle = 'none', label = 'MLS')
ax[0,0].legend()

ax[0,1].set_xlabel('l')
ax[0,1].set_ylabel('F')
ax[0,1].grid()
ax[0,1].set_title('Fvy')
ax[0,1].plot(data["t"], data["y"], '-o', fillstyle = 'none', label = 'MLS')
ax[0,1].legend()

ax[1,0].set_xlabel('l')
ax[1,0].set_ylabel('F')
ax[1,0].grid()
ax[1,0].set_title('Fpx')
ax[1,0].plot(data["t"], data["u"], '-o', fillstyle = 'none', label = 'MLS')
ax[1,0].legend()

ax[1,1].set_xlabel('l')
ax[1,1].set_ylabel('F')
ax[1,1].grid()
ax[1,1].set_title('Fpy')
ax[1,1].plot(data["t"], data["v"], '-o', fillstyle = 'none', label = 'MLS')
ax[1,1].legend()

plt.show()
plt.close()



fig, ax = plt.subplots(ncols = 2, nrows = 2)

ax[0,0].set_xlabel('l')
ax[0,0].set_ylabel('F')
ax[0,0].grid()
ax[0,0].set_title('Fvx')
ax[0,0].plot(data["t"], data["Fsx"], '-o', fillstyle = 'none', label = 'MLS')
ax[0,0].legend()

ax[0,1].set_xlabel('l')
ax[0,1].set_ylabel('F')
ax[0,1].grid()
ax[0,1].set_title('Fvy')
ax[0,1].plot(data["t"], data["Fsy"], '-o', fillstyle = 'none', label = 'MLS')
ax[0,1].legend()

ax[1,0].set_xlabel('l')
ax[1,0].set_ylabel('F')
ax[1,0].grid()
ax[1,0].set_title('Fpx')
ax[1,0].plot(data["t"], data["Fpx"], '-o', fillstyle = 'none', label = 'MLS')
ax[1,0].legend()

ax[1,1].set_xlabel('l')
ax[1,1].set_ylabel('F')
ax[1,1].grid()
ax[1,1].set_title('Fpy')
ax[1,1].plot(data["t"], data["Fpy"], '-o', fillstyle = 'none', label = 'MLS')
ax[1,1].legend()

plt.show()
plt.close()