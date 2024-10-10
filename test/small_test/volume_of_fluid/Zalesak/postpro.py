import sys
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

###############################################################################
# Load data
###############################################################################
ini1 = np.reshape(np.fromfile('ini_050.raw'), ( 50, 50), order = 'F')
ini2 = np.reshape(np.fromfile('ini_100.raw'), (100,100), order = 'F')
ini3 = np.reshape(np.fromfile('ini_200.raw'), (200,200), order = 'F')
fin1 = np.reshape(np.fromfile('fin_050.raw'), ( 50, 50), order = 'F')
fin2 = np.reshape(np.fromfile('fin_100.raw'), (100,100), order = 'F')
fin3 = np.reshape(np.fromfile('fin_200.raw'), (200,200), order = 'F')

# Plot
fig, ax = plt.subplots(ncols = 3, nrows = 1, figsize = (15,5))

X, Y = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
ax[0].set_title('50x50')
ax[0].set_xlim([0.3,0.7])
ax[0].set_ylim([0.55,0.95])
ax[0].set_aspect('equal')
ax[0].grid(True)
ax[0].contour(X, Y, np.transpose(ini1), [0.5], colors = 'black', linestyles = 'dashed')
ax[0].contour(X, Y, np.transpose(fin1), [0.5], colors =   'red', linestyles =  'solid')


X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
ax[1].set_title('100x100')
ax[1].set_xlim([0.3,0.7])
ax[1].set_ylim([0.55,0.95])
ax[1].set_aspect('equal')
ax[1].grid(True)
ax[1].contour(X, Y, np.transpose(ini2), [0.5], colors = 'black', linestyles = 'dashed')
ax[1].contour(X, Y, np.transpose(fin2), [0.5], colors =   'red', linestyles =  'solid')

X, Y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, 1, 200))
ax[2].set_title('200x200')
ax[2].set_xlim([0.3,0.7])
ax[2].set_ylim([0.55,0.95])
ax[2].set_aspect('equal')
ax[2].grid(True)
ax[2].contour(X, Y, np.transpose(ini3), [0.5], colors = 'black', linestyles = 'dashed')
ax[2].contour(X, Y, np.transpose(fin3), [0.5], colors =   'red', linestyles =  'solid')

plt.tight_layout()
plt.savefig("fig.png")
if (display): plt.show()
plt.close(fig)
