import numpy as np
import matplotlib.pyplot as plt
import pandas
from matplotlib.lines import Line2D
import glob

####################################################################################################
# Load data
####################################################################################################
# Load Basilisk data
refD1b = pandas.read_csv('reference/Re1Ca02b.csv')
refD2b = pandas.read_csv('reference/Re1Ca04b.csv')
refD3b = pandas.read_csv('reference/Re1Ca09b.csv')
sCa02 = np.genfromtxt('reference/shapeRe1Ca02b.txt')
sCa04 = np.genfromtxt('reference/shapeRe1Ca04b.txt')
sCa09 = np.genfromtxt('reference/shapeRe1Ca09b.txt')

# Load my data
C1b1 = pandas.read_csv('case_1.csv')
C2b1 = pandas.read_csv('case_2.csv')
C3b1 = pandas.read_csv('case_3.csv')

Nx, Ny = 64, 64
file1 = sorted(glob.glob('data_1/vof_*'))[-1]
vof1 = np.reshape(np.fromfile(file1), (Nx,Ny), order = 'F')
file2 = sorted(glob.glob('data_2/vof_*'))[-1]
vof2 = np.reshape(np.fromfile(file2), (Nx,Ny), order = 'F')
file3 = sorted(glob.glob('data_3/vof_*'))[-1]
vof3 = np.reshape(np.fromfile(file3), (Nx,Ny), order = 'F')
delta = 2./Nx
x = np.linspace(-1 + delta / 2.,1 - delta / 2.,Nx)
y = np.linspace(-1 + delta / 2.,1 - delta / 2.,Ny)
X, Y = np.meshgrid(x,y)

####################################################################################################
# Plot
####################################################################################################
fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize = (12,6))

# deformation index plot
ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[0].set_xlim([0, 3.5])
ax[0].set_xlabel(r'$t$', fontsize = 20)
ax[0].set_ylim([0, 0.6])
ax[0].set_ylabel(r'$D$', fontsize = 20)

ax[0].plot(refD1b["t"], refD1b["D"],  'or', fillstyle = 'none', ms = 10, label = 'Ca = 0.2, basilisk')
ax[0].plot(  C1b1["t"],   C1b1["D"],  '-', linewidth = 2, color = 'black')

ax[0].plot(refD2b["t"], refD2b["D"],  'ob', fillstyle = 'none', ms = 10, label = 'Ca = 0.4, basilisk')
ax[0].plot(  C2b1["t"],   C2b1["D"],  '-', linewidth = 2, color = 'black')

ax[0].plot(refD3b["t"], refD3b["D"],  'og', fillstyle = 'none', ms = 10, label = 'Ca = 0.9, basilisk')
ax[0].plot(  C3b1["t"],   C3b1["D"],  '-', linewidth = 2, color = 'black', label = r'FEN')

ax[0].legend(loc = 'upper left', fontsize = 14)

# Shape plot
ax[1].set_xlim([-1,1])
ax[1].set_ylim([-1,1])
ax[1].set_aspect('equal')
styles = ['solid', 'solid', 'solid']
colors = ['r', 'b', 'g']
lines = []
ax[1].plot(sCa02[:,0][::8], sCa02[:,1][::8], 'or', fillstyle = 'none')
ax[1].contour(X, Y, np.transpose(vof1), [0.5], colors = 'red')
ax[1].plot(sCa04[:,0][::8], sCa04[:,1][::8], 'ob', fillstyle = 'none')
ax[1].contour(X, Y, np.transpose(vof2), [0.5], colors = 'blue')
ax[1].plot(sCa09[:,0][::8], sCa09[:,1][::8], 'og', fillstyle = 'none')
ax[1].contour(X, Y, np.transpose(vof3), [0.5], colors = 'green')
for i, style in enumerate(styles):
    lines.append(Line2D([0], [0], color = colors[i], linewidth=1, linestyle = style))
labels = ['Ca = 0.2', 'Ca = 0.4', 'Ca = 0.9']
ax[1].legend(lines, labels, fontsize = 14)
ax[1].tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)

plt.tight_layout()
plt.savefig('plot.png')
plt.show()
plt.close()
