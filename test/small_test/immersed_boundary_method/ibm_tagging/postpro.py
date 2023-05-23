import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import json

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

##########################################################################################
# READ SETUP FILE
##########################################################################################
setup_file = open('setup.json')
setup = json.load(setup_file)

# Setup Grid
Nx = setup["Grid"]["Nx"]
Ny = setup["Grid"]["Ny"]
Nz = setup["Grid"]["Nz"]
Lx = setup["Grid"]["Lx"]
Ly = setup["Grid"]["Ly"]
Lz = setup["Grid"]["Lz"]
x0 = setup["Grid"]["origin"][0]
y0 = setup["Grid"]["origin"][1]
dx = Lx/Nx
dy = Ly/Ny
assert(dx == dy)
X = np.linspace(dx/2, Lx - dx/2, Nx)
Y = np.linspace(dy/2, Ly - dy/2, Ny)

##########################################################################################
# LOAD DATA FROM FILES
##########################################################################################
data = np.fromfile('tag_c.raw')
tag_c = np.reshape(data, (Nx,Ny), order = 'F')
data = np.fromfile('tag_x.raw')
tag_x = np.reshape(data, (Nx,Ny), order = 'F')
data = np.fromfile('tag_y.raw')
tag_y = np.reshape(data, (Nx,Ny), order = 'F')

##########################################################################################
# DEFINE THE CIRCLE
##########################################################################################
case_file = open('case.json')
case = json.load(case_file)
xc = case["Case"]["xc"]
yc = case["Case"]["yc"]
rc = case["Case"]["rc"]
NN = 100
phi = np.zeros((NN,NN))
xx = np.linspace(0, Lx, NN)
yy = np.linspace(0, Ly, NN)
for j in range(NN):
  for i in range(NN):
    phi[i,j] = math.sqrt( (xx[i] - xc)**2 + (yy[j] - yc)**2) - rc

##########################################################################################
# PLOT
##########################################################################################
fig, ax = plt.subplots(nrows = 1, ncols = 3)

# Cell center points
points = np.zeros((Nx*Ny,2))
c = 0
for j in range(Ny):
    for i in range(Nx):
        points[c,:] = [j, i]
        c = c + 1
ticks = np.zeros(Nx, dtype = int)
for i in range(Nx):
    ticks[i] = i
labels = [str(int(x + 1)) for x in ticks]
ax[0].set_xticks(ticks)
ax[0].set_xticklabels(labels)
ax[0].set_yticks(ticks)
ax[0].set_yticklabels(labels)
ax[0].set_xlim([-0.5, Nx-0.5])
ax[0].set_ylim([-0.5, Ny-0.5])
ax[0].set_xlabel('i')
ax[0].set_ylabel('j')
ax[0].set_title('cell center points')
ax[0].imshow(tag_c, origin = 'lower', cmap = 'bwr')
for c in range(Nx*Ny):
    ax[0].plot(points[c,0], points[c,1], 'o', color = 'black')
ax[0].contour(yy*Ny - 0.5, xx*Nx - 0.5, phi, [0.0], colors = 'black')
ax[0].set_aspect('equal')

# x-face points
ax[1].set_xticks(ticks)
ax[1].set_xticklabels(labels)
ax[1].set_yticks(ticks)
ax[1].set_yticklabels(labels)
ax[1].set_xlim([-0.5, Nx-0.5])
ax[1].set_ylim([-0.5, Ny-0.5])
ax[1].set_xlabel('i')
ax[1].set_ylabel('j')
ax[1].set_title('x-face points')
ax[1].imshow(tag_x, origin = 'lower', cmap = 'bwr')
for c in range(Nx*Ny):
    ax[1].plot(points[c,0], points[c,1], '>', color = 'black')
ax[1].contour(yy*Ny - 0.5, xx*Nx -1, phi, [0.0], colors = 'black')
ax[1].set_aspect('equal')


# y-face points
ax[2].set_xticks(ticks)
ax[2].set_xticklabels(labels)
ax[2].set_yticks(ticks)
ax[2].set_yticklabels(labels)
ax[2].set_xlim([-0.5, Nx-0.5])
ax[2].set_ylim([-0.5, Ny-0.5])
ax[2].set_xlabel('i')
ax[2].set_ylabel('j')
ax[2].set_title('y-face points')
cm = ax[2].imshow(tag_y, origin = 'lower', cmap = 'bwr')
for c in range(Nx*Ny):
    ax[2].plot(points[c,0], points[c,1], '^', color = 'black')
ax[2].contour(yy*Ny - 1, xx*Nx - 0.5, phi, [0.0], colors = 'black')
ax[2].set_aspect('equal')

plt.colorbar(cm, ax=ax, location = 'bottom')
plt.savefig('plot.png')
if (display): plt.show()
plt.close()
