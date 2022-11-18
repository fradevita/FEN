import numpy as np
import matplotlib.pyplot as plt
import math
import sys
from scipy.optimize import curve_fit

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

# Define reference circle
xc = 0.5
yc = 0.5
r = 0.25

NN = 100
phi = np.zeros((NN,NN))
xx = np.linspace(0, 1, NN)
yy = np.linspace(0, 1, NN)
for j in range(NN):
  for i in range(NN):
    phi[j,i] = math.sqrt( (xx[i] - xc)**2 + (yy[j] - yc)**2) - r

N = 8
dx = 1/N
x = np.linspace(0+dx/2, 1-dx/2, N)
y = np.linspace(0+dx/2, 1-dx/2, N)
data = np.fromfile('tag_c')
tag_c = np.reshape(data, (N,N), order = 'F')
data = np.fromfile('tag_x')
tag_x = np.reshape(data, (N,N), order = 'F')
data = np.fromfile('tag_y')
tag_y = np.reshape(data, (N,N), order = 'F')


fig, ax = plt.subplots()
ax.imshow(tag_c)
ax.grid()
ax.set_aspect('equal')
plt.show()
plt.close()

# fig, ax = plt.subplots(nrows = 1, ncols = 3)
# ax[0].set_xlim([0.5, N-0.5])
# ax[0].set_ylim([0.5, N-0.5])
# ax[0].set_xticks(np.linspace(0.5, N-1 + 0.5, N))
# ax[0].set_yticks(np.linspace(0.5, N-1 + 0.5, N))
# ax[0].imshow(tag_c)
# ax[0].contour(yy*N, xx*N, phi, [0.0], colors = 'red')
# ax[0].grid()
# ax[0].set_aspect('equal')

# ax[1].set_xticks(np.arange(0, 1.0, 1.0/N))
# ax[1].set_yticks(np.arange(0, 1.0, 1.0/N))
# ax[1].contour(yy, xx, phi, [0.0])
# ax[1].contourf(x, y, tag_x)
# ax[1].grid()
# ax[1].set_aspect('equal')

# ax[2].set_xticks(np.arange(0, 1.0, 1.0/N))
# ax[2].set_yticks(np.arange(0, 1.0, 1.0/N))
# ax[2].contour(yy, xx, phi, [0.0])
# ax[2].contourf(x, y, tag_y)
# ax[2].grid()
# ax[2].set_aspect('equal')

# plt.show()
# plt.close()