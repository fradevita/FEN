import numpy as np
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

Nx = 16
Ny = 16
Nz = 1
Sf = Nx*Ny*Nz

Lx = 1
Ly = 1
dx = Lx/Nx
dy = Ly/Ny

x0 = 0.5
y0 = 0.5
r0 = 0.433

# Read data
data = np.fromfile('h.dat')
h_raw = data[0:Sf]
h = np.reshape(h_raw,(Nx,Ny,Nz), order='F')

X = np.arange(0, Lx, dx)
Y = np.arange(0, Ly, dy)

# Reference solution
sol = np.zeros((int(Nx),int(Ny)))
for j in range(0,np.size(Y)):
  for i in range(0,np.size(X)):
    sol[i][j] = ( (i+0.5)*dx - x0)**2 + ( (j+0.5)*dy - y0)**2 - r0**2

# Figure
fig, ax1 = plt.subplots()
plt.contour(X+dx/2, Y+dy/2, h[:,:,0] ,[0.5], colors = 'red')
plt.plot(0,0,'red', label = 'reconstructed')
plt.contour(X+dx/2, Y+dy/2, sol,[0], colors='black', linestyles = 'dashed')
plt.plot(0,0,'black', label = 'analytical')

# Label every other level using strings
ax1.set_aspect('equal', adjustable='box')
plt.xlabel('x')
plt.ylabel('y')
plt.title('VoF reconstruction')
plt.xlim([0,1])
plt.ylim([0,1])
plt.grid(True)
plt.legend(loc = 10)
if(display): plt.show()
plt.close()
