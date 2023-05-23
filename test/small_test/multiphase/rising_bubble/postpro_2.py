import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import sys

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

Nx = 64*2
Ny = 128*2
Nz = 1
Sf = Nx*Ny*Nz

# Create the grid
Lx = 1
Ly = 2
dx = Lx/Nx
dy = Ly/Ny
X = np.arange(0, Lx, dx)
Y = np.arange(0, Ly, dy)

# Read simulation file
data = np.fromfile('vof.dat')
data_raw = data[0:Sf]
vof = np.reshape(data_raw,(Nx,Ny,Nz), order = 'F')

# Read center of massa data
com_sim = np.genfromtxt('output.txt')

# Read reference data
shape_ref = np.genfromtxt('shape_ref_2.txt')
com_ref = np.genfromtxt('com_ref_2.txt')

# Plot
plt.figure(1)
plt.ylim([0.6, 1.4])
plt.xlabel("x")
plt.ylabel("y")
plt.title("Interface shape at T = 3")
plt.contour(X + 0.5*dx, Y + 0.5*dy, np.transpose(vof[:,:,0]), [0.5])
plt.plot(shape_ref[:,0], shape_ref[:, 1], label = 'ref')
plt.legend()
plt.savefig("shape_2.png")
if (display): plt.show()
plt.close(1)

plt.figure(2)
plt.xlim([0, 3])
plt.xlabel("t")
plt.ylabel("$y_c$")
plt.title("Center of mass vertical position")
plt.plot(com_sim[:,0], com_sim[:, 1], label = 'sim')
plt.plot(com_ref[:,0], com_ref[:, 3], label = 'ref')
plt.legend()
plt.savefig("pos_2.png")
if (display): plt.show()
plt.close(2)

plt.figure(3)
plt.xlim([0, 3])
plt.xlabel("x")
plt.ylabel("$v_c$")
plt.title("Center of mass vertical velocity")
plt.plot(com_sim[:,0], com_sim[:, 2], label = 'sim')
plt.plot(com_ref[:,0], com_ref[:, 4], label = 'ref')
plt.legend()
plt.savefig("vel_2.png")
if (display): plt.show()
plt.close(3)
