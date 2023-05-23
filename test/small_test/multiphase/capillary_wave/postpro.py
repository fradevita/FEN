import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import math
import sys
import json
import glob

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

prosperetti = np.genfromtxt('prosperetti.csv', delimiter=',')

#########################################################################################
# READ SETUP FILE
#########################################################################################
grid = json.load(open('grid.json'))
Lx = grid["Grid"]["Lx"]
Ly = grid["Grid"]["Ly"]
Lz = grid["Grid"]["Lz"]
x0 = grid["Grid"]["origin"][0]
y0 = grid["Grid"]["origin"][1]

#########################################################################################
# READ CASE FILE
#########################################################################################
case_file = open('case.json')
case = json.load(case_file)
sigma = case["Case"]["sigma"]
l = case["Case"]["lambda"]
rho = case["Case"]["rho"]
k = 2*math.pi/l
omega_0 = math.sqrt(sigma*k**3/(2*rho))
Tref = omega_0
Tmax = 25.0/Tref

# Cycle over resolutions
res = [8, 16, 32, 64]

timestep = np.genfromtxt('timestep.dat')
cc = 0

L1 = np.zeros(4)
L2 = np.zeros(4)
for r in res:

    # Setup Grid
    Nx = r
    Ny = 3*r
    Nz = 1
    dx = Lx/Nx
    dy = Ly/Ny
    assert(dx == dy)
    X = np.linspace(x0 + dx/2, x0 + Lx - dx/2, Nx)
    Y = np.linspace(y0 + dy/2, y0 + Ly - dy/2, Ny)

    dt = timestep[cc]
    filelist = glob.glob('data_'+str(Nx).zfill(2)+'/*')
    nfile = len(filelist)
    max_amp = np.zeros(nfile)
    t = np.zeros(nfile)
    counter = 0
    
    for filename in sorted(filelist):

        # Read simulation file
        step = filename[12:19]
        data = np.fromfile('data_'+str(r).zfill(2)+'/vof_'+step)
        vof = np.reshape(data, (Nx,Ny,Nz), order = 'F')

        # Reconstruct the interface
        amp = np.zeros(X.size)
        for i in range(X.size):
            vofy = np.zeros(Y.size)
            for j in range(Y.size):
                vofy[j] = vof[i,j,0]
                amp[i] = np.interp(0.5,vofy,Y)
        
        # Maximum amplitude
        max_amp[counter] = np.amax(abs(amp))
        t[counter] = int(step)*dt
        counter = counter + 1

    for i in range(nfile):
        e = abs(max_amp[i] - np.interp(t[i]*Tref,prosperetti[:,0],prosperetti[:,1]))
        if (e > L1[cc]): L1[cc] = e
        L2[cc] = L2[cc] + e**2
    L2[cc] = math.sqrt(L2[cc])

    cc = cc + 1

# Compare with prosperetti solution
plt.xlim([0, 25])
plt.ylim([0, 0.01])
plt.xlabel("t")
plt.ylabel("a")
plt.title("Capillary wave test case")
plt.plot(t*Tref, max_amp, 'o', label = 'Simulation')
plt.plot(prosperetti[:, 0], prosperetti[:, 1], label = 'Prosperetti')
plt.legend()
plt.savefig("comparison.png")
if (display): plt.show()
plt.close()

# Convergence rate
scaling1 = np.zeros(4)
scaling2 = np.zeros(4)
i = 0
for r in res:
    scaling1[i] = 0.4/r
    scaling2[i] = 4.0/r**2
    i = i + 1

plt.figure()
plt.xlabel(r'$N$')
plt.ylabel(r'$e$')
plt.loglog(res, L1, 'o', label = r'$L_1$')
plt.loglog(res, L2, 's', label = r'$L_2$')
plt.loglog(res, scaling1, '-', color = 'black', label = r'$N^{-1}$')
plt.loglog(res, scaling2, '--', color = 'black', label = r'$N^{-2}$')
plt.legend()
plt.savefig("convergence.png")
if (display): plt.show()
plt.close()
