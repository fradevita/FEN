import sys
import numpy as np
import glob
import matplotlib.pyplot as plt

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1


dir_list = sorted(glob.glob('data_*'))

W = np.zeros((6,3))
D = np.zeros((6,3))
ks = np.linspace(1.0, 1.5, 6)

for dir in dir_list:
    dat = np.genfromtxt(dir + '/energy.dat')
    i1 = int(dir[-1])
    i2 = int(dir[-3])
    W[i2,i1] = dat[-1,2]
    dat = np.genfromtxt(dir + '/D.dat')
    D[i2,i1] = dat[-1]

Eka0 = np.genfromtxt('E_ka0.txt')
Eka2 = np.genfromtxt('E_ka2.txt')
Eka4 = np.genfromtxt('E_ka4.txt')

L0 = 20.0e-3
D0 = 5.0e-3
t = 1.0e-3
V = L0*D0*t
epsx = 0.02
factor = 2.0/epsx**2/V

fig, ax = plt.subplots(figsize = (10,5))

ax.set_xlabel('Ks [N/m]')
ax.set_ylabel('Elastic energy, W [J]')
ax.plot(ks, W[:,0], '-o')#, label = 'k_a = 0.0 N/m')
ax.plot(ks, W[:,1], '-^')#, label = 'K_a = 0.2 N/m')
ax.plot(ks, W[:,2], '-s')#, label = 'K_a = 0.4 N/m')

ax2 = ax.twinx()
ax2.set_ylabel("Equivalent Young's modulus [Pa]")
ax2.plot(ks, W[:,0]*factor, 'o', label = 'k_a = 0.0 N/m')
ax2.plot(ks, W[:,1]*factor, '^', label = 'K_a = 0.2 N/m')
ax2.plot(ks, W[:,2]*factor, 's', label = 'K_a = 0.4 N/m')

ax2.plot(Eka0[:,0], Eka0[:,1], 'ok', fillstyle = 'none', label = 'Tanaka et al (2012), K_a = 0.0 N/m')
ax2.plot(Eka2[:,0], Eka2[:,1], '^k', fillstyle = 'none', label = 'Tanaka et al (2012), K_a = 0.2 N/m')
ax2.plot(Eka4[:,0], Eka4[:,1], 'sk', fillstyle = 'none', label = 'Tanaka et al (2012), K_a = 0.4 N/m')

plt.legend()
plt.tight_layout()
if (display):
    plt.show()
else:
    plt.savefig("energy.png")
plt.show()
plt.close()

Dka0 = np.genfromtxt('D_ka0.txt')
Dka2 = np.genfromtxt('D_ka2.txt')
Dka4 = np.genfromtxt('D_ka4.txt')

plt.figure(figsize = (10,5))
plt.xlabel('Ks [N/m]')
plt.ylabel("Poisson'r ratio")
plt.ylim([0.15, 0.5])
plt.plot(ks, D[:,0]/epsx, '-o', label = 'k_a = 0.0 N/m')
plt.plot(ks, D[:,1]/epsx, '-^', label = 'K_a = 0.2 N/m')
plt.plot(ks, D[:,2]/epsx, '-s', label = 'K_a = 0.4 N/m')
plt.plot(Dka0[:,0], Dka0[:,1], 'ok', fillstyle = 'none', label = 'Tanaka et al (2012), K_a = 0.0 N/m')
plt.plot(Dka2[:,0], Dka2[:,1], '^k', fillstyle = 'none', label = 'Tanaka et al (2012), K_a = 0.2 N/m')
plt.plot(Dka4[:,0], Dka4[:,1], 'sk', fillstyle = 'none', label = 'Tanaka et al (2012), K_a = 0.4 N/m')
plt.legend()
plt.tight_layout()
if (display):
    plt.show()
    print('        Test completd.')
else:
    plt.savefig("Poisson.png")
    print('        Test completd, check image.')
plt.close()
