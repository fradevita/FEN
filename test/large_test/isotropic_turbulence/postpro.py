import numpy as np
import matplotlib.pyplot as plt

D1 = np.genfromtxt('output.txt')
D2 = np.genfromtxt('isotropic.hit3d')
D3 = np.genfromtxt('isotropic.basilisk')

plt.figure()
plt.xlim([0, 300])
plt.xlabel('t')
plt.yscale('log')
plt.ylabel(r"$K^{'}_e$")
plt.plot(D1[:,0], D1[:,1]    , label = 'data')
plt.plot(D2[:,0], D2[:,2]*1.5, label = 'spectral')
plt.plot(D3[:,0], D3[:,2]    , label = 'basilisk')
plt.legend()
plt.savefig('ke.png')
plt.show()
plt.close()
