import numpy as np
import matplotlib.pyplot as plt

Heltai = np.genfromtxt('Heltai.txt')
data = np.genfromtxt('data/tip.dat')

plt.figure(figsize = (10,5))
plt.xlabel('t')
plt.xlim([0, 3])
plt.ylabel(r"$\delta$")
plt.ylim([-0.35, 0.35])
plt.plot(Heltai[:,0], Heltai[:,1], 'o', label = 'Heltai et al 2017')
plt.plot(data[:,0], data[:,1], '-', label = 'simulation, theta')
plt.legend()
plt.tight_layout()
plt.savefig('Heltai.png')
plt.show()
plt.close()
