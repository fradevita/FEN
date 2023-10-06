import numpy as np
import matplotlib.pyplot as plt
import math
import pandas

data = pandas.read_csv('data/S')
a = 1.0
x_ref = np.genfromtxt('xref.txt')
plt.figure()
plt.xlim([0, 24])
plt.xlabel(r'$y/a$')
plt.ylabel(r'$x/a$')
plt.plot(data["y"]/a, data["x"]/a, label = r'$a/\Delta x = 32$')
plt.plot(x_ref[:,0], x_ref[:,1], 'o', label = 'Xia et al. (2009)')
plt.legend(loc = 3)
plt.savefig('yx.png')
plt.show()
plt.close()

theta_ref = np.genfromtxt('thetaref.txt')
plt.figure()
plt.xlim([0, 24])
plt.xlabel(r'$y/a$')
plt.ylabel(r'$\theta/\pi$')
plt.plot(data["y"]/a, 0.25 + data["thetaz"]/math.pi, label = r'$a/\Delta x = 32$')
plt.plot(theta_ref[:,0], theta_ref[:,1], 'o', label = 'Xia et al. (2009)')
plt.legend()
plt.savefig('yteta.png')
plt.show()
plt.close()

# Physical scaling
Re = 12.5
Fr = 0.126
a = 0.1
nu = 0.01
Vt = 12.5*nu/a
Vref = Vt
L = 4*a
T = a/Vref
v_ref = np.genfromtxt('vref.txt')
plt.figure()
plt.xlim([0, 1])
plt.ylim([-1.5, 0])
plt.xlabel(r'$t(s)$')
plt.ylabel(r'$V_t(cm/s)$')
plt.plot(data["t"]*T, data["v"]*Vref, label = r'$a/\Delta x = 32$')
plt.plot(v_ref[:,0], v_ref[:,1], 'o', label = 'Xia et al. (2009)')
plt.legend()
plt.savefig('Vt.png')
plt.show()
