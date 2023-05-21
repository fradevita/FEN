import numpy as np
import matplotlib.pyplot as plt

s = np.fromfile('s.raw')
plt.figure()
plt.contourf(s)
plt.show()
plt.close()
