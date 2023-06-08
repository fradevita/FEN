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

# Reference data
ref = np.genfromtxt('Huang.txt')

# Sorting
ref = ref[ref[:,0].argsort()]

# Simulation data
data = np.genfromtxt('out.txt')

plt.figure()
plt.xlim([0, 20])
plt.xlabel('t')
plt.ylim([-0.4, 0.4])
plt.ylabel('d')
plt.plot(ref[:,0], ref[:,1], '-o', label = 'Huang et. al. 2007')
plt.plot(data[:,0], data[:,1], label = 'simulation')
plt.legend()
if (display):
    plt.show()
    print('Test completd.')
else:
    plt.savefig("plot.png")
    print('Test completd, check image.')
