import numpy as np
import pandas
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

# Read data
pan = np.genfromtxt('pan_12')
data = pandas.read_csv('data/Circle')

l = len(data["t"])
x0 = data["x"]
x2 = data["y"]

# Add periodicity
x1 = np.zeros(l)
c = 0
x1[0] = x0[0]
for i in range(1,l):
    if (x0[i] < x0[i-1]):
        c = c + 1
    x1[i] = x0[i] + c

# Plot
plt.figure()
plt.plot(x1-0.5,x2,label='data')
plt.plot(pan[:, 0], pan[:, 1], 'o', label='Pan')
plt.xlim([0, 50])
plt.ylim([0.25, 0.43])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend()
plt.savefig("plot.png")
if (display): plt.show()
plt.close()
