import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

# Particle radius
R = 0.125

# Read data
data = pandas.read_csv('out.csv')
pan = np.genfromtxt('pan_12')

l = len(data["t"])
x0 = data["x"].values
y = data["y"].values

x = np.zeros_like(x0)
c = 0
x[0] = x0[0]
for i in range(1,len(x)):
    if (x0[i] < x0[i-1]):
        c = c + 1
    x[i] = x0[i] + c

plt.figure()
plt.plot(x - 0.5, y, label='data')
plt.plot(pan[:, 0], pan[:, 1], 'o', label='Pan')
plt.xlim([0, 50])
plt.ylim([0.25, 0.43])
plt.legend()
plt.savefig("plot.png")
if(display): plt.show()
plt.close()

fig, ax = plt.subplots()
ax.set_xlim([0, 30])
ax.plot(data["t"], x0)
ax.plot(data["t"], x0 + R, '--g')
ax.plot(data["t"], x0 + R + 1./96., '--y')
ax.plot(data["t"], np.ones_like(x0), '-k')
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$x$")
ax2 = ax.twinx()
ax2.plot(data["t"], data["Fx"], '-o', color = 'red')
ax2.set_ylabel(r'$F_x$')
plt.grid()
plt.tight_layout()
plt.show()
