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

data = pandas.read_csv('data/C')
l = len(data["t"])
pan = np.genfromtxt('pan_12')
x0 = data["x"].values
print(x0)
x2 = data["y"].values
x1 = np.zeros(l)
c = 0
x1[0] = x0[0]
for i in range(1,l):
    if (x0[i] < x0[i-1]):
        c = c + 1
    x1[i] = x0[i] + c
 
plt.figure()
plt.plot(x1-0.5,x2,label='data')
plt.plot(pan[:, 0], pan[:, 1], 'o', label='Pan')
plt.xlim([0, 50])
plt.ylim([0.25, 0.43])
plt.legend()
plt.savefig("plot.png")
if(display): plt.show()
plt.close()

fig, ax = plt.subplots()
ax.plot(data["t"], data["Fx"], color = 'red')
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$u$", color = 'red')
ax2 = ax.twinx()
ax2.plot(data["t"], data["x"], color = 'blue')
ax2.set_ylabel(r"$v$", color = 'blue')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots()
ax.plot(data["t"], data["x"], color = 'red')
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$x$", color = 'red')
plt.tight_layout()
plt.show()
