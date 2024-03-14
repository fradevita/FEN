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
pan = np.genfromtxt('ref_12.txt')
data = pandas.read_csv('out.csv')

l = len(data["t"])
x0 = data["x"].values
x2 = data["y"].values

# Add periodicity
x1 = np.zeros(l)
c = 0
x1[0] = x0[0]
for i in range(1,l):
    if (x0[i] < x0[i-1]):
        c = c + 1
    x1[i] = x0[i] + c

omega_ref = -0.05467
x2c_ref = 0.273
omega = data["omega"].values
print (f'        x2c  : {   x2[-1]:16.8e}, reference value: {  x2c_ref:16.8e}') 
print (f'        omega: {omega[-1]:16.8e}, reference value: {omega_ref:16.8e}') 

# Plot
plt.figure()
plt.plot(x1-0.5,x2,label='data')
plt.plot(pan[:, 0], pan[:, 1], 'o', label='Pan')
plt.xlim([0, 50])
plt.ylim([0.25, 0.43])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend()
if (display):
    plt.show()
    print('        Test completd.')
else:
    plt.savefig("plot.png")
    print('        Test completd, check image.')
plt.close()

