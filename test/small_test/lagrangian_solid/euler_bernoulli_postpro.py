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

# Beam parameters
l = 1.0
w = 0.01
h = 1.0e-3
E = 210.1e+10
I = w*h**3/12.0
rho = 7850
F = -1.0e-4

# Euler Bernoulli solution
def Euler_Bernoulli(x):
    return F*x**2*(3*l-x)/6/E/I

# Number of resolutions
N = [32,64,128,256]
emax = np.zeros(4)

# Cycle over output file
c = 0

for r in N:

    # Load data
    outfile = np.genfromtxt('out_'+str(r).zfill(3)+'.txt')

    # Compute maximum error
    emax[c] = 0.0
    sol = np.zeros_like(outfile[:,1])
    for i in range(len(outfile[:,1])):
        sol[i] = Euler_Bernoulli(outfile[i,0])
        e = abs(outfile[i,1] - Euler_Bernoulli(outfile[i,0]))
        if (e > emax[c]): emax[c] = e
    c = c + 1

scaling1 = np.zeros(len(N))
scaling2 = np.zeros(len(N))
for i in range(len(N)):
    scaling1[i] = 3.0e-5/N[i]
    scaling2[i] = 1.0e-3/N[i]**2

plt.figure()
plt.xlabel('N')
plt.ylabel('e')
plt.loglog(N, emax, 'o')
plt.loglog(N, scaling1, label = 'N^{-1}')
plt.loglog(N, scaling2, label = 'N^{-2}')
plt.legend()
if (display):
    plt.show()
    print('        Test completd.')
else:
    plt.savefig("euler_bernoulli.png")
    print('        Test completd, check image.')
plt.close()
