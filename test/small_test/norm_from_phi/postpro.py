import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from scipy.optimize import curve_fit

if (len(sys.argv) > 1):
    if (sys.argv[1] == 'silent'):
        display = 0
    else:
        print('wrong option in the run_all.sh script.')
        display = 0
else:
    display = 1

###################################################################################################
# Interface reconstruction plot
###################################################################################################
f2 = glob.glob('fort.2*')
f4 = glob.glob('fort.4*')
assert(len(f2) == len(f4))

# Create a high resolution surface 
N = 200
delta = 1./float(N)
x = np.linspace(delta/2., 1.0 - delta/2., N)
y = np.linspace(delta/2., 1.0 - delta/2., N)
phi = np.zeros((N,N), order = 'F')
for j in range(N):
    for i in range(N):
        phi[i,j] = 0.368 - np.sqrt( (x[i] - 0.525)**2 + (y[j] - 0.464)**2  )

fig, axes = plt.subplots(ncols = 2, nrows = 1)
for ax in axes: 
    ax.set_xlabel(r'$x$')
    ax.set_xticks(np.linspace(0, 1, 17))
    ax.set_xticklabels('')
    ax.set_yticks(np.linspace(0, 1, 17))
    ax.set_ylabel(r'$y$')
    ax.set_yticklabels('')
    ax.contour(x, y, np.transpose(phi), [0], colors = 'black', linewidths = 4)

for id in range(len(f2)):

    f = np.genfromtxt(sorted(f2)[id])
    pr = np.zeros((5,5))
    xr = np.zeros((5))
    yr = np.zeros((5))
    xr = f[0:5, 0]
    yr = f[::5,1]
    c = 0
    for j in range(5):
        for i in range(5):
            pr[i,j] = f[c, 2]
            c += 1
    axes[0].contour(xr, yr, np.transpose(pr), [0.], colors = 'red')
        
    f = np.genfromtxt(sorted(f4)[id])
    pr = np.zeros((5,5))
    xr = f[0:5, 0]
    yr = f[::5,1]
    c = 0
    for j in range(5):
        for i in range(5):
            pr[i,j] = f[c, 2]
            c += 1
    axes[1].contour(xr, yr, np.transpose(pr), [0.], colors = 'red')

for i, ax in enumerate(axes):
    ax.set_title('order = {:1d}'.format(((i+1)*2)))
    ax.set_aspect('equal')
    ax.grid()
plt.savefig('fig1.png')
if (display): plt.show()
plt.close()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
X, Y = np.meshgrid(x,y)
ax.set_ylabel(r'$y$')
ax.set_xlabel(r'$x$')
ax.set_zlabel(r'$\phi$')
ax.set_xlim([0.5, 1.])
ax.set_ylim([0.5, 1.])
ax.plot_wireframe(X, Y, np.transpose(phi))
for id in range(len(f2)):
    f = np.genfromtxt(sorted(f2)[id])
    pr = np.zeros((5,5))
    xr = np.zeros((5))
    yr = np.zeros((5))
    xr = f[0:5, 0]
    yr = f[::5,1]
    c = 0
    for j in range(5):
        for i in range(5):
            pr[i,j] = f[c, 2]
            c += 1
    Xr, Yr = np.meshgrid(xr, yr)
    ax.plot_wireframe(Xr, Yr, np.transpose(pr), colors = 'red')
ax.contour(X, Y, np.transpose(phi), [0.0])
ax.view_init(43, -63)
plt.tight_layout()
plt.savefig('fig2.png')
if (display): plt.show()
plt.close()

###################################################################################################
# Convergence rate plot
###################################################################################################
ecd = np.genfromtxt('ecd.dat')
epr2 = np.genfromtxt('epr2.dat')
epr4 = np.genfromtxt('epr4.dat')

# Fit the error
def objective(x, a, b):
  return a*x + b

pars, cov = curve_fit(f = objective, xdata = np.log(ecd[:,0]), ydata = np.log(ecd[:,1]))
a, b = pars
if (-a > 1.8):
    print('Convergence rate with 2nd order central differencing (Young) is: %4.2f' % abs(a))

pars, cov = curve_fit(f = objective, xdata = np.log(ecd[:,0]), ydata = np.log(epr2[:,1]))
a, b = pars
if (-a > 1.8):
    print('Convergence rate with 2nd order polynomial reconstruction is: %4.2f' % abs(a))

pars, cov = curve_fit(f = objective, xdata = np.log(ecd[:,0]), ydata = np.log(epr4[:,1]))
a, b = pars
if (-a > 1.8):
    print('Convergence rate with 4th order polynomial reconstruction is: %4.2f' % abs(a))


scaling1 = 0.1/ecd[:,0]
scaling2 = 1.0/ecd[:,0]**2
scaling4 = 100.0/ecd[:,0]**4

plt.figure()
plt.xlabel('N')
plt.ylabel('|e|')
plt.loglog(ecd[:,0], ecd[:,1], 'o', label = 'central differencing')
plt.loglog(epr2[:,0], epr2[:,1], 'D', label = 'polynomial reconstruction 2nd order')
plt.loglog(epr4[:,0], epr4[:,1], 'D', label = 'polynomial reconstruction 4th order')
plt.loglog(ecd[:,0], scaling1, '-k', label = '1/N')
plt.loglog(ecd[:,0], scaling2, '--k', label = '1/N^2')
plt.loglog(ecd[:,0], scaling4, '-.k', label = '1/N^4')
plt.legend()
plt.savefig('fig3.png')
if (display): plt.show()
plt.close()

