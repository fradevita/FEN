import numpy as np
from scipy import fft
import matplotlib.pyplot as plt

# Load Heltai reference solution
Heltai = np.genfromtxt('Heltai.txt')
Heltai = Heltai[np.argsort(Heltai[:,0])]

# Interpolate Heltai solution on a costant spacing interval
N = 1000
T = 3.0
dt = T/(N-1)
time = np.linspace(0., T, N)
zH = np.zeros_like(time)
for i, t in enumerate(time):
    zH[i] = np.interp(t, Heltai[:,0], Heltai[:,1])

# Construct analytical solution
fan = 2.6418
a = 0.338    
sol = np.zeros_like(time)
for i, t in enumerate(time):
    sol[i] = -a*np.cos(fan*2.0*np.pi*t)

# Load my numerical solution
data = np.genfromtxt('data/tip.dat')

# Compare solutions in physical space
plt.figure(figsize = (10, 5))
plt.xlabel('t [s]')
plt.ylabel('a [m]')
plt.title('Solution in phyiscal space')
plt.plot(time, zH, '-k', label = 'Heltai et al.')
plt.plot(time, sol, 'sr', fillstyle = 'none', label = 'Analytical solution')
plt.plot(data[::100,0], data[::100,1], 'ob', fillstyle = 'none', label = 'numerical solution')
plt.legend()
plt.tight_layout()
plt.savefig('Physical_space.png', dpi = 300)
plt.show()
plt.close()

# FFT of the analytical solution
fftsol = fft.fft(sol)
freqsol = fft.fftfreq(N, dt)[:N//2]
print('Analytical frequency: ', fan, 'Hz')

# FFT of the numerical solution
Nnum = len(data[:,0])
fftnum = fft.fft(data[:,1])
freqnum = fft.fftfreq(Nnum, data[1,0] - data[0,0])[:Nnum//2]
print('Numerical frequency: ', freqnum[np.argmax(abs(fftnum[:Nnum//2]))], 'Hz')

# FFT of interpolated Heltai result
fftzH = fft.fft(zH)
freqzH = fft.fftfreq(N, dt)[:N//2]
print('Heltai frequency: ', freqzH[np.argmax(abs(fftzH[:N//2]))], 'Hz')

# Comparison in frequency space
plt.figure(figsize = (10, 5))
plt.xlim([0, 10])
plt.xlabel('f [Hz]')
plt.ylabel('a [m]')
plt.title('Solution in frequency space')
plt.plot(freqsol, np.abs(fftsol[0:N//2])*2.0/N      , label = 'analytical solution')
plt.plot(freqnum, np.abs(fftnum[0:Nnum//2])*2.0/Nnum, label = 'numerical solution')
plt.plot(freqzH , np.abs(fftzH[0:N//2])*2.0/N       , label = 'interpolated Heltai solution')
plt.plot(fan, a, 'or', fillstyle = 'none')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Physical_space.png', dpi = 300)
plt.show()
plt.close()
