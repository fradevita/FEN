import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

####################################################################################################
#                                       Load Data                                                  #
####################################################################################################
ABC = pd.read_csv('data/ABC.csv')
Huang200Y = pd.read_csv('references/Huang_Re200_Y.csv')
Lee200Y = pd.read_csv('references/Lee_Re200_Y.csv')
F = pd.read_csv('data/Forces.csv')
LeeCd = pd.read_csv('references/Lee_Re200_Cd.csv')
LeeCl = pd.read_csv('references/Lee_Re200_Cl.csv')

####################################################################################################
#                                       Perform FFT                                                #
####################################################################################################
# To perfrom FFT we need equispaced signals, hence we perform interpolation on a given equispaced 
# time array
N = 1000
t = np.linspace(0, 20, N)
y1 = np.interp(t, ABC["t"], ABC["XB"])
y2 = np.interp(t, Huang200Y["t"], Huang200Y["Y"])
y3 = np.interp(t, Lee200Y["t"], Lee200Y["Y"]) 

# We select only the portion of the signal from t = 10 to t = 20
istart = np.argmin(abs(t - 10.))

# From this time, search the first minimum and the last one
for i in range(istart + 1, N-1):
    if y1[i] < y1[i+1] and y1[i] < y1[i-1]:
        i0 = i
        break
for i in range(N-2, istart + 1, -1):
    if y1[i] < y1[i+1] and y1[i] < y1[i-1]:
        i1 = i
        break
t1  = t[i0:i1]
y1e = y1[i0:i1]

for i in range(istart + 1, N-1):
    if y2[i] < y2[i+1] and y2[i] < y2[i-1]:
        i0 = i
        break
for i in range(N-2, istart + 1, -1):
    if y2[i] < y2[i+1] and y2[i] < y2[i-1]:
        i1 = i
        break
t2  = t[i0:i1]
y2e = y2[i0:i1]

for i in range(istart + 1, N-1):
    if y3[i] < y3[i+1] and y3[i] < y3[i-1]:
        i0 = i
        break
for i in range(N-2, istart + 1, -1):
    if y3[i] < y3[i+1] and y3[i] < y3[i-1]:
        i1 = i
        break
t3  = t[i0:i1]
y3e = y3[i0:i1]

# Perform FFT of the signals
fft1 = np.fft.fft(y1e - y1e.mean())*2/len(y1e)
fft2 = np.fft.fft(y2e - y2e.mean())*2/len(y2e)
fft3 = np.fft.fft(y3e - y3e.mean())*2/len(y3e)

# Setup frequency arrays
Tw = t1[-1] - t1[0]
freq1 = np.zeros_like(t1)
for i in range(len(t1)):
    freq1[i] = 2.*np.pi*i/Tw

Tw = t2[-1] - t2[0]
freq2 = np.zeros_like(t2)
for i in range(len(t2)):
    freq2[i] = 2.*np.pi*i/Tw

Tw = t3[-1] - t3[0]
freq3 = np.zeros_like(t3)
for i in range(len(t3)):
    freq3[i] = 2.*np.pi*i/Tw

####################################################################################################
#                                         Make Plot                                                #
####################################################################################################
fig, ax = plt.subplots(ncols = 2, nrows = 2, figsize = (15,10))

ax[0,0].set_title('Physical domain')
ax[0,0].set_xlim([0, 20])
ax[0,0].set_xlabel(r'$tU/L$')
ax[0,0].set_ylabel(r'$x/L$')
ax[0,0].plot(  Lee200Y["t"],   Lee200Y["Y"], 'or', label = 'Lee et al (2014)', fillstyle = 'none')
ax[0,0].plot(Huang200Y["t"], Huang200Y["Y"], 'sk', label = 'Huang & Sung (2010)', fillstyle = 'none')
ax[0,0].plot(      ABC["t"],      ABC["XB"], '-b', label = 'FEN')
ax[0,0].add_patch(Rectangle((10, -0.4), 8.2, 0.8, alpha = 0.2, color = 'grey'))
ax[0,0].grid(True)
ax[0,0].legend()

ax[0,1].set_title('Frequency domain')
ax[0,1].set_xlim([0, 20])
ax[0,1].set_xlabel(r'$\omega$')
ax[0,1].set_ylabel(r'$a(\omega)$')
ax[0,1].plot(freq1, abs(fft1),  '-b', label = 'FEN')
ax[0,1].plot(freq2, abs(fft2), '--k', label = 'Huang & Sung (2010)')
ax[0,1].plot(freq3, abs(fft3), '-.r', label = 'Lee et al (2014)')
ax[0,1].grid(True)
ax[0,1].legend()

ax[1,0].set_xlim([0, 20])
ax[1,0].set_xlabel(r'$tU/L$')
ax[1,0].set_ylim([0, 1.4])
ax[1,0].set_ylabel(r'$C_d$')
ax[1,0].plot(    F["t"],   2*F["Fz"],  '-b', label = 'FEN')
ax[1,0].plot(LeeCd["t"], LeeCd["Cd"], '-.r', label = 'Lee et al (2014)', fillstyle = 'none')
ax[1,0].legend()
ax[1,0].grid(True)

ax[1,1].set_xlim([0, 20])
ax[1,1].set_xlabel(r'$tU/L$')
ax[1,1].set_ylim([-1.5, 1.5])
ax[1,1].set_ylabel(r'$C_l$')
ax[1,1].plot(    F["t"],   2*F["Fx"],  '-b', label = 'FEN')
ax[1,1].plot(LeeCl["t"], LeeCl["Cl"], '-.r', label = 'Lee et al (2014)', fillstyle = 'none')
ax[1,1].legend()
ax[1,1].grid(True)

plt.tight_layout()
plt.show()
plt.close()
