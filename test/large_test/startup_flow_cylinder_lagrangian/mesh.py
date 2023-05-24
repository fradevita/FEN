import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.patches import Ellipse
import sys

# Ellipse axis
a = 0.5
b = a
# Area
h = (a - b)**2/(a + b)**2
Area = math.pi*(a + b)*(1.0 + 3.0*h/(10.0 + math.sqrt(4.0 - 3.0*h)))

# Griglia Euleriana
Lx = 12.0
Nx = 512*12/2.0
delta = Lx/Nx

# Rapporto tra le griglie
ratio = 0.7

# Spaziatura lagrangiana
rl = ratio*delta

# Numero di punti lagrangiani in base al perimetro approssimato
# dell'ellisse
nv = int(Area/rl) + 1 

# Spaziatura lagrangiana corrispondente ad nv intero
rl = Area/nv

# Coordinate finali del centro di massa
xce = 6.0
yce = 6.0
theta = 0.0

# Function to compute intersection between circle and ellipse
def get_intersections(a, b, xc, yc, r):

    # Coefficienti
    C1 = -2*xc
    C2 = -2*yc
    C3 = xc**2 + yc**2 - r**2
    C4 = a**2/b**2
    C5 = -a**2
    C6 = -(1 - C4)/C1
    C7 = -C2/C1
    C8 = -C3/C1 + C5/C1

    # Coefficienti dell'equazione di quarto grado
    p = [C6**2, 2*C7*C6, 2*C6*C8 + C7**2 + C4, 2*C7*C8, C8**2 + C5]

    # Radici del polinomio di quarto grado
    roots = np.roots(p)

    # Controllo quante radici reali ho
    nreal = 0
    ncomplex = 0
    for n in range(len(roots)):
        if (np.iscomplex(roots[n])):
            ncomplex = ncomplex + 1
        else:
            nreal = nreal + 1

    # Per le radici reali calcolo l'intersezione
    points = np.zeros((nreal,2))
    c = 0
    for n in range(len(roots)):
        if (np.iscomplex(roots[n]) == False):
            x = np.real(roots[n])**2*C6 + C7*np.real(roots[n]) + C8
            points[c,0] = x
            points[c,1] = np.real(roots[n])
            c = c + 1

    return points

# Centro della prima circonferenza
xc = a
yc = 0

# Raggio della circonferenza
r = rl
final_points = np.zeros((nv,2))
final_points[0,0] = xc
final_points[0,1] = yc

avanzo = True
c = 1
while(avanzo):
    points = get_intersections(a, b, xc, yc, r)

    # Seleziono la radice che mi interessa
    if (points[1,1] > points[0,1]):
        xc = points[1,0]
        yc = points[1,1]
    else:
        xc = points[0,0]
        yc = points[0,1]

    final_points[c,:] = [xc, yc] 
    if (xc < 0.0):
        avanzo = False
    c = c + 1

avanzo = True
while(avanzo):
    points = get_intersections(a, b, xc, yc, r)

    # Seleziono la radice che mi interessa
    if (points[1,1] < points[0,1]):
        xc = points[1,0]
        yc = points[1,1]
    else:
        xc = points[0,0]
        yc = points[0,1]

    final_points[c,:] = [xc, yc] 
    if (yc < 0.0):
        avanzo = False
    c = c + 1

avanzo = True
while(avanzo):
    points = get_intersections(a, b, xc, yc, r)

    # Seleziono la radice che mi interessa
    if (points[1,1] < points[0,1]):
        xc = points[1,0]
        yc = points[1,1]
    else:
        xc = points[0,0]
        yc = points[0,1]

    final_points[c,:] = [xc, yc] 
    if (xc > 0.0):
        avanzo = False
    c = c + 1

avanzo = True
while(avanzo):
    points = get_intersections(a, b, xc, yc, r)

    # Seleziono la radice che mi interessa
    if (points[1,1] > points[0,1]):
        xc = points[1,0]
        yc = points[1,1]
    else:
        xc = points[0,0]
        yc = points[0,1]

    final_points[c,:] = [xc, yc] 
    if (yc > 0.0):
        avanzo = False
    c = c + 1
    if (c > nv-1):
        avanzo = False
    
#fig, ax = plt.subplots()
#ellipse = Ellipse((0, 0), 2*a, 2*b)
#ax.add_artist(ellipse)
#ax.plot(final_points[:,0], final_points[:,1], 'ro-')
#ax.set_xlim([-0.6, 0.6])
#ax.set_ylim([-0.3, 0.3])
#plt.show()
#plt.close()

# Calcolo distanza media tra i punti
d = np.zeros(nv)
dmean = 0.0
rms = 0.0
for n in range(1,nv):
    d[n] = math.sqrt( (final_points[n,0] - final_points[n-1,0])**2 + (final_points[n,1] - final_points[n-1,1])**2 )
    dmean = dmean + d[n]
    rms = rms + d[n]**2
dmean = dmean / nv
rms = rms / nv
#print(dmean, math.sqrt(rms), rl)

# Traslo l'ellisse
for n in range(nv):
    final_points[n,0] = final_points[n,0] + xce
    final_points[n,1] = final_points[n,1] + yce

# # Ruoto l'ellisse
# theta = theta*math.pi/180.
# R = np.zeros((2,2))
# R[0,0] =  math.cos(theta)
# R[0,1] = -math.sin(theta)
# R[1,0] =  math.sin(theta)
# R[1,1] =  math.cos(theta)
# d = np.zeros(2)
# for n in range(nv):
#     d[0] = final_points[n,0] - xce
#     d[1] = final_points[n,1] - yce
#     rot = np.matmul(R,d)
#     final_points[n,0] = xce + rot[0]
#     final_points[n,1] = yce + rot[1]

# Salvo la griglia
np.savetxt('mesh.txt', final_points, delimiter = ' ')

# fig, ax = plt.subplots()
# ellipse = Ellipse((0, 0), 2*a, 2*b)
# ax.add_artist(ellipse)
# ax.plot(final_points[:,0], final_points[:,1], 'ro-')
# plt.show()
# plt.close()
