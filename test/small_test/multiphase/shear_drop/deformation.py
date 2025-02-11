import glob
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from alive_progress import alive_bar

####################################################################################################
# Setup
####################################################################################################
setup = json.load(open('grid.json'))
Nx = setup["Grid"]["Nx"]
Ny = setup["Grid"]["Ny"]
Nz = setup["Grid"]["Nz"]
Lx = setup["Grid"]["Lx"]
Ly = setup["Grid"]["Ly"]
Lz = setup["Grid"]["Lz"]
x0 = setup["Grid"]["origin"][0]
y0 = setup["Grid"]["origin"][1]
dx = Lx/Nx
dy = Ly/Ny
assert(dx == dy)
X = np.linspace(dx/2, Lx - dx/2, Nx)
Y = np.linspace(dy/2, Ly - dy/2, Ny)
dt = 0.122070E-03
iout = 100

####################################################################################################
# Evaluate deformataion parameter for each case
####################################################################################################
for case in range(1,4):

    datadir = './data_'+str(case) + '/'
    filelist = glob.glob(datadir + '/vof_*')
    Tend = len(filelist)*dt*iout
    t = np.linspace(0,Tend,len(filelist))

    D = []
    Dr = []

    with alive_bar(total = len(filelist), title = 'working on case ' + str(case) + ' ...') as bar:
        for file in sorted(filelist):
            #print('working on file: ' + file)
            # Load VoF function
            vof = np.reshape(np.fromfile(file), (Nx,Ny), order = 'F')
            # Set points on the contour at VoF = 0.5
            fig, ax = plt.subplots()
            cs = ax.contour(X, Y, np.transpose(vof), [0.5])
            for item in cs.collections:
                for i in item.get_paths():
                    v = i.vertices
                    x = v[:, 0]
                    y = v[:, 1]

            # Evaluate minor and major axis
            dmax = 0.0
            dmin = 1000.0
            for i in range(len(x)):
                d = math.sqrt((x[i] - 1.0)**2 + (y[i] - 1.0)**2)
                if d > dmax:
                    dmax = d
                    xmax = x[i]
                    ymax = y[i]
                if d < dmin:
                    dmin = d
                    xmin = x[i]
                    ymin = y[i]
            D.append((dmax-dmin)/(dmax + dmin))
            plt.close()
            bar()
    
    DF = pd.DataFrame({'t': t, 'D': D})
    DF.to_csv('case_'+str(case)+'.csv',index = False)
