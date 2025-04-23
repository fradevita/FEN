import sys
import glob
import numpy as np
import pandas as pd
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
from alive_progress import alive_bar

# Extract time data
data = pd.read_csv("energy.csv")
t = data["t"].values
#A = data["A"].values
#V = data["V"].values
Ws = data["Ws"].values
Wb = data["Wb"].values
Wa = data["Wa"].values
WA = data["WTA"].values
Wv = data["WV"].values

# Global potential
W = Ws + Wb + Wa + WA + Wv

dt = 1.0e-7
iout = 100
filelist = sorted(glob.glob('data/*.stl'))

with alive_bar(len(filelist), title = 'Making figures ...') as bar:

    for i, filename in enumerate(filelist):
        fig = plt.figure(figsize = (10,5))
        
        # =============
        # First subplot
        # =============
        ax = fig.add_subplot(1, 2, 1)
        ax.set_xlabel(r'$t$ [s]')
        ax.set_xlim([0, t[-1]])
        ax.set_ylabel(r'$E$ [J]')
        ax.plot(t, Ws     , label = 'Ws'     , linestyle = (0, (3, 5, 1, 5)), linewidth = 2)
        ax.plot(t, Wb     , label = 'Wb'     , linestyle =          'dashed', linewidth = 2)
        ax.plot(t, Wa + WA, label = 'Wa + WA', linestyle =          'dotted', linewidth = 2)
        ax.plot(t, Wv     , label = 'Wv'     , linestyle =         'dashdot', linewidth = 2)
        ax.plot(t, W      , label = 'W'      , linestyle =           'solid', linewidth = 2)
        ax.legend()
        time = i*iout*dt
        Wt = np.interp(time, t, W)
        ax.plot(time, Wt, 'or')
        
        # =============
        # Second subplot
        # =============
        ax = fig.add_subplot(1, 2, 2, projection=   '3d')
        ax.set_xlim([-4.0e-6, 4.0e-6])
        ax.set_ylim([-4.0e-6, 4.0e-6])
        ax.set_zlim([-4.0e-6, 4.0e-6])
        
        # Load the STL files and add the vectors to the plot
        STLmesh = mesh.Mesh.from_file(filelist[i])
        ax.add_collection3d(mplot3d.art3d.Poly3DCollection(STLmesh.vectors, edgecolor = 'k'))

        # Auto scale to the mesh size
        scale = STLmesh.points.flatten()
        ax.auto_scale_xyz(scale, scale, scale)
        
        plt.tight_layout()
        plt.savefig('data/figure_' + str(i).zfill(3) + '.png', dpi = 300)
        plt.close()

        bar()

import os
os.system("cp " + filelist[-1] + " end.stl")