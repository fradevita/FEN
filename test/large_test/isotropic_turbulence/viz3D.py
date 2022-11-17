import numpy as np
import math
import sys
import plotly.graph_objects as go
import matplotlib.pyplot as plt

def get_the_slice(x, y, z, surfacecolor):
    return go.Surface(x=x,
                      y=y,
                      z=z,
                      surfacecolor=surfacecolor,
                      coloraxis='coloraxis')

def get_lims_colors(surfacecolor):# color limits for a slice
    return np.min(surfacecolor), np.max(surfacecolor)

Nx = 128
Ny = 128
Nz = 128
c = 0
#for n in range(0,500000,1000):
for n in range(1,10,1):

    print(n)

    filename = 'data/vx_'+str(n).zfill(7)+'.raw'
    data = np.fromfile(filename)
    u = np.reshape(data, (Nx,Ny,Nz), order = 'F')

    filename = 'data/vy_'+str(n).zfill(7)+'.raw'
    data = np.fromfile(filename)
    v = np.reshape(data, (Nx,Ny,Nz), order = 'F')

    filename = 'data/vz_'+str(n).zfill(7)+'.raw'
    data = np.fromfile(filename)
    w = np.reshape(data, (Nx,Ny,Nz), order = 'F')

    y = np.linspace(0, 2*math.pi, Ny)
    z = np.linspace(0, 2*math.pi, Nz)
    y, z = np.meshgrid(y,z)
    x = np.zeros(y.shape)
    sminz, smaxz = get_lims_colors(u[int(Nz/2),:,:])
    slice_x = get_the_slice(x, y, z, np.transpose(u[0,:,:]))
    
    x = np.linspace(0, 2*math.pi, Nx)
    z = np.linspace(0, 2*math.pi, Nz)
    x, z = np.meshgrid(x,z)
    y = np.zeros(x.shape)
    sminy, smaxy = get_lims_colors(v[:,0,:])
    slice_y = get_the_slice(x, y, z, np.transpose(v[:,0,:]))
    
    x = np.linspace(0, 2*math.pi, Nx)
    y = np.linspace(0, 2*math.pi, Ny)
    x, y = np.meshgrid(x,y)
    z = np.zeros(x.shape)
    sminz, smaxz = get_lims_colors(w[:,:,0])
    slice_z = get_the_slice(x, y, z, np.transpose(w[:,:,0]))
    
    def colorax(vmin, vmax):
        return dict(cmin=vmin,
                    cmax=vmax)

    fig1 = go.Figure(data=[slice_x, slice_y, slice_z])
    fig1.update_layout(
        title_text='Velocity Contour', 
        #title_x=0.5,
        width=700,
        height=700,
        #scene_zaxis_range=[0,2*math.pi], 
        coloraxis=dict(colorscale='BrBG',
                        colorbar_thickness=25,
                        colorbar_len=0.75,
                        **colorax(sminz, smaxz)))
    fig1.write_image('image/img_'+str(c).zfill(4)+'.png')
    #fig1.show()
    c = c +1
