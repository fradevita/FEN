from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import glob

file_list = sorted(glob.glob('data/S_*.stl'))

for i, filename in enumerate(file_list):
    print(i)
    
    # Create a new plot
    figure = pyplot.figure(figsize = (10,6))
    axes = mplot3d.Axes3D(figure)

    axes.set_xlim3d(left=0, right=1) 
    axes.set_ylim3d(bottom=0, top=0.1) 
    axes.set_zlim3d(bottom=-0.35, top=0.35) 
    axes.set_box_aspect((1, 0.2, 1))

    axes.set_xlabel('x')
    axes.set_ylabel('y')
    axes.set_zlabel('z')

    # ax.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
    # ax.zaxis.set_ticklabels([])

    # for line in ax.xaxis.get_ticklines():
    #     line.set_visible(False)
    for line in axes.yaxis.get_ticklines():
        line.set_visible(False)
        # for line in ax.zaxis.get_ticklines():
        #     line.set_visible(False)

    # Load the STL files and add the vectors to the plot
    #your_mesh = mesh.Mesh.from_file('data/S_0207000.stl')
    your_mesh = mesh.Mesh.from_file(filename)
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))

    # Auto scale to the mesh size
    #scale = your_mesh.points.flatten()
    axes.auto_scale_xyz(1.0, 0.1, 0.1)

    # Show the plot to the screen
    #pyplot.show()
    pyplot.savefig('images/img_'+str(i).zfill(4)+'.png')
    pyplot.close()
