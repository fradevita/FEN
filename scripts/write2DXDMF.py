import json
import sys
import glob

# Load setup File
grid_file = 'grid.json'
grid = json.load(open(grid_file))

# Setup Grid
Nx = grid["Grid"]["Nx"]
Ny = grid["Grid"]["Ny"]
Nz = grid["Grid"]["Nz"]
Lx = grid["Grid"]["Lx"]
Ly = grid["Grid"]["Ly"]
Lz = grid["Grid"]["Lz"]
dx = Lx/Nx
dy = Ly/Ny
assert(dx == dy)

# Output XDMF file
filename = '2Dfields.xmf'
f = open(filename, 'w')

# Header for xml file
f.write('''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">
<Domain>
<Grid Name="Box" GridType="Collection" CollectionType="Temporal">
''')

# Loop over the files in the data directory
nf = len(glob.glob('data/vx*'))

iout = int(sys.argv[1])

for n in range(nf):

    step = n*iout
    
    f.write('''
    <!-- time step -->
    <Grid Name="Box %d" GridType="Uniform"> # 
    <Topology TopologyType="2DCORECTMesh" Dimensions="1 %d %d"/>
    <Geometry GeometryType="ORIGIN_DXDYDZ">
      <DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="4" Format="XML">
        0.0 0.0 0.0
      </DataItem>
      <DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="4" Format="XML">
	%f %f %f
      </DataItem>
    </Geometry>
    <Time Value="%d" />
    '''%(n, Ny+1, Nx+1, dx, dx, dy, n))

    # First velocity component
    f.write('''\n
    <Attribute Name="U" AttributeType="Scalar" Center="Cell">
      <DataItem Dimensions="%d %d" NumberType="Float" Precision="8" Endian="Little" Format="Binary">
        %s
      </DataItem>
    </Attribute>
    '''%(Ny, Nx, 'data/vx_'+str(step).zfill(7)+'.raw'))

    # Second velocity component
    f.write('''\n
    <Attribute Name="V" AttributeType="Scalar" Center="Cell">
      <DataItem Dimensions="%d %d" NumberType="Float" Precision="8" Endian="Little" Format="Binary">
        %s
      </DataItem>
    </Attribute>
    '''%(Ny, Nx, 'data/vy_'+str(step).zfill(7)+'.raw'))

    # Pressure
    f.write('''\n
    <Attribute Name="P" AttributeType="Scalar" Center="Cell">
      <DataItem Dimensions="%d %d" NumberType="Float" Precision="8" Endian="Little" Format="Binary">
        %s
      </DataItem>
    </Attribute>
    </Grid>
    '''%(Ny, Nx, 'data/p_'+str(step).zfill(7)+'.raw'))
  
    
# End the xmf file
f.write('''
   </Grid>
</Domain>
</Xdmf>
''')
