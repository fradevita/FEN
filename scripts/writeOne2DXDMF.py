import json
import sys

# Load Grid File
grid_file = open('setup.json')
Grid = json.load(grid_file)

# Setup Grid
Nx = Grid["Grid"]["Nx"]
Ny = Grid["Grid"]["Ny"]
Nz = Grid["Grid"]["Nz"]
Lx = Grid["Grid"]["Lx"]
Ly = Grid["Grid"]["Ly"]
Lz = Grid["Grid"]["Lz"]
dx = Lx/Nx
dy = Ly/Ny
assert(dx == dy)

step = int(sys.argv[1])

# Output XDMF file
filename = '2Dfields.xmf'
f = open(filename, 'w')

# Header for xml file
f.write('''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">
  <Domain>
    <Grid GridType="Uniform">
      <Topology TopologyType="2DCORECTMesh" Dimensions="%d %d"/>
      <Geometry GeometryType="ORIGIN_DXDY">
      <DataItem Name="Origin" Dimensions="2" NumberType="Float" Precision="4" Format="XML">
        0 0
      </DataItem>
      <DataItem Name="Spacing" Dimensions="2" NumberType="Float" Precision="4" Format="XML">
	%f %f
      </DataItem>
    </Geometry>
'''%(Ny, Nx, dx, dy))

# First velocity component
f.write('''\n
    <Attribute Name="U" AttributeType="Scalar" Center="Node">
      <DataItem Dimensions="%d %d" NumberType="Float" Precision="8" Endian="Little" Format="Binary">
        %s
      </DataItem>
    </Attribute>
'''%(Ny, Nx, 'data/vx_'+str(step).zfill(7)+'.raw'))

# Second velocity component
f.write('''\n
    <Attribute Name="V" AttributeType="Scalar" Center="Node">
      <DataItem Dimensions="%d %d" NumberType="Float" Precision="8" Endian="Little" Format="Binary">
        %s
      </DataItem>
    </Attribute>
'''%(Ny, Nx, 'data/vy_'+str(step).zfill(7)+'.raw'))

# Pressure
f.write('''\n
    <Attribute Name="P" AttributeType="Scalar" Center="Node">
      <DataItem Dimensions="%d %d" NumberType="Float" Precision="8" Endian="Little" Format="Binary">
        %s
      </DataItem>
    </Attribute>
'''%(Ny, Nx, 'data/p_'+str(step).zfill(7)+'.raw'))

# End the xmf file
f.write('''
   </Grid>
</Domain>
</Xdmf>
''')
