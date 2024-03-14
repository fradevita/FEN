import numpy as np

def add_block(fp, x0: np.ndarray, delta: float):
    # T1
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P1
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P2
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T2
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P2
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P4
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T3
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P4
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P5
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T4
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P5
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P6
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T5
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P6
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P7
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T6
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P7
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P8
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T7
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P8
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P9
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T8
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P9
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P1
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    return

def add_half_block_x(fp, x0: np.ndarray, delta: float):
    # T1
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P1
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P2
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T2
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P2
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P4
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T3
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P4
    fp.writelines(f'        vertex {x0[0] + 2.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P5
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T4
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P5
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P6
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    return

def add_half_block_y(fp, x0: np.ndarray, delta: float):
    # T1
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P1
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P2
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T6
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 0.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P2
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P4
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T7
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P4
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P5
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    # T8
    fp.writelines('facet normal 0.0 0.0 1.0\n')
    fp.writelines('    outer loop\n')
    fp.writelines(f'        vertex {x0[0] + 1.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P5
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 2.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P6
    fp.writelines(f'        vertex {x0[0] + 0.0*delta:16.8f} {x0[1] + 1.0*delta:16.8f} {x0[2] + 0.0*delta:16.8f}\n') # P3
    fp.writelines('    endloop\n')
    fp.writelines('endfacet\n')
    return

Lx = 1.0
Nbx = 12
delta = Lx/Nbx/2.
Nby = 12

fp = open('mesh.stl', 'w')
fp.writelines('solid mesh.stl\n')
for j in range(0,Nby):
    for i in range(Nbx):
        add_block(fp, [0.0 + i*2*delta, 0.0 + j*2*delta, 0.0], delta)
# Add one half block for ghost nodes on x = -dx
#i = 0
#for j in range(Nby):
#       add_half_block_y(fp, [-delta + i*2*delta, 0.0 + (j)*2*delta, 0.0], delta)
fp.writelines('endsolid mesh.stl')
