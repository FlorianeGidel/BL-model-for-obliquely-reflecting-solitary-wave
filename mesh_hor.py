
# Define the domain on which the soliton propagates. This domain is characterized by an incident channel of length Lc, in which the initial soliton propagates in the x-direction. At the end of this channel, the side wall makes an angle 'psi' with the x-direction. This wall has length Lw. To avoid very sharp angles in the mesh's quadrialeterals, the right-end of the channel is cut with a wall in the y-direction of length 'delta' (see Fig. 5 in the paper). The function 'half_domain' generates a '.geo' file containing the domain boundaries and mesh specifications, that will be compiled in Gmsh to obtain the '.msh' file, with refinements dx and dy in largest areas, and dxx and dyy in finest areas.

from firedrake import *

def half_domain(psi_exp, Ls,  Li, delta, dx, dy, dxx, dyy):
    
    # --- Boundary coordinates --- #
    shift_s = delta/sin(psi_exp)
    x1 = 0.0
    x2 = Li
    x3 = x2 + (Ls-shift_s)*cos(psi_exp)
    x4 = x3
    x5 = x4
    x6 = x2
    x7 = x1
    x8 = x1
    x9 = x6
    
    y1 = 0.0
    y2 = y1
    y3 = y2 + (Ls-shift_s)*sin(psi_exp)
    y4 = y3 + 0.9*delta
    y5 = y3 + delta
    y6 = y5
    y7 = y6
    y8 = y2 + 0.9*delta
    y9 = y8
    
    
    # --- Writing the '.geo' file --- #
    script_file = open('horizontal.geo', 'w')
    script_file.write('%1s%1s%1s\n\n' %('ls = ', ls, ';'))
    ls = 1.0
    
    # Placing the points
    script_file.write('%1s%1s%1s%1s%1s\n' %('Point(1) = {', x1, ', ' , y1, ', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n' %('Point(2) = {',x2,', ', y2 ,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n' %('Point(3) = {',x3,', ', y3,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n' %('Point(4) = {',x4,', ', y4,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n' %('Point(5) = {',x5,', ', y5,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n\n' %('Point(6) = {',x6,', ', y6,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n\n' %('Point(7) = {',x7,', ', y7,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n\n' %('Point(8) = {',x8,', ', y8,', 0.0, ls};'))
    script_file.write('%1s%1s%1s%1s%1s\n\n' %('Point(9) = {',x9,', ', y9,', 0.0, ls};'))
    
    # tracing the lines
    for i in range(1,8):
        script_file.write('%1s%1s%1s%1s%1s%1s%1s\n' %('Line(',i,') = {', i,',', i+1, '};'))
    
    script_file.write('%1s%1s%1s%1s%1s%1s%1s\n' %('Line(',8,') = {',8,',', 1, '};'))
    script_file.write('%1s%1s%1s%1s%1s%1s%1s\n\n' %('Line(',9,') = {', 6,',', 9, '};'))
    script_file.write('%1s%1s%1s%1s%1s%1s%1s\n\n' %('Line(',10,') = {', 9,',', 2, '};'))
    script_file.write('%1s%1s%1s%1s%1s%1s%1s\n\n' %('Line(',11,') = {', 8,',', 9, '};'))
    script_file.write('%1s%1s%1s%1s%1s%1s%1s\n\n' %('Line(',12,') = {', 9,',', 4, '};'))

    # Define the surface
    script_file.write('%1s\n' %('Line Loop(13) = {6, 7, 11, -9};'))
    script_file.write('%1s\n\n' %('Plane Surface(14) = {13};'))
    
    script_file.write('%1s\n' %('Line Loop(15) = {9, 12, 4, 5};'))
    script_file.write('%1s\n\n' %('Plane Surface(16) = {15};'))
    
    script_file.write('%1s\n' %('Line Loop(17) = {8, 1, -10, -11};'))
    script_file.write('%1s\n\n' %('Plane Surface(18) = {17};'))

    script_file.write('%1s\n' %('Line Loop(19) = {12, -3, -2, -10};'))
    script_file.write('%1s\n\n' %('Plane Surface(20) = {19};'))

    # Specify the number of nodes on each boundary
    N = int(Li/dx)
    script_file.write('%1s%1s%1s\n' %('Transfinite Line{1, 6, 11} = ',N, ';'))

    N = int((Ls-shift_s)/dxx)
    script_file.write('%1s%1s%1s\n' %('Transfinite Line{2, 12, 5} = ',N, ';'))

    N = int((y9-y2)/dyy)
    script_file.write('%1s%1s%1s\n' %('Transfinite Line{3, 8, 10} = ',N, ';'))

    N = int((y6-y9)/dy)
    script_file.write('%1s%1s%1s\n' %('Transfinite Line{4, 7, 9} = ',N, ';'))

    # Get quadrilaterals
    script_file.write('%1s\n' %('Transfinite Surface{14, 16, 18, 20};'))
    script_file.write('%1s\n\n' %('Recombine Surface{14, 16, 18, 20};'))
    
    script_file.close()


## PARAMETERS ##
Lc = 5.0        # Length of the channel
Lw = 500        # Length of the wall
delta = 6.0     # Size of the end wall
psi = pi/6      # Angle between the wall and the x-direction
dx = 0.4        # x-refinement in largest areas
dy = 1.5        # y-refinement in largest areas
dxx = 0.25      # x-refinement in finest areas
dyy = 0.25      # y-refinement in finest areas

half_domain(psi, Lw,  Lc, delta, dx, dy, dxx, dyy)
