# BL-model-for-obliquely-reflecting-solitary-wave
This repository contains all the files needed to run simulations of incident solitary wave interacting with an oblique wall.
One first needs to define the domain, with the file 'mesh_hor.py' (main parameters to set are the lengths of the incident channel and oblique wall, and the angle of incidence).
Running this file will create the 'horizontal.geo' file, that can be open in Gmsh.
With Gmsh, one can generate a 2D mesh in the domain (mesh -> 2D -> save), which creates the 'horizontal.msh' file.
Finally, one can run the main file 'BL_soliton.py', making sure that the length of the channel, the length of the oblique wall and the angle of incidence are the same as in the 'mesh_hor.py' file. We highly recommend to run this file in parallel to reduce the computational time.
Repository owner: Floriane Gidel Paper authors : Floriane Gidel (mmfg@leeds.ac.uk), Onno Bokhove (O.Bokhove@leeds.ac.uk) and Anna Kalogirou (A.Kalogirou@leeds.ac.uk)
