from dolfin import Mesh
import os

if not os.path.isfile("cylpipe1.xml"):
    try:
        os.system("gmsh cylpipe.geo -3 -o cylpipe.msh")
        os.system("dolfin-convert cylpipe.msh cylpipe.xml")
        os.system("rm cylpipe.msh")
    except RuntimeError:
        raise "Gmsh is required to run this demo"

# Create a mesh
mesh = Mesh("mesh/cylpipe.xml")
