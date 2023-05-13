# Nicolò Grilli
# University of Bristol
# 19 Febbraio 2023

# A teaching FE code for linear elasticity

import meshio
from node2d import Node2D
from truss2d2 import Truss2D2

mesh = meshio.read("Job-1.inp")

nodo1 = Node2D(0,mesh.points[0][:])
nodo2 = Node2D(1,mesh.points[1][:])

truss = Truss2D2(nodo1,nodo2,1,1)

print(truss.stiffness_matrix)
