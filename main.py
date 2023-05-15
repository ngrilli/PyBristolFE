# Nicolò Grilli
# University of Bristol
# 19 Febbraio 2023

# A teaching FE code for linear elasticity

import meshio
from node2d import Node2D
from truss2d2 import Truss2D2
from feproblem import FEProblem
from material import Material
import numpy as np

from beam21 import Beam21

mesh = meshio.read("Job-1.inp")

steel = Material(200e3,0.3)

fe = FEProblem('Beam21',mesh,steel,'None')

fe.calculate_number_of_nodes()
fe.find_element_type()
fe.calculate_number_of_elements()
fe.calculate_number_of_dofs_per_node()
fe.calculate_global_stiffness_matrix()

print(np.matrix(fe.K))

print('\n')
print('\n')
print('\n')

nodo1 = Node2D(0,mesh.points[0][:])
nodo2 = Node2D(0,mesh.points[1][:])
nodo3 = Node2D(0,mesh.points[2][:])

print(nodo1.coords)
print(nodo2.coords)
print(nodo3.coords)

elem1 = Beam21([nodo1,nodo2],steel,1,1)
elem2 = Beam21([nodo2,nodo3],steel,1,1)

print(elem1.stiffness_matrix)


