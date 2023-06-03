# Nicolò Grilli
# University of Bristol
# 19 Febbraio 2023

# A teaching FE code for linear elasticity

# Requirements: python 3
# pip install meshio
# pip install netCDF4

import meshio
from node2d import Node2D
from truss2d2 import Truss2D2
from feproblem import FEProblem
from material import Material
from abaqusparser import AbaqusParser
from boundaryconditions import BoundaryConditions
from vtkwriter import VTKwriter
import numpy as np

from beam21 import Beam21

mesh = meshio.read("Job-1.inp")

steel = Material(200e3,0.3)

bc = BoundaryConditions()

fe = FEProblem('Beam21',mesh,steel,bc)

fe.calculate_number_of_nodes()
fe.find_element_type()
fe.calculate_number_of_elements()
fe.calculate_number_of_dofs_per_node()
fe.calculate_global_stiffness_matrix()

abaqus_input_file = AbaqusParser('Job-1.inp',bc)
abaqus_input_file.ReadBC()

fe.apply_BC()

print('\n')
print('\n')
print('\n')

nodo1 = Node2D(0,mesh.points[0][:])
nodo2 = Node2D(0,mesh.points[1][:])
nodo3 = Node2D(0,mesh.points[2][:])

print(nodo1.coords)
print(nodo2.coords)
print(nodo3.coords)

print('\n')
print('\n')
print('\n')

elem1 = Beam21([nodo1,nodo2],steel,1,1)
elem2 = Beam21([nodo2,nodo3],steel,1,1)

print(elem1.stiffness_matrix)

print('\n')
print('\n')
print('\n')

fe.solve()

vtksol = VTKwriter('out.vtk',mesh,fe)
vtksol.load_output_field()
vtksol.write_output()


#output_mesh = meshio.Mesh(mesh.points,mesh.cells,point_data={"u": [soluzione[0],soluzione[3],soluzione[6]]},)

#output_mesh.write("out.vtk",file_format="vtk")
