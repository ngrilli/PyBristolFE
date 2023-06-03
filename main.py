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

mesh = meshio.read("Bridge.inp")

steel = Material()

bc = BoundaryConditions()

abaqus_input_file = AbaqusParser('Bridge.inp',steel,bc)
abaqus_input_file.ReadBC()
abaqus_input_file.ReadMaterial()

fe = FEProblem('Truss2D2',mesh,steel,bc)

fe.calculate_number_of_nodes()
fe.find_element_type()
fe.calculate_number_of_elements()
fe.calculate_number_of_dofs_per_node()
fe.calculate_global_stiffness_matrix()
fe.apply_BC()

print('\n')
print('\n')
print('\n')

fe.solve()

vtksol = VTKwriter('Bridge.vtk',mesh,fe)
vtksol.load_output_field()
vtksol.write_output()
