# Nicolò Grilli
# University of Bristol
# 19 Febbraio 2023

# A teaching FE code for linear elasticity

import sys
import meshio
from node2d import Node2D
from truss2d2 import Truss2D2
from feproblem import FEProblem
from material import Material
from abaqusparser import AbaqusParser
from boundaryconditions import BoundaryConditions
from vtkwriter import VTKwriter
import numpy as np

input_file_name = str(sys.argv[1])
input_file_name = input_file_name.lstrip('job=')

mesh = meshio.read(input_file_name)

steel = Material()

bc = BoundaryConditions()

abaqus_input_file = AbaqusParser(input_file_name,steel,bc)
abaqus_input_file.ReadBC()
abaqus_input_file.ReadMaterial()
abaqus_input_file.ReadElementType()

fe = FEProblem(mesh,steel,bc,abaqus_input_file)

fe.calculate_number_of_nodes()
fe.find_element_type()
fe.calculate_number_of_elements()
fe.calculate_number_of_dofs_per_node()
fe.calculate_global_stiffness_matrix()
fe.apply_BC()
fe.solve()

output_file_name = input_file_name.rstrip('.inp')
output_file_name += '.vtk'

vtksol = VTKwriter(output_file_name,mesh,fe)
vtksol.load_output_field()
vtksol.write_output()
