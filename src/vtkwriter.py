# Nicolò Grilli
# University of Bristol
# 2 Giugno 2023

# VTK file writer class

import meshio
import numpy as np

class VTKwriter:
	
	def __init__(self,filename,mesh,feproblem):
		self.filename = filename # output file name
		self.mesh = mesh
		self.feproblem = feproblem
		
	# create output mesh with displacement nodal field
	def load_output_field(self):
		point_solution_dict = {} # solution dictionary for point output
		for dof in range(self.feproblem.dofs_per_node):
			point_solution_array = []
			for node in range(self.feproblem.number_of_nodes):
				point_solution_array.append(self.feproblem.u[node*self.feproblem.dofs_per_node + dof])
			if (dof == 0):
				point_solution_dict['ux'] = point_solution_array
			elif (dof == 1):
				point_solution_dict['uy'] = point_solution_array
			elif (dof == 2):
				point_solution_dict['utheta'] = point_solution_array
		cell_solution_dict = {} # solution dictionary for cell output
		# an additional square bracket is needed because meshio uses it to separate different element types
		cell_solution_dict['Exx'] = [self.feproblem.Exx]
		cell_solution_dict['Eyy'] = [self.feproblem.Eyy]
		cell_solution_dict['Exy'] = [self.feproblem.Exy]
		cell_solution_dict['Ezz'] = [self.feproblem.Ezz]
		cell_solution_dict['Sxx'] = [self.feproblem.Sxx]
		cell_solution_dict['Syy'] = [self.feproblem.Syy]
		cell_solution_dict['Sxy'] = [self.feproblem.Sxy]
		cell_solution_dict['Szz'] = [self.feproblem.Szz]
		# convert to 3D points to avoid VTK warning
		points3D = np.column_stack([self.mesh.points[:,0], self.mesh.points[:,1], np.zeros(self.mesh.points.shape[0])])
		self.output_mesh = meshio.Mesh(points3D,self.mesh.cells,point_data=point_solution_dict,cell_data=cell_solution_dict)

	def write_output(self):
		self.output_mesh.write(self.filename,file_format="vtk")
		
