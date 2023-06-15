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
		solution_dict = {} # solution dictionary for output
		for dof in range(self.feproblem.dofs_per_node):
			solution_array = []
			for node in range(self.feproblem.number_of_nodes):
				solution_array.append(self.feproblem.u[node*self.feproblem.dofs_per_node + dof])
			if (dof == 0):
				solution_dict['ux'] = solution_array
			elif (dof == 1):
				solution_dict['uy'] = solution_array
			elif (dof == 2):
				solution_dict['utheta'] = solution_array
		# convert to 3D points to avoid VTK warning
		points3D = np.column_stack([self.mesh.points[:,0], self.mesh.points[:,1], np.zeros(self.mesh.points.shape[0])])
		self.output_mesh = meshio.Mesh(points3D,self.mesh.cells,point_data=solution_dict,)

	def write_output(self):
		self.output_mesh.write(self.filename,file_format="vtk")
		
