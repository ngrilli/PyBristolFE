# Nicolò Grilli
# University of Bristol
# 10 Giugno 2023

# a plane strain triangle element
# Notation CPE3 used in Abaqus

import numpy as np

class PlainStrainTriangle(ConstantStrainTriangle):
	
	def __init__(self,nodes):
		super().init(nodes)
