# Nicolò Grilli
# University of Bristol
# 14 Maggio 2023

# a material
# includes also beam characteristics

class Material:
	
	def __init__(self):
		self.young_modulus = 0.0
		self.poisson_ratio = 0.0
		self.beam_cross_section = 0.0
		self.moment_of_area = 0.0
