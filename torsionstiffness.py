import math
import numpy as np

import sectionproperties as SP
import shearcentre as SC

def tosionalstiffness(section_prop):
	T = 1  # unit force
	r = section_prop.Ha / 2  # radius of semi circle
	ha = section_prop.Ha  # height of the aileron
	tspar = section_prop.tspar
	tskin = section_prop.tskin
	l = math.sqrt((ha / 2) ** 2 + (section_prop.Ca - ha / 2) ** 2)  # length of the skin in triangular part
	a1 = math.pi * ha ** 2 / 8  # enclosed areas
	a2 = ha * (section_prop.Ca - r)/2

	A = np.zeros((3, 3))
	b = np.zeros((3, 1))
	A[0, :] = [2*a1, 2*a2, 0]  # T = T1 + T2 = 1
	b[0] = 1
	A[1, :] = [(r*math.pi/tskin + 2*r/tspar) / (2*a1), (-2*r/tspar)/(2*a1), -1]  # redundant shear flow in cell 1
	A[2, :] = [(-2*r/tspar)/(2*a2), (2*l/tskin + 2 * r/tspar)/(2*a2), -1]  # redundant shear flow in cell 2

	result = np.linalg.solve(A,b)
	J = 1/result[-1, 0]
	return result[0, 0], result[1, 0], J
