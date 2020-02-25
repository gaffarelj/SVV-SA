import math
import numpy as np

import sectionproperties as SP
import shearcentre as SC

section_prop = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
# The following parameters can be used out of section_prop (example: print(section_prop.z_centroid)
# Iyy, Izz, Am, z_centroid, stiff_area, boomcoords, boomcoords_hinge,
# tskin, tspar, hstiff, tstiff, wstiff, Ha, Ca, r, l_topskin, beta, l_spar_to_end
SC.set_sect(section_prop)

def tosionalstiffness():
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
