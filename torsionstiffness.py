import math
import numpy as np

import sectionproperties as SP
import shearcentre as SC

section_prop = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
# The following parameters can be used out of section_prop (example: print(section_prop.z_centroid)
# Iyy, Izz, Am, z_centroid, stiff_area, boomcoords, boomcoords_hinge,
# tskin, tspar, hstiff, tstiff, wstiff, Ha, Ca, r, l_topskin, beta, l_spar_to_end
SC.set_sect(section_prop)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
print("shear center:", xi)

print(q1(0.03),"hello")
print(q5(0.03))

def torsionstiffness():
	"""
	Calculate the torsional stiffness using shear flow in left cell
	Formulas:
	G*dtheta/dz = 1/2/A1 * f qds/t
	J = T/(G*dtheta/dz)
	"""
	# 4 regions to calculate for left cell
	T = 1
	r = section_prop.Ha/2
	ha = section_prop.Ha
	tspar = section_prop.tspar
	tskin = section_prop.tskin
	a1 = math.pi*ha**2/8  # enclosed area

	n = 1000
	int_q1 = 0
	int_q2 = 0#(q2(0) + q2(ha / 2)) * ha / 4 / tspar
	int_q5 = 0#(q5(0) + q5(ha / 2)) * ha / 4 / tspar
	int_q6 = 0
	for i in range(n+1):
		int_q1 += q1(i * np.pi / (2 * n)) * np.pi * r / (2 * n) / tskin
		int_q6 += q6(i * np.pi / (2 * n)) * np.pi * r / (2 * n) / tskin
		int_q2 += q2(r - i * r / n) * r / n / tspar
		int_q5 += q5(r - i * r / n) * r / n / tspar
	print(int_q2,(q2(0) + q2(ha / 2)) * ha / 4 / tspar)
	print(int_q5,(q5(0) + q5(ha / 2)) * ha / 4 / tspar)

	G_dtheta_dz = (int_q1+int_q2+int_q5+int_q6)/2/a1
	J = T/G_dtheta_dz
	return J

a = torsionstiffness()
print(a)

