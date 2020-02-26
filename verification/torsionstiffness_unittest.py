import math
import numpy as np
import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import stress as STR
print('1')
import sys
import os.path
print ('1')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
print('1')

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
print('1')
SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
#xi = -0.007513567161803937
print('1')

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

#calculate torsional stiffness with original input values
res1, res2, J = tosionalstiffness(sect)
print ('sectionproperty:tskin ', sect.tskin)
print ('original torsstif', J)
print (res1, res2)
J_original = J

#skin thickness multiplication by two
sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
sect.tskin = 2*sect.tskin
res1, res2, J = tosionalstiffness(sect)
print ('sectionprop', sect.tskin)
print ('torsstif', J)
discr = ((J-J_original)/J_original)*100
diff = J/J_original
print ('increase see' diff)