import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import  torsionstiffness as TS
import stress as STR
import max_stress as MS
import numpy as np
import math

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)

SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)

_, _, J = TS.torsionalstiffness(sect)

MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()
MC.do_plots()

for x in [0.418, 0.544]:
	Mz, My = MC.Mz(x), MC.My(x)
	Sz, Sy = MC.Sz(x), MC.Sy(x)
	T = MC.T(x)
	s = STR.stress(Mz, My, Sz, Sy, T, sect, q1, q2, q3, q4, q5, q6)
	s.section_stress()
	s.plot_all(x)

input("Press ENTER to find the maximum stress location...")
# Find maximum stress location
max_s, max_s_x, max_q, max_q_x, max_vm, max_vm_x = MS.get_max(MC.Mz, MC.My, MC.Sz, MC.Sy, sect, q1, q2, q3, q4, q5, q6, do_plot=False)
print(max_s, max_s_x)
print(max_q, max_q_x)
print(max_vm, max_vm_x)

max_s, max_q, max_vm = 8.61069683e+08, 512682.06822037, 979999190.1066628
max_s_x, max_q_x, max_vm_x = 0.544, 0.418, 0.418