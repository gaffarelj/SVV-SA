import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import  torsionstiffness as TS
import stress as STR
import numpy as np
import math

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)

SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)

_, _, J = TS.tosionalstiffness(sect)

#MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
#Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()
#MC.do_plots()
#MC.plot_result(MC.My, "My_b")

x = 0.554
#Mz, My = MC.Mz(x)[0], MC.My(x)[0]
Mz, My = -34645.28897198964, 119630.64930303364
Sz, Sy = 300000, -90000

s = STR.stress(Mz, My, Sz, Sy, sect, q1, q2, q3, q4, q5, q6)
s.section_stress()
s.plot_all()

# FIND MAXIMUM STRESS COORDINATES
#max_s, max_q, max_vm = float("-Inf"), float("-Inf"), float("-Inf")
#max_s_x, max_q_x, max_vm_x = None, None, None
#for x in np.arange(0, 1.691, 0.005):
#	print(f"Computing stress for x = {x}", end="\r")
#	Mz, My = MC.Mz(x), MC.My(x)
#	s = STR.stress(Mz, My, sect, q1, q2, q3, q4, q5, q6)
#	s.section_stress()
#	data = np.array(s.stresses)
#	str = (data[:,2].flatten()).max()
#	data = np.array(s.vm_stresses)
#	vm = (data[:,2].flatten()).max()
#	s = str
#	#s, q, vm = np.fabs(s.stresses).max(), np.fabs(s.shear_flows).max(), np.fabs(s.vm_stresses).max()
#	if s > max_s:
#		max_s = s
#		max_s_x = x
#	#if q > max_q:
#	#	max_q = q
#	#	max_q_x = x
#	if vm > max_vm:
#		max_vm = vm
#		max_vm_x = x

#print()
#print(max_s, max_s_x)
##print(max_q, max_q_x)
#print(max_vm, max_vm_x)