import stress as STR
import numpy as np
from matplotlib import pyplot as plt

def get_max(Mz_f, My_f, Sz_f, Sy_f, sect, q1, q2, q3, q4, q5, q6, do_plot=False, dx=0.001):
	max_s, max_q, max_vm = float("-Inf"), float("-Inf"), float("-Inf")
	max_s_x, max_q_x, max_vm_x = None, None, None
	xs = np.arange(0, 1.691, dx)
	stresses, qs, vms = [], [], []
	for x in xs:
		print(f"Computing stress for x = {round(x, 3)}", end="\r")
		Mz, My = Mz_f(x), My_f(x)
		Sz, Sy = Sz_f(x), Sy_f(x)
		s = STR.stress(Mz, My, Sz, Sy, sect, q1, q2, q3, q4, q5, q6)
		s.section_stress()
		data = np.array(s.stresses)
		str = (data[:,2].flatten()).max()
		data = np.array(s.vm_stresses)
		vm = (data[:,2].flatten()).max()
		data = np.array(s.shear_flows)
		q = (data[:,2].flatten()).max()
		stresses.append(str)
		vms.append(vm)
		qs.append(q)
		s = str
		if s > max_s:
			max_s = s
			max_s_x = x
		if q > max_q:
			max_q = q
			max_q_x = x
		if vm > max_vm:
			max_vm = vm
			max_vm_x = x
	print()
	if do_plot:
		fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
		ax1.plot(xs, stresses)
		ax1.set_ylabel("Normal stress [Pa]")
		ax1.set_xlabel("x [m]")
		ax2.plot(xs, qs)
		ax2.set_ylabel("Shear flow [N/m]")
		ax2.set_xlabel("x [m]")
		ax3.plot(xs, vms)
		ax3.set_ylabel("von Mises stress [Pa]")
		ax3.set_xlabel("x [m]")
		plt.savefig(f"plots/stresses/stress-along-x.pdf")
		plt.show()
	return max_s, max_s_x, max_q, max_q_x, max_vm, max_vm_x