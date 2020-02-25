import math
from interpolation import get_load
from shearcentre import shear_centre

x1 = 0.149  # m
x2 = 0.554  # m
x3 = 1.541  # m
xa = 0.272   # m
theta = math.radians(26)  # rad
Ca = 0.484  # m
la = 1.691  # m
ha = 0.173  # m
P = 37.9*1000  # N

def my_x_matrix(p,x,my_x):
	my_x_LHS = [-max(x-x1, 0), -max(x-x2, 0), -max(x-x3, 0), 0, 0, 0, -max(x-(x2-xa/2), 0), 0, 0, 0, 0, 0]
	my_x_RHS = [my_x - p*math.cos(theta)*max(x-(x2-xa/2), 0)]
	return my_x_LHS, my_x_RHS
# lhs,rhs = my_x_matrix(P,la,0)
# print(lhs,rhs)


def q_x(x,n_chord=200,n_span=400):  # Simple integration of AD on each chord
	interval_chord = Ca/n_chord
	interval_span = la/n_span

	load = get_load(n_chord=n_chord, n_span=n_span)
	load_chord = {}

	for eachchord in load:
		if abs(eachchord[0][0]) <= x:
			s = 0
			for item in eachchord:
				s += item[-1]
			load_chord[round(abs(item[0]), 5)] = s*interval_chord

	sort_qx = dict(sorted(load_chord.items()))
	return sort_qx
# a = q_x(0.3)
# print(a)

def mz_x_matrix(x,mz_x):
	mz_x_LHS = [0,0,0,-max(x-x1,0), -max(x-x2,0), -max(x-x3,0), 0, 0,0,0,0,0]
	n_chord = 200
	n_span = 400
	load_spanwise = q_x(x,n_chord,n_span)
	interval_span = la / n_span

	ff = 0
	for x_coord in load_spanwise.keys():
		ff += x_coord*load_spanwise[x_coord]*interval_span  # equals to the double integration in formula

	p_part = P*math.sin(theta)*(x2+xa/2)

	mz_x_RHS = [mz_x-ff-p_part]
	return mz_x_LHS, mz_x_RHS
# lhs,rhs = mz_x_matrix(la,0)
# print(lhs,rhs)


def tau_x(x):
	n_chord = 200
	n_span = 400
	interval_chord = Ca/n_chord
	interval_span = la/n_span

	load = get_load(n_chord=n_chord, n_span=n_span)
	load_span = {}
	for eachchord in load:
		for data in eachchord:
			if data[0] <= x:
				if round(data[1],5) in load_span:
					load_span[round(data[1],5)] += data[-1]*interval_span*interval_chord
				else:
					load_span[round(data[1],5)] = data[-1]*interval_span*interval_chord
	sort_taux = dict(sorted(load_span.items()))
	return sort_taux
# a = tau_x()
# print(a)

def torque_x(x,t_x):
	xi = shear_centre(200)
	print(xi)
	lhs = [0,0,0,-xi,-xi,-xi,ha/2,0,0,0,0,0]
	tau = tau_x(x)
	f_tau = 0
	for z_coord in tau.keys():
		arm = (xi+ha/2) -z_coord
		f_tau += tau[z_coord]*arm
	rhs = [t_x-P*math.sin(theta)*xi+P*math.cos(theta)*(ha/2)-f_tau]
	return lhs,rhs

# a = torque_x(la,0)
# print(a)
#
# def v_x(x):








