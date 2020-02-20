import math
from interpolation import get_load

x1 = 0.149  # m
x2 = 0.554  # m
x3 = 1.541  # m
xa = 0.272   # m
theta = math.radians(26)  # rad
Ca = 0.484  # m
la = 1.691  # m
ha = 0.173  # m

### im not sure about this
def my_x_matrix(x,p,my_x):
	my_x_LHS = [-max(x-x1,0), -max(x-x2,0), -max(x-x3,0), 0, 0, 0, -max(x-(x2-xa/2))]
	my_x_RHS = [my_x - p*math.cos(theta)*max(x-(x2-xa/2))
	return my_x_LHS,my_x_RHS

# AD load
def q_x():  # Simple integration of AD on each chord
	n_chord = 200
	n_span = 400
	interval_chord = Ca/n_chord

	load = get_load(n_chord=n_chord, n_span=n_span)
	load_chord = {}

	for eachchord in load:
		s = 0
		for item in eachchord:
			s += item[-1]
		load_chord[round(abs(item[0]), 5)] = s*interval_chord

	sort_d = dict(sorted(load_chord.items()))
	return sort_d

a = q_x()
print(a)

def mz_x_matrix(x, mz_x):
	mz_x_LHS = [-max(x-x1,0), -max(x-x2,0), -max(x-x3,0), 0, 0, 0, 0)]
	mz_x_RHS = [mz_x - )
	return mz_x_LHS,mz_x_RHS


