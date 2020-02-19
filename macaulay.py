import math
from interpolation import get_load

x1 = 0.149  # m
x2 = 0.554  # m
x3 = 1.541  # m
xa = 0.272   # m
theta = math.radians(26)  # rad
Ca = 0.484  # m
la = 1.691  # m

# my_x matrix ### im not sure about this
# def my_x_matrix(x):
# 	my_x_LHS = [-max(x-x1,0), -max(x-x2,0), -max(x-x3,0), 0, 0, 0, -max(x-(x2-xa/2))]
# 	my_x_RHS = [my_x - P*max(x-(x2-xa/2))
# 	return my_x_LHS,my_x_RHS

# AD load
# def q_x(x):
#
load = get_load(n_chord=, n_span=150)
print(len(load))
for item in load:
	print(len(item))