import math
import numpy as np


def get_load(C_a=0.484, l_a=1.691, n_span=None, n_chord=None):
	f = open("data/aerodynamicloadcrj700.dat").readlines()
	data = []
	for i, row in enumerate(f): # row is span-wise
		t = []
		N_z = len(f)
		N_x = len(row.split(","))
		z = z_coord(i, N_z, C_a)
		for j, point in enumerate(row.split(",")):
			x = x_coord(j, N_x, l_a)
			load = float(point)
			t.append([x, z, load])
		if n_chord is not None:
			t = interp(t, n_chord, N_x, l_a)
		data.append(t)
	return data

def comp_theta(i, N):
	return (i - 1) / N * math.pi

def z_coord(i, N_z, C_a):
	theta = comp_theta(i, N_z)
	theta_1 = comp_theta(i + 1, N_z)
	z_i = -0.5 * (C_a / 2 * (1 - math.cos(theta)) + C_a/2 * (1 - math.cos(theta_1)))
	return z_i

def x_coord(i, N_x, l_a):
	theta = comp_theta(i, N_x)
	theta_1 = comp_theta(i + 1, N_x)
	x_i = 0.5 * (l_a / 2 * (1 - math.cos(theta)) + l_a/2 * (1 - math.cos(theta_1)))
	return x_i

def interp(data, n, N, l, coord=0):
	n = max(2, n)
	data = np.array(data)
	line = [data[0]]
	coordinates = data[:,coord].tolist()
	mesh_size = N / n
	mesh_coords = [coordinates[0]]
	for i in range(1, n-1):
		mesh_i = i * mesh_size
		mesh_coord = x_coord(mesh_i, N, l)
		mesh_coords.append(mesh_coord)
	mesh_coords.append(coordinates[-1])
	print("mesh", len(mesh_coords), mesh_coords)
	print("given", len(coordinates), coordinates)
	input()
	for i in range(1, n - 1):
		x, z, load = data[i]
		x_1, z_1, load_1 = data[i+1]
	line.append(data[-1])
	print(line)
	return data

def interp_two_points(n, point_a, point_b):
	line = [point_a]
	for i in range(1, n + 1):

		print(i)

#interp_two_points(3, 1, 3)
#print(z_coord(10, 81, 0.484))
#print(x_coord(4, 41, 1.691))
get_load(n_chord=80)