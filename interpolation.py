import math
import numpy as np
from matplotlib import pyplot as plt


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
			#print("non-interp: ", t)
			load = np.array(t)[:,2].tolist()
			x = np.array(t)[:,0].tolist()
			t = interp(t, n_chord, N_x, l_a, z)
			#print()
			#print("interp: ", t)
			#input()
			new_load = np.array(t)[:,2].tolist()
			new_x = np.array(t)[:,0].tolist()
			plt.plot(x, load, linestyle='-', marker='v', linewidth=1)
			plt.plot(new_x, new_load, linestyle='-', marker='x', linewidth=1)
			plt.show()
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

def interp(data, n, N, l, z_coord, coord=0):
	new_data = []
	coordinates = np.array(data)[:,coord].tolist()
	interp_coordinates = interp_coords(coordinates, n, N, l, coord=0)
	point_a, point_b = coordinates[0], coordinates[1]
	last_coord_b = 1
	points_a, points_b = [point_a], [point_b]
	for new_coord in interp_coordinates:
		if new_coord > point_b:
			if last_coord_b + 1 < len(coordinates):
				point_a, point_b = point_b, coordinates[last_coord_b+1]
				last_coord_b += 1
		load = interp_two_points(new_coord, point_a, point_b, data[last_coord_b-1][2], data[last_coord_b][2])
		new_data.append([new_coord, z_coord, load])
	return new_data

def interp_coords(coordinates, n, N, l, coord=0):
	n = max(2, n)
	mesh_size = N / n
	mesh_coords = [coordinates[0]]
	for i in range(1, n-1):
		mesh_i = i * mesh_size
		if coord == 0:
			mesh_coord = x_coord(mesh_i, N, l)
		else:
			mesh_coord = z_coord(mesh_i, N, l)
		mesh_coords.append(mesh_coord)
	mesh_coords.append(coordinates[-1])
	return mesh_coords

def interp_two_points(new_coord, point_a, point_b, load_a, load_b):
	if point_a == point_b:
		return load_a
	return load_a + (new_coord - point_a) * (load_b - load_a) / (point_b - point_a)

#interp_two_points(3, 1, 3)
#print(z_coord(10, 81, 0.484))
#print(x_coord(4, 41, 1.691))
get_load(n_chord=200)