import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_load(C_a=0.484, l_a=1.691, n_span=None, n_chord=None):
	f = open("data/aerodynamicloadcrj700.dat").readlines()
	data = []
	for i, row in enumerate(f): # row is span-wise
		t = []
		N_z, N_x = len(f), len(row.split(","))
		z = comp_coord(i, N_z, C_a, "z")
		for j, point in enumerate(row.split(",")):
			x = comp_coord(j, N_x, l_a, "x")
			load = float(point)
			t.append([x, z, load])
		if n_chord is not None:
			t = interp(t, n_chord, N_x, l_a, z, "x")
		data.append(t)
	plot_data(data)
	return data

def plot_data(data):
	fig = plt.figure()
	data = np.array(data)
	x, z, load = data[:,:,0].flatten(), data[:,:,1].flatten(), data[:,:,2].flatten()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, z, load, marker="x")
	ax.set_xlabel('X')
	ax.set_ylabel('Z')
	ax.set_zlabel('Load')
	plt.show()

def comp_theta(i, N):
	"""
	Compute angle theta, as specified in the Assignment PDF.
	"""
	return (i - 1) / N * math.pi

def comp_coord(i, N, l, c_type):
	"""
	Compute the coordinate, accordingly to the Assignment PDF.
	Inputs:
		- i: index of the row/column
		- N: total number of rows/columns
		- l: chord/span lenght
	Both computation for x and z coordinates have been combined in this function.
	"""
	theta = comp_theta(i, N)
	theta_1 = comp_theta(i + 1, N)
	f = -1 if c_type == "z" else 1
	coord = f / 2 * (l / 2 * (1 - math.cos(theta)) + l / 2 * (1 - math.cos(theta_1)))
	return coord

def interp(data, n, N, l, other_coord, c_type):
	new_data = []
	coord = 0 if c_type == "x" else 1
	coordinates = np.array(data)[:,coord].tolist()
	interp_coordinates = interp_coords(coordinates, n, N, l, c_type=c_type)
	point_a, point_b = coordinates[0], coordinates[1]
	last_coord_b = 1
	points_a, points_b = [point_a], [point_b]
	for new_coord in interp_coordinates:
		if new_coord > point_b:
			if last_coord_b + 1 < len(coordinates):
				point_a, point_b = point_b, coordinates[last_coord_b+1]
				last_coord_b += 1
		load = interp_two_points(new_coord, point_a, point_b, data[last_coord_b-1][2], data[last_coord_b][2])
		if c_type == "z":
			new_data.append([new_coord, other_coord, load])
		else:
			new_data.append([other_coord, new_coord, load])
	return new_data

def interp_coords(coordinates, n, N, l, c_type):
	"""
	Interpolate the new coordinates for interpolation (in one direction), between the existing ones.
	Inputs:
		- coordinates: existing coordinates
		- n: number of new coordinates specified for interpolation
		- N: number of existing coordinates
		- l: length of the span or chord
		- c_type: coordinate type/direction: "x" or "z"
	"""
	n = max(2, n)						# Make sure there is at least two points for the interpolation
	mesh_size = N / n					# Compute the mesh size: existing number of coord. / specified for interp.
	mesh_coords = [coordinates[0]]		# Set first interp. coord. as the first existing one
	for i in range(1, n-1):				# Compute new interp. coord. 
		mesh_i = i * mesh_size			# Fake row/column index = mesh size * step
		mesh_coord = comp_coord(mesh_i, N, l, c_type)	# Compute the actual coordinate, as from row/column index
		mesh_coords.append(mesh_coord)
	mesh_coords.append(coordinates[-1])	# Set last interp. coord. as the last existing one
	return mesh_coords

def interp_two_points(new_coord, point_a, point_b, load_a, load_b):
	"""
	Interpolate the load between two points.
	Point a is before the new coordinate, point b after.
	"""
	if point_a == point_b:	# Make sure point a is different to point b (zero division)
		return load_a
	# The load at the new coordinate is interpolated by fitting a line to the known
	# load before and after that new coordinate
	return load_a + (new_coord - point_a) * (load_b - load_a) / (point_b - point_a)

get_load(n_chord=100)