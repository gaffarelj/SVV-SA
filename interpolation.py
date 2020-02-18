import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_load(C_a=0.484, l_a=1.691, n_span=None, n_chord=None):
	"""
	Get the aerodynamic load.
	Inputs:
		- C_a: chord length
		- l_a: span length
		- n_span: number of points to which to interpolate the load, span-wise
		- n_chord: number of points to which to interpolate the load, chord-wise
	"""
	f = open("data/aerodynamicloadcrj700.dat").readlines()	# Open the aerodynamic load
	data = []
	for i, row in enumerate(f):						# Go trough all rows (span-wise)
		row_d = []									# List to contain data from the row
		N_z, N_x = len(f), len(row.split(","))		# The number of z/x coordinates = number of rows/columns
		if n_chord is None:							# THIS IS A TEMPORARY FIX, NOT IDEAL TO KEEP
			n_chord = N_x
		z = comp_coord(i, N_z, C_a, "z")			# Compute the z coordinate
		for j, point in enumerate(row.split(",")):	# Go trough all columns of the row
			x = comp_coord(j, N_x, l_a, "x")		# Compute the x coordinate
			load = 1000 * float(point)				# Convert the string load to a float, in [N] instead of [kN]
			row_d.append([x, z, load])				# Save the x, z coord. and the load to a temp list
		if n_chord is not None:						# If n_chord had been specified for interpolation...
			# ...interpolate new loads along the x-axis (chord-wise)
			row_d = interp(row_d, n_chord, N_x, l_a, z, "x")
		data.append(row_d)							# Save the data of the row
	if n_span is not None:							# If n_span had been specified for interpolation...
		interp_data = []
		i = 0
		for column in zip(*data):					# Go through every column
			column = np.array(column)
			x, z, load = column[:,0], column[:,1], column[:,2]	# Extract data from the column
			try:
				# Apply linear interpolation on the column
				column_interp = interp(column.tolist(), n_span, N_z, C_a, x[i], "z")
				interp_data.append(column_interp)	# Save the interpolated data
			except IndexError:						# THIS IS A TEMPORARY FIX, NOT IDEAL TO KEEP
				pass
			i += 1
		data = list(zip(*interp_data))
	plot_data(data)									# Plot all data in 3D
	return data

def plot_data(data):
	"""
	Plot the resulting data in 3D
	"""
	fig = plt.figure()
	data = np.array(data) # Convert the list to a numpy array
	# Extract the x, z, and load from the multidimensional data array
	x, z, load = data[:,:,0].flatten(), data[:,:,1].flatten(), data[:,:,2].flatten()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, z, load, marker="x")
	ax.set_xlabel('X [m]')
	ax.set_ylabel('Z [m]')
	ax.set_zlabel('Load [N]')
	ax.set_title('Aerodynamic load over the aileron')
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
	"""
	Interpolate the dataset to n new points
	Inputs:
		- data: given dataset, to be interpolated
		- n: number of points to which to interpolate
		- N: number of rows/columns in the existing dataset
		- other_coord: coordinate list that won't be interpolated
		- c_type: coordinate type/direction: "x" or "z"
	"""
	new_data = []
	coordinates = np.array(data)[:,0].tolist()								# Extract the coordinate (to be interpolated) from the data
	interp_coordinates = interp_coords(coordinates, n, N, l, c_type)		# Interpolate the existing coord. to n new coord.
	point_a, point_b = coordinates[0], coordinates[1]						# Set the first point before and after the new inter. coord. 
	last_coord_b = 1														# Save the index of the point b
	for new_coord in interp_coordinates:									# Go through all interpolated coordinates
		# If the interpolated coord. is above point b...
		if (c_type == "x" and new_coord > point_b) or (c_type == "z" and new_coord < point_b) and last_coord_b + 1 < len(coordinates):
			try:
				point_a, point_b = point_b, coordinates[last_coord_b + 1]	# ...set point a as previous point b, and b as next
				last_coord_b += 1											# Save the index of the new point b
			except IndexError:												# THIS IS A TEMPORARY FIX, NOT IDEAL TO KEEP
				pass
		# Interpolate the load at the interp. coordinate, between the loads at points a and b
		load = interp_two_points(new_coord, point_a, point_b, data[last_coord_b-1][2], data[last_coord_b][2])
		# Save the load, new interp. coord., and unchanged coord. in the new data array
		new_data.append([new_coord, other_coord, load]) if c_type == "z" else new_data.append([other_coord, new_coord, load])
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


get_load(n_chord=None, n_span=150)