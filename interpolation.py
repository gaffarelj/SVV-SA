import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_load(C_a=0.484, l_a=1.691, n_span=None, n_chord=None, do_plot=False, fixed_load=None):
	"""
	Get the aerodynamic load.
	Inputs:
		- C_a: chord length
		- l_a: span length
		- n_span: number of points to which to interpolate the load, span-wise
		- n_chord: number of points to which to interpolate the load, chord-wise
	"""
	if fixed_load is not None:
		span_x_coords = interp_coords(list(np.linspace(0, l_a, 20)), n_span, 20, l_a, "x")
		sections_data = []
		for x in span_x_coords:
			sections_data.append([x, [[0.25*C_a, fixed_load]]])
		if do_plot:
			plot_data(sections_data)
		return span_x_coords, sections_data
	

	f = open("data/aerodynamicloadcrj700.dat").readlines()	# Open the aerodynamic load
	data = []
	N_z, N_x = len(f), len(f[0].split(","))			# The number of z/x coordinates = number of rows/columns

	for i, row in enumerate(f):						# Go trough all rows (chord-wise)
		z = comp_coord(i+1, N_z, C_a, "z")			# Compute z-coordinate base on row number
		col = []
		for j, point in enumerate(row.split(",")):	# Go trough every column
			x = comp_coord(j+1, N_x, l_a, "x")		# Compute z-coordinate base on row number
			load = 1000 * float(point)				# Extract the load at that point (in N)
			col.append([z, x, load])				# Save the z and x coords + the load at that point
		data.append(col)

	# Flip the dataset to cut span-wise instead of chord-wise
	span_x_coords = []	# Span-wise x coordinates
	sections_data = []	# z, load per section
	for span_section in zip(*data):	# Flip the dataset here, so go trough each column first
		x = span_section[0][1]		# Extract the x coord. of the section
		span_x_coords.append(x)		# Save the x coord. of the section
		span_cut = []
		for point in span_section:	# For each point on the chord (in this section)
			z, load = point[0], point[2]	# Save the z-coordinate and load on the chord
			span_cut.append([z, load])
		sections_data.append([x, span_cut])	# Save the data as x-coord. : [z-coords:loads]
	# Make sure to make the span-wise interpolation even if it not asked,
	# because the chord-wise interp. is build within the span-wise one
	if n_chord is not None and n_span is None:
		n_span = N_x
	if n_span is not None:
		interp_sections_data = []
		interp_span_x_coords = interp_coords(span_x_coords, n_span, N_x, l_a, "x")	# Interpolate the x-coordinates on n new locations
		point_a, point_b = span_x_coords[0], span_x_coords[1]						# Set the first points before and after the new inter. coord. 
		last_coord_b = 1															# Save the index of the point b
		for new_coord in interp_span_x_coords:										# Go through all interpolated coordinates
			if new_coord > point_b and last_coord_b + 1 < len(span_x_coords):		# If the interpolated coord. is above point b...
				point_a, point_b = point_b, span_x_coords[last_coord_b + 1]			# ...set point a as previous point b, and b as next
				last_coord_b += 1													# Save the index of the new point b
			span_cut = []
			for i, chord_point in enumerate(sections_data[last_coord_b-1][1]):		# For each chord point of the section
				z = chord_point[0]													# Get the z-coordinate on the chord
				load = interp_two_points(new_coord, point_a, point_b, 
					chord_point[1], sections_data[last_coord_b][1][i][1])			# Interpolate the load !SPAN-WISE!, on that chord point
				span_cut.append([z, load])
			if n_chord is not None:													# Also apply interpolation chord-wise...
				interp_span_cut = []
				z = list(zip(*span_cut))[0]											# Extract the z-coordinates of evey chord-point in this section
				loads = list(zip(*span_cut))[1]										# Extract the load at evey chord-point in this section
				interp_z_coords = interp_coords(z, n_chord, N_z, C_a, "z")			# Interpolate the z-coordinates on n new locations
				point_c_a, point_c_b = z[0], z[1]									# Set the first points before and after the new inter. coord. 
				last_coord_c_b = 1													# Save the index of the point b
				for new_c_coord in interp_z_coords:									# Go trough all interpolated chord coordinates
					if new_c_coord < point_c_b and last_coord_c_b + 1 < len(z):		# If the interpolated coord. is below point b...
						point_c_a, point_c_b = point_c_b, z[last_coord_c_b + 1]		# ...set point a as previous point b, and b as next
						last_coord_c_b += 1											# Save the index of the new point b
					load_c = interp_two_points(new_c_coord, point_c_a, point_c_b, 	
					   loads[last_coord_c_b-1], loads[last_coord_c_b])				# Interpolate the load CHORD-WISE 
					interp_span_cut.append([new_c_coord, load_c])
				span_cut = interp_span_cut											# Save the section loads
			interp_sections_data.append([new_coord, span_cut])
		span_x_coords, sections_data = interp_span_x_coords, interp_sections_data
	# If specified, plot the load
	if do_plot:
		plot_data(sections_data)
	# Return the span coordinates (the x coordinate at every cut)
	# And the sections data, as follows:
	#	[
	#		[x1, [z1, load], [z2, load], [z3, load], ...],
	#		[x2, [z1, load], [z2, load], ...],
	#		...
	#	]
	return span_x_coords, sections_data

def plot_data(sections):
	"""
	Plot the resulting data in 3D
	"""
	fig = plt.figure()
	# Extract the x, z, and load from the multidimensional data array
	x, z, loads = [], [], []
	for section in sections:
		x_section = section[0]
		for chord_point in section[1]:
			x.append(x_section)
			z.append(chord_point[0])
			loads.append(chord_point[1])
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, z, loads, marker="x")
	ax.set_xlabel('X (span) [m]')
	ax.set_ylabel('Z (chord) [m]')
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
	if point_a == point_b:	# Make sure point a is different to point b (no zero division)
		return load_a
	# The load at the new coordinate is interpolated by fitting a line to the known
	# load before and after that new coordinate
	return load_a + (new_coord - point_a) * (load_b - load_a) / (point_b - point_a)

## Get load by opening a file
#get_load(n_span=100, n_chord=150, do_plot=True)
## Get load with fixed load at 25% chord
#get_load(C_a=0.605, l_a=2.661, n_span=100, do_plot=True, fixed_load=5540)