import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

unittest_interp = interp_two_points(1., 0., 2., 50., 100.)        #testing interpolation function with point (0,50) and (2,100) and checking whether it gives 75
if unittest_interp != 75.0:
    raise AssertionError
    
unittest_interp = interp_two_points(-50, -100, 100., 100., -100.)   #testing interpolation function with positive and negative values and testing whether the correct analytically found answer is given
if unittest_interp != 50.:
    raise AssertionError

print("No AssertionErrors have arised, and the interp_two_points function works correctly")


def comp_theta(i, N):
	"""
	Compute angle theta, as specified in the Assignment PDF.
	"""
	return (i - 1) / N * math.pi

unittest_theta = comp_theta(5, -5)        #testing theta function with point i = 5 and N = -5 and checking whether it gives -2.51
if round(unittest_theta, 2) != -2.51:
    raise AssertionError
    
unittest_theta = comp_theta(10, 1200)        #testing theta function with point i = 5 and N = -5 and checking whether it gives -2.51
if round(unittest_theta, 2) != 0.02:
    raise AssertionError

print("No AssertionErrors have arised, and the comp_theta function works correctly")

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

unittest_coord = comp_coord(10, 1200, 10, "z")        #testing coordinate function with point i = 10 and N = 1200, l = 10 and c_type = z and checking whether it gives -0.0016
if round(unittest_coord, 4) != -0.0016:
    raise AssertionError
    
unittest_coord = comp_coord(10, 1200, 10, 5)        #testing coordinate function with point i = 10 and N = 1200, l = 10 and c_type = else and checking whether it gives 0.0016
if round(unittest_coord, 4) != 0.0016:
    raise AssertionError

print("No AssertionErrors have arised, and the comp_coord function works correctly")




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


unittest_interp_coords = interp_coords([1,2], 1, 1, 1, "z")        
    raise AssertionError
    
unittest_interp_coords = interp_coords([1,2], 3, 3, 3, "x")        
if unittest_interp_coords != [1, 0.3749999999999999, 2]:
    raise AssertionError

print("No AssertionErrors have arised, and the comp_coord function works correctly")
