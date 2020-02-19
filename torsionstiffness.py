import math
from sectionproperties import tskin,tspar,Ha
def spacinggenerator(distance, ha):
	"""
	Return the list of length of 5 parts in left cell
	Input:
	- distance: the distance calculated between evenly distributed booms
	- ha: height of the aileron
	"""
	skin_spar = math.pi*ha/4 - distance # the tiny part on semi circle on the left to spar
	result = [distance, skin_spar, ha, skin_spar, distance]
	print(result)
	return result

def roundintergration(sf_ab, sf_bnb, distance, thickness):
	"""
	Round integration between booms
	Inputs:
	- sf_ab: shear flow after the boom
	- sf_bnb: shear flow before next boom
	"""
	return (sf_ab+sf_bnb)*distance/2/thickness

def torsionstiffness(T, shearflows, t_skin, t_spar, distance_between_booms, ha):
	"""
	Calculate the torsional stiffness using shear flow in left cell
	Formulas:
	G*dtheta/dz = 1/2/A1 * f qds/t
	J = T/(G*dtheta/dz)

	Inputs:
	- T: Torsion
	- shearflows: shear_floe per section in left cell
	- t_skin: skin thickness
	- t_spar: spar thickness
	- distance_between_booms: evenly distributed distance between booms
	- ha: height of the aileron
	"""
	# 5 regions to calculate for left cell
	spacing = spacinggenerator(distance_between_booms,ha)
	a1 = math.pi*ha**2/8  # enclosed area
	intergrations = []

	for i in [0, 1, 3, 4]:  # semi circle
		sfab, sfbnb = shearflows[i]
		distance = spacing[i]
		r = roundintergration(sfab, sfbnb, distance,t_skin)
		intergrations.append(r)

	a, b = shearflows[2]  # spar
	r_spar = roundintergration(a, b, spacing[2], t_spar)
	intergrations.append(r_spar)

	print(intergrations)

	G_dtheta_dz = sum(intergrations)/2/a1
	J = T/G_dtheta_dz
	return J


T = 1  # torsion
ha = Ha  # height of aileron, m
distance = 0.07752484  # distance between booms, m
t_skin = tskin
t_spar = tspar

shearflows_leftcell = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
print(len(shearflows_leftcell))

a = torsionstiffness(T, shearflows_leftcell, t_skin, t_spar, distance, ha)
print(a)