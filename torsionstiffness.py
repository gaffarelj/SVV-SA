import math
# G*dtheta/dz = 1/2/A1 * f qds/t
# j = T/(G*dtheta/dz)

def spacinggenerator(distance, ha):
	skin_spar = math.pi*ha/4 - 2*distance
	result = [distance, distance, skin_spar, ha, skin_spar, distance, distance]
	print(result)
	return result

def roundintergration(sf_ab, sf_bnb, distance, thickness): # shearflow after boom and before next boom for each section
	return (sf_ab+sf_bnb)*distance/2/thickness

def torsionstiffness(T,shearflows,t_skin,t_spar,distance_between_booms,ha):
	# 7 regions to calculate for left cell
	spacing = spacinggenerator(distance_between_booms,ha)
	a1 = math.pi*ha*ha/8
	intergrations = []

	for i in range(0,3): # semi circle above
		sfab, sfbnb = shearflows[i]
		distance = spacing[i]
		r = roundintergration(sfab, sfbnb, distance,t_skin)
		intergrations.append(r)
	a,b = shearflows[3]
	r_spar = roundintergration(a,b,spacing[3],t_spar)
	intergrations.append(r_spar)
	for i in range(4,len(shearflows)):
		sfab, sfbnb = shearflows[i]
		distance = spacing[i]
		r = roundintergration(sfab, sfbnb, distance, t_skin)
		intergrations.append(r)

	print(intergrations)
	G_dtheta_dz = sum(intergrations)/2/a1

	J = T/G_dtheta_dz
	return J


T = 1  # torsion
ha = 2  # height of aileron
distance = 0.07752484  # m
t_skin = 0.0011  # thickness
t_spar = 0.0025

shearflows =[(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
print(len(shearflows))

a = torsionstiffness(T, shearflows, t_skin, t_spar, distance, ha)
print(a)