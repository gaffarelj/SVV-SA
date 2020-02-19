import numpy as np


def integrate(f, a, b, threshold=1e-6):
	# TO CHANGE
	n = 1e4
	xs = np.arange(a, b, 1/n)
	y = [f(x) for x in xs]
	return np.trapz(y, xs)
	# TO CHANGE
	
def integrate_b(f, a, b, threshold=1e-6):
	pass
	
integ = integrate(np.sin, 0, np.pi/2)
integ_b = integrate_b(np.sin, 0, np.pi/2)
print(integ, integ_b)