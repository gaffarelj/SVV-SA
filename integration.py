import numpy as np

def integrate(f, xs, a, b, threshold=1e-6):
	# TO CHANGE
	n = 1e4
	y = [f(x) for x in xs]
	return np.trapz(y, xs)
	# TO CHANGE


integ = integrate(np.sin, np.arange(0, np.pi/2, 1/1e4), 0, np.pi/2)
print(integ)