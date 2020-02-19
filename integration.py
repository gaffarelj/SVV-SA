def integrate(f, a, b, dx=1e-4):
	n = int((b - a) / dx)
	area = 0
	for i in range(1, n):
		x = a + i * dx
		area += f(x)
	return dx * (area + (f(a) + f(b)) / 2)