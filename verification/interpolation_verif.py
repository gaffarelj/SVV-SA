import sys
import os.path
from matplotlib import pyplot as plt
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import interpolation as I

# Unit tests

print("interp_two_points")

unittest_interp = I.interp_two_points(1., 0., 2., 50., 100.)
print(75, unittest_interp)

unittest_interp = I.interp_two_points(-50, -100, 100., 100., -100.)
print(50, unittest_interp)


print("comp_theta")

unittest_theta = I.comp_theta(5, 5)
analytical_theta = 2.513274123
print(analytical_theta, unittest_theta, analytical_theta-unittest_theta)

unittest_theta = I.comp_theta(10, 1200)
analytical_theta = 0.0235619449
print(analytical_theta, unittest_theta, analytical_theta-unittest_theta)


print("comp_coord")

unittest_coord = I.comp_coord(10, 1200, 10, "z")
analytical_coord = -0.0015506
print(analytical_coord, unittest_coord, analytical_coord-unittest_coord)

unittest_coord = I.comp_coord(10, 500, 10, "x")
analytical_coord = 0.0089293
print(analytical_coord, unittest_coord, analytical_coord-unittest_coord)


print("interp_coords")

N, l, c_type = 5, 3, "z"
coords = [1, 2, 3, 4, 5]
coords = [I.comp_coord(i, N, l , c_type) for i in coords]
unittest_interp_coords = I.interp_coords(coords, 20, N, l, c_type)
plt.plot(coords, np.zeros(len(coords)), marker="x", linewidth=0.5, label="Given")
plt.plot(unittest_interp_coords, np.ones(len(unittest_interp_coords)), marker="x", linewidth=0.5, label="Interpolated")
plt.legend()
plt.xlabel("z [m]")
plt.yticks([], [])
plt.show()

N, l, c_type = 5, 10, "x"
coords = [1, 2, 3, 4, 5]
coords = [I.comp_coord(i, N, l , c_type) for i in coords]
unittest_interp_coords = I.interp_coords(coords, 15, N, l, c_type)
plt.plot(coords, np.zeros(len(coords)), marker="x", linewidth=0.5)
plt.plot(unittest_interp_coords, np.ones(len(unittest_interp_coords)), marker="x", linewidth=0.5)
plt.show()


# Subsystem test

I.get_load(n_span=None, n_chord=None, do_plot=True)
I.get_load(n_span=100, n_chord=150, do_plot=True)