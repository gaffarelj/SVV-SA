import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

class stress():
	def __init__(self, Mz, My, sect):
		self.Iyy = sect.Iyy
		self.Izz = sect.Izz
		self.Mz = Mz
		self.My = My
		self.stresses = []
		self.r = sect.r
		self.beta = sect.beta
		self.l_topskin = sect.l_topskin
		self.zc = sect.z_centroid

	def section_stress(self):
		# Hinge
		z = 0
		ys = np.arange(-self.r, self.r, 0.005)
		for y in ys:
			self.add_res(z, y)
		# Circle
		for theta in np.arange(0, np.pi, 0.05):
			z = np.sin(theta) * self.r
			y = np.cos(theta) * self.r
			self.add_res(z, y)
		# Top sheet
		for p in np.arange(0, self.l_topskin, 0.005):
			z = -p * np.cos(self.beta)
			y = self.r - p * np.sin(self.beta)
			self.add_res(z, y)
		# Bottom sheet
		for p in np.arange(0, self.l_topskin, 0.005):
			z = -p * np.cos(self.beta)
			y = -self.r + p * np.sin(self.beta)
			self.add_res(z, y)

	def add_res(self, z, y):
		z -= self.zc
		strs = self.stress_norm_xx(y, z)
		self.stresses.append([z, y, strs])

	def stress_norm_xx(self, y, z):
		return self.Mz / self.Izz * y + self.My / self.Iyy * z

	def plot_stress(self):
		data = np.array(self.stresses)
		z, y, s = data[:,0], data[:,1], data[:,2]
		plt.xlim(0.225, -0.325)
		plt.ylim(-0.09, 0.09)
		plt.xlabel("z [m]")
		plt.ylabel("y [m]")
		plt.scatter(z, y, c=s)
		plt.colorbar()
		plt.show()