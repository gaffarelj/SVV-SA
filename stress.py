import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import math

class stress():
	def __init__(self, Mz, My, Sz, Sy, T, sect, q1, q2, q3, q4, q5, q6, show_plot=True):
		self.Iyy = sect.Iyy
		self.Izz = sect.Izz
		self.Mz = Mz
		self.My = My
		self.Sz = Sz
		self.Sy = Sy
		self.T = T
		self.stresses = []
		self.shear_flows = []
		self.shear_stress = []
		self.vm_stresses = []
		self.r = sect.r
		self.beta = sect.beta
		self.l_topskin = sect.l_topskin
		self.l_spar_to_end = sect.l_spar_to_end
		self.zc = sect.z_centroid
		self.q1, self.q2, self.q3, self.q4, self.q5, self.q6 = q1, q2, q3, q4, q5, q6
		self.tskin = sect.tskin
		self.tspar = sect.tspar
		self.show_plot = show_plot

	def section_stress(self):
		# Hinge
		z = 0
		ys = np.arange(-self.r, self.r, 0.0025)
		for y in ys:
			sigma = self.add_stress(z, y)
			if y > 0:
				s = y
				tau = self.add_shear(z, y, s, self.q2, self.Sz, self.Sy, self.T)
			else:
				s = -y
				tau = self.add_shear(z, y, s, self.q5, self.Sz, self.Sy, self.T)
			self.add_vonMises(y, z, sigma, tau, self.tspar)
		# Circle
		for theta in np.arange(0, np.pi, 0.025):
			z = np.sin(theta) * self.r
			y = -np.cos(theta) * self.r
			sigma = self.add_stress(z, y)
			s = theta - np.pi / 2
			if theta < np.pi / 2:
				tau = self.add_shear(z, y, s, self.q6, self.Sz, self.Sy, self.T)
			else:
				tau = self.add_shear(z, y, s, self.q1, self.Sz, self.Sy, self.T)
			self.add_vonMises(y, z, sigma, tau, self.tskin)
		# Top sheet
		for s in np.arange(0, self.l_topskin, 0.0025):
			z = -s * np.cos(self.beta)
			y = self.r - s * np.sin(self.beta)
			sigma = self.add_stress(z, y)
			tau = self.add_shear(z, y, s, self.q3, self.Sz, self.Sy, self.T)
			self.add_vonMises(y, z, sigma, tau, self.tskin)
		# Bottom sheet
		for s in np.arange(0, self.l_topskin, 0.0025):
			z = -s * np.cos(self.beta)
			y = -self.r + s * np.sin(self.beta)
			sigma = self.add_stress(z, y)
			tau = self.add_shear(-z - self.l_spar_to_end, -y - self.r, s, self.q4, self.Sz, self.Sy, self.T)
			self.add_vonMises(y, z, sigma, tau, self.tskin)

	def add_stress(self, z, y):
		z -= self.zc
		strs = self.stress_norm_xx(y, z)
		self.stresses.append([z, y, strs])
		return strs

	def add_shear(self, z, y, s, q, Sz, Sy, T):
		z -= self.zc
		shear = q(s, Sz, Sy, T)
		self.shear_flows.append([z, y, shear])
		return shear

	def stress_norm_xx(self, y, z):
		return self.Mz / self.Izz * y + self.My / self.Iyy * z

	def add_vonMises(self, y, z, stress, q, t):
		z -= self.zc
		shear_stress = q / t
		vm = math.sqrt(stress ** 2  + 3 * shear_stress ** 2)
		self.vm_stresses.append([z, y, vm])
		self.shear_stress.append([z, y, shear_stress])

	def plot(self, data, fname, label, x):
		z, y, s = data[:,0], data[:,1], data[:,2]
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_aspect('equal')
		plt.xlim(0.225, -0.325)
		plt.ylim(-0.125, 0.125)
		plt.xlabel("z [m]")
		plt.ylabel("y [m]")
		plt.scatter(z, y, c=s, s=3)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="3%", pad=0.1)
		plt.colorbar(cax=cax, label=label)
		plt.set_cmap("jet")
		plt.savefig(f"plots/stresses/{fname}-{str(x).replace('.', '-')}.pdf", bbox_inches='tight')
		#z_line = np.arange(-0.2, 0.2, 0.001)
		#alpha = np.arctan(-self.My*self.Izz/self.Mz/self.Iyy)
		#y_line = [z_c * np.sin(alpha) for z_c in z_line]
		#plt.plot(z_line, y_line)
		if self.show_plot:
			plt.show()
		else:
			plt.close()

	def plot_all(self, x):
		self.plot(np.array(self.stresses), "s", "Normal stress [Pa]", x)
		self.plot(np.array(self.vm_stresses), "vm", "von Mises stress [Pa]", x)
		self.plot(np.array(self.shear_flows), "sf", "Shear flow [N/m]", x)