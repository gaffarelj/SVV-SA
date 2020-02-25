import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

class stress():
	def __init__(self, Mz, My, sect, q1, q2, q3, q4, q5, q6):
		self.Iyy = sect.Iyy
		self.Izz = sect.Izz
		self.Mz = Mz
		self.My = My
		self.stresses = []
		self.shear_flows = []
		self.r = sect.r
		self.beta = sect.beta
		self.l_topskin = sect.l_topskin
		self.l_spar_to_end = sect.l_spar_to_end
		self.zc = sect.z_centroid
		self.q1, self.q2, self.q3, self.q4, self.q5, self.q6 = q1, q2, q3, q4, q5, q6

	def section_stress(self):
		# Hinge
		z = 0
		ys = np.arange(-self.r, self.r, 0.0025)
		for y in ys:
			self.add_stress(z, y)
			if y > 0:
				s = y
				self.add_shear(z, y, s, self.q2)
			else:
				s = -y
				self.add_shear(z, y, s, self.q5)
		# Circle
		for theta in np.arange(0, np.pi, 0.025):
			z = np.sin(theta) * self.r
			y = -np.cos(theta) * self.r
			self.add_stress(z, y)
			s = theta - np.pi/2
			if theta < np.pi/2:
				self.add_shear(z, y, s, self.q6)
			else:
				self.add_shear(z, y, s, self.q1)
		# Top sheet
		for s in np.arange(0, self.l_topskin, 0.005):
			z = -s * np.cos(self.beta)
			y = self.r - s * np.sin(self.beta)
			self.add_stress(z, y)
			self.add_shear(z, y, s, self.q3)
		# Bottom sheet
		for s in np.arange(0, self.l_topskin, 0.005):
			z = -s * np.cos(self.beta)
			y = -self.r + s * np.sin(self.beta)
			self.add_stress(z, y)
			self.add_shear(-z - self.l_spar_to_end, -y - self.r, s, self.q4)

	def add_stress(self, z, y):
		z -= self.zc
		strs = self.stress_norm_xx(y, z)
		self.stresses.append([z, y, strs])

	def add_shear(self, z, y, s, q):
		z -= self.zc
		shear = q(s)
		self.shear_flows.append([z, y, shear])

	def stress_norm_xx(self, y, z):
		return self.Mz / self.Izz * y + self.My / self.Iyy * z
	
	def plot_everything(self):
		data = np.array(self.stresses)
		z, y, s = data[:,0], data[:,1], data[:,2]
		plt.xlim(0.225, -0.325)
		plt.ylim(-0.09, 0.09)
		plt.xlabel("z [m]")
		plt.ylabel("y [m]")
		plt.scatter(z, y, c=s)
		plt.colorbar()
		plt.set_cmap("jet")
		#z_line = np.arange(-0.2, 0.2, 0.001)
		#alpha = np.arctan(-self.My*self.Izz/self.Mz/self.Iyy)
		#y_line = [z_c * np.sin(alpha) for z_c in z_line]
		#plt.plot(z_line, y_line)
		plt.show()	

	def plot_stress(self):
		data = np.array(self.stresses)
		z, y, s = data[:,0], data[:,1], data[:,2]
		plt.xlim(0.225, -0.325)
		plt.ylim(-0.09, 0.09)
		plt.xlabel("z [m]")
		plt.ylabel("y [m]")
		plt.scatter(z, y, c=s)
		plt.colorbar()
		plt.set_cmap("jet")
		#z_line = np.arange(-0.2, 0.2, 0.001)
		#alpha = np.arctan(-self.My*self.Izz/self.Mz/self.Iyy)
		#y_line = [z_c * np.sin(alpha) for z_c in z_line]
		#plt.plot(z_line, y_line)
		plt.show()
			
	def plot_shear_flows(self):
		data = np.array(self.shear_flows)
		z, y, s = data[:,0], data[:,1], data[:,2]
		plt.xlim(0.225, -0.325)
		plt.ylim(-0.09, 0.09)
		plt.xlabel("z [m]")
		plt.ylabel("y [m]")
		plt.scatter(z, y, c=s)
		plt.colorbar()
		plt.set_cmap("jet")
		plt.show()