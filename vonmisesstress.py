import math

def vmstress(shearstress,mz,my):
	"""
	Calculate the Von Mises stress with formula which can be simplified as:
	sigma_vm = sqrt(sigma_xx^2 + 3*tau_yz^2)
	Inputs:
	- shearstress:
	- mz:
	- my:
	"""