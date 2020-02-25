import sectionproperties as SP
import shearcentre as SC
import torsionstiffness as TS
import macaulay_new as MA

section_prop = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
# The following parameters can be used out of section_prop (example: print(section_prop.z_centroid)
# Iyy, Izz, Am, z_centroid, stiff_area, boomcoords, boomcoords_hinge,
# tskin, tspar, hstiff, tstiff, wstiff, Ha, Ca, r, l_topskin, beta, l_spar_to_end
SC.set_sect(section_prop)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
print("shear center:", xi)

x = 0.5
def vmstress(x=0.5,y,z):
	"""
	Calculate the Von Mises stress with formula which can be simplified as:
	sigma_vm = sqrt(sigma_xx^2 + 3*tau_yz^2)
	Inputs:

	"""
	sigma =MA.My(x)/section_prop.Iyy*z + 
	vm = math.sqrt(sigma_xx**2+3*tau_yz**2)
	return vm