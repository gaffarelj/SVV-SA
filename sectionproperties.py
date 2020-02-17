import math
import numpy as np

Ha = 0.173
Ca = 0.484
tskin = 0.0011
beta = math.radians(12.5687)

#stiffener positions (origin at hinge line, from bottom left going clockwise, in meters)
#[z, y]
zy1 = [-0.321833, -0.016870]
zy2 = [-0.246166, -0.033740]
zy3 = [-0.170499, -0.050611]
zy4 = [-0.094832, -0.067481]
zy5 = [-0.019165, -0.084351]
zy6 = [0.054024, -0.067555]
zy7 = [0.0865, 0]
zy8 = [0.054024, 0.067555]
zy9 = [-0.019165, 0.084351]
zy10 = [-0.094832, 0.067481]
zy11 = [-0.170499, 0.050611]
zy12 = [-0.246166, 0.033740]
zy13 = [-0.321833, 0.016870]
boomcoords = np.array([zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11, zy12, zy13])

#stiffener area (in meters^2)
stiff_area = 3/78125

#enclosed area (in meters^2)
Am = math.pi*(Ha/2)**2+0.3975*0.0865

#centroid (z-coordinate from hinge line)
sumA_skin = (0.5*math.pi*Ha*tskin)+2*(0.4068*tskin)
sumAz_skin = (0.5*math.pi*Ha*tskin)*(Ha/math.pi)+2*(0.4068*tskin)*(-0.2034*math.cos(beta))
sumz_boom = 0
for i in range(13):
    sumz_boom += boomcoords[i,0]
sumAz_boom = sumz_boom*stiff_area
z_centroid = (sumAz_skin + sumAz_boom)/(sumA_skin + 13*stiff_area)

#moments of inertia (Iyy and Izz about the hinge line, in m^4)
#Iyy
Iyy_thin = ((tskin*0.4068**3*(math.cos(beta))**2)/12)+tskin*0.4068*(0.2034*math.cos(beta))**2




