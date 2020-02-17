import math
import numpy as np

Ha = 0.173
Ca = 0.484
tskin = 0.0011

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
boomcoord = np.array(zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11, zy12, zy13)

#stiffener area (in meters^2)
stiff_area = 3/78125

#centroid
sumAzskin = (0.5*math.pi*Ha*tskin)*(Ha/math.pi)+2*(0.4068*tskin)*(-0.2034*math.cos(12.5687))
sumAzboom = 0
