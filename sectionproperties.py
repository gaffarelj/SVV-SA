import math
import numpy as np
import matplotlib.pyplot as plt

Ha = 0.173
Ca = 0.484
tskin = 0.0011
beta = math.radians(12.5687)

# stiffener positions (origin at hinge line, from bottom left going clockwise, in meters)
# [z, y]
zy1 = [-0.35671025, -0.00887626]
zy2 = [-0.27513076, -0.02662878]
zy3 = [-0.19355126, -0.0443813]
zy4 = [-0.11197177, -0.06213382]
zy5 = [-0.03039227, -0.07988634]
zy6 = [0.04924123, -0.07111646]
zy7 = [0.0865, 0]
zy8 = [0.04924123, 0.07111646]
zy9 = [-0.03039227, 0.07988634]
zy10 = [-0.11197177, 0.06213382]
zy11 = [-0.19355126, 0.0443813]
zy12 = [-0.27513076, 0.02662878]
zy13 = [-0.35671025, 0.00887626]
boomcoords = np.array([zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11, zy12, zy13])

# stiffener area (in meters^2)
stiff_area = 3 / 78125

# enclosed area (in meters^2)
Am = math.pi * (Ha / 2) ** 2 + 0.3975 * 0.0865

# centroid (z-coordinate from hinge line, in meters)
sumA_skin = (0.5 * math.pi * Ha * tskin) + 2 * (0.4068 * tskin)
sumAz_skin = (0.5 * math.pi * Ha * tskin) * (Ha / math.pi) + 2 * (0.4068 * tskin) * (-0.2034 * math.cos(beta))
sumz_boom = 0
for i in range(13):
    sumz_boom += boomcoords[i, 0]
sumAz_boom = sumz_boom * stiff_area
z_centroid = (sumAz_skin + sumAz_boom) / (sumA_skin + 13 * stiff_area)
boomcoords_centroid = np.zeros((13, 2))
for i in range(13):
    boomcoords_centroid[i, 0] = boomcoords[i, 0] - z_centroid
    boomcoords_centroid[i, 1] = boomcoords[i, 1]

# moments of inertia (Iyy and Izz about the centroid, in meters^4)
# Iyy
Iyy_thin = ((tskin * 0.4068 ** 3 * (math.cos(beta)) ** 2) / 12) + tskin * 0.4068 * (
        0.2034 * math.cos(beta) + z_centroid) ** 2
Iyy_spar = (Ha * tskin ** 3) / 12
Iyy_ring = tskin * (Ha / 2) ** 3 * ((math.pi ** 2 - 8) / (2 * math.pi)) + (0.5 * math.pi * Ha * tskin) * (
        (Ha / math.pi) - z_centroid) ** 2
Iyy_booms = 0
for i in range(13):
    Iyy_booms += (boomcoords_centroid[i, 0]) ** 2 * stiff_area
Iyy_total = Iyy_thin + Iyy_spar + Iyy_ring + Iyy_booms

# Izz
Izz_thin = ((tskin * 0.4068 ** 3 * (math.sin(beta)) ** 2) / 12) + tskin * 0.4068 * (0.2034 * math.sin(beta)) ** 2
Izz_spar = (tskin * Ha ** 3) / 12
Izz_ring = ((Ha / 2) ** 3 * tskin * math.pi) / 2
Izz_booms = 0
for i in range(13):
    Izz_booms += (boomcoords_centroid[i, 1]) ** 2 * stiff_area
Izz_total = Izz_thin + Izz_spar + Izz_ring + Izz_booms
