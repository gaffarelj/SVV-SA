import math
import numpy as np
"""
Section Properties - idealized aileron

Computes the centroid, enclosed area and moments of inertia (Iyy and Izz) for the idealized section.
All units are meters unless otherwise specified.
Stiffener positions are hardcoded (frame of reference is the hinge line), more can be added.
Inputs are below, as given in Aircraft Data, coded for CRJ700.
"""

def get_geometry(Nstiffeners=13, Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025,
                hstiff=0.014, tstiff=0.0012, wstiff=0.018, beta=math.atan(8.65/39.75)):
    boomcoords = get_boomcoords()
    boomcoords_hinge = boomcoords
    # stiffener area (in meters^2)
    stiff_area = (hstiff + wstiff) * tstiff     # small angle approximation is used
    # enclosed area (in meters^2)
    Am = math.pi * (Ha / 2) ** 2 + 0.3975 * 0.0865
    z_centroid, boomcoords = comp_centroid(Ha, tskin, tspar, boomcoords, stiff_area)
    # moments of inertia (Iyy and Izz about the centroid, in meters^4)
    Iyy = comp_Iyy(tskin, beta, z_centroid, tspar, Ha, boomcoords, stiff_area)
    Izz = comp_Izz(tskin, beta, tspar, Ha, boomcoords, stiff_area)
    return Am, z_centroid, Iyy, Izz, tskin, boomcoords, boomcoords_hinge, stiff_area

def get_boomcoords():
    # stiffener positions (origin at hinge line, from bottom right going clockwise, in meters)
    # [z, y]
    zy1 = [-0.35671, -0.00888]
    zy2 = [-0.27513, -0.02663]
    zy3 = [-0.19355, -0.04438]
    zy4 = [-0.11197, -0.06213]
    zy5 = [-0.03039, -0.07989]
    zy6 = [0.04924, -0.07112]
    zy7 = [0.0865, 0]
    zy8 = [0.04924, 0.07112]
    zy9 = [-0.03039, 0.07989]
    zy10 = [-0.11197, 0.06213]
    zy11 = [-0.19355, 0.0444]
    zy12 = [-0.27513, 0.02663]
    zy13 = [-0.35671, 0.00888]
    return np.array([zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11, zy12, zy13])

def comp_centroid(Ha, tskin, tspar, boomcoords, stiff_area):
    # centroid (z-coordinate from hinge line, in meters)
    sumA_skin = (0.5 * math.pi * Ha * tskin) + 2 * (0.4068 * tskin) + (tspar * Ha)
    sumAz_skin = (0.5 * math.pi * Ha * tskin) * (Ha / math.pi) + 2 * (0.4068 * tskin) * (-0.19875)
    sumAz_boom = 0
    for i in range(len(boomcoords)):
        sumAz_boom += boomcoords[i, 0] * stiff_area
    z_centroid = (sumAz_skin + sumAz_boom) / (sumA_skin + 13 * stiff_area)
    boomcoords_centroid = np.zeros((13, 2))
    for i in range(len(boomcoords)):        # shift boom origin from hinge line to centroid
        boomcoords_centroid[i, 0] = boomcoords[i, 0] - z_centroid
        boomcoords_centroid[i, 1] = boomcoords[i, 1]
    return z_centroid, boomcoords_centroid

def comp_Iyy(tskin, beta, z_centroid, tspar, Ha, boomcoords, stiff_area):
    # Iyy
    Iyy_thin = ((tskin * (0.4068 ** 3) * ((math.cos(beta)) ** 2)) / 12) + tskin * 0.4068 * (0.19875 + z_centroid) ** 2
    Iyy_spar = (Ha * (tspar ** 3)) / 12 + tspar * Ha * z_centroid ** 2
    Iyy_ring = tskin * (Ha / 2) ** 3 * ((math.pi ** 2 - 8) / (2 * math.pi)) + (0.5 * math.pi * Ha * tskin) * ((Ha / math.pi) - z_centroid) ** 2
    Iyy_booms = 0
    for i in range(len(boomcoords)):
        Iyy_booms += ((boomcoords[i, 0]) ** 2) * stiff_area
    Iyy_total = 2 * Iyy_thin + Iyy_spar + Iyy_ring + Iyy_booms
    return Iyy_total

def comp_Izz(tskin, beta, tspar, Ha, boomcoords, stiff_area):
    # Izz
    Izz_thin = ((tskin * (0.4068 ** 3) * ((math.sin(beta)) ** 2)) / 12) + tskin * 0.4068 * (0.19875 * math.sin(beta)) ** 2
    Izz_spar = (tspar * (Ha ** 3)) / 12 # This is correct
    Izz_ring = ((Ha / 2) ** 3 * tskin * math.pi) / 2    # This is correct
    print(Izz_thin, Izz_spar, Izz_ring)
    Izz_booms = 0
    for i in range(len(boomcoords)):
        Izz_booms += (boomcoords[i, 1]) ** 2 * stiff_area
    Izz_total = 2 * Izz_thin + Izz_spar + Izz_ring + Izz_booms
    return Izz_total

Am, z_centroid, Iyy, Izz, tskin, boomcoords, boomcoords_hinge, stiff_area = get_geometry()
#print(Am, z_centroid, Iyy, Izz, tskin, stiff_area, boomcoords, "\n\n", boomcoords_hinge, stiff_area)
print(Izz)