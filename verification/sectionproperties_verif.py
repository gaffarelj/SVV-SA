import math
import numpy as np


def get_boomcoords():
    # stiffener positions (origin at hinge line, from bottom right going clockwise, in meters)
    # [z, y]
    zy1 = [-6.75, 1.25]
    zy2 = [2.5, 0.]
    zy3 = [-6.75, -1.25]
    return np.array([zy1, zy2, zy3])



Ha = 2.5
tskin = 0.2
tspar = 0.2
stiff_area = 0.5
boomcoords = get_boomcoords()
beta=math.atan(2.5/12.5)


def comp_centroid(Ha, tskin, tspar, boomcoords, stiff_area):
    # centroid (z-coordinate from hinge line, in meters)
    sumA_skin = (0.5 * math.pi * Ha * tskin) + 2 * (12.74 * tskin) + (tspar * Ha)
    sumAz_skin = (0.5 * math.pi * Ha * tskin) * (Ha / math.pi) + 2 * (12.74 * tskin) * (-6.25)
    sumAz_boom = 0
    for i in range(len(boomcoords)):
        sumAz_boom += boomcoords[i,0] * stiff_area
    z_centroid = (sumAz_skin + sumAz_boom) / (sumA_skin + 3 * stiff_area)
    boomcoords_centroid = np.zeros((3, 2))
    for i in range(len(boomcoords)):        # shift boom origin from hinge line to centroid
        boomcoords_centroid[i, 0] = boomcoords[i, 0] - z_centroid
        boomcoords_centroid[i, 1] = boomcoords[i, 1]
    if round(z_centroid, 2) != -4.66:         #Here the centroid coordinate is rounded to 2 digits. 
        raise AssertionError                    #If the z_centroid coordinate is not equal to algebraically calculated value, the assertionerror is raised
    return z_centroid, boomcoords_centroid

def comp_Iyy(tskin, beta, z_centroid, tspar, Ha, boomcoords, stiff_area):
    # Iyy
    Iyy_thin = ((tskin * (12.74 ** 3) * ((math.cos(beta)) ** 2)) / 12) + tskin * 12.74 * (6.25 + z_centroid) ** 2
    Iyy_spar = (Ha * (tspar ** 3)) / 12 + tspar * Ha * z_centroid ** 2
    Iyy_ring = tskin * (Ha / 2) ** 3 * ((math.pi ** 2 - 8) / (2 * math.pi)) + (0.5 * math.pi * Ha * tskin) * ((Ha / math.pi) - z_centroid) ** 2
    Iyy_booms = 0
    for i in range(len(boomcoords)):
        Iyy_booms += ((boomcoords[i, 0]) ** 2) * stiff_area
    Iyy_total = 2 * Iyy_thin + Iyy_spar + Iyy_ring + Iyy_booms
    if round(Iyy_total, 3) != 307.997:         #Here the Iyy is rounded to 3 digits. 
        print(Iyy_total)
        raise AssertionError                    #If the Iyy_total is not equal to algebraically calculated value, the assertionerror is raised
    return Iyy_total

def comp_Izz(tskin, beta, tspar, Ha, boomcoords, stiff_area):
    # Izz
    Izz_thin = ((tskin * (0.4068 ** 3) * ((math.sin(beta)) ** 2)) / 12) + tskin * 0.4068 * (0.19875 * math.sin(beta)) ** 2
    Izz_spar = (tspar * (Ha ** 3)) / 12 # This is correct
    Izz_ring = ((Ha / 2) ** 3 * tskin * math.pi) / 2    # This is correct
    #print(Izz_thin, Izz_spar, Izz_ring)
    Izz_booms = 0
    for i in range(len(boomcoords)):
        Izz_booms += (boomcoords[i, 1]) ** 2 * stiff_area
    Izz_total = 2 * Izz_thin + Izz_spar + Izz_ring + Izz_booms
    if round(Izz_total, 3) != 2.437:         #Here the Izz is rounded to 3 digits. 
        print(Izz_total)
        raise AssertionError                    #If the Izz_total is not equal to algebraically calculated value, the assertionerror is raised
    return Izz_total


z_centroid, boomcoords_centroid = comp_centroid(Ha, tskin, tspar, boomcoords, stiff_area)
print("No AssertionErrors have arised, and the comp_centroid function works correctly")
comp_Iyy(tskin, beta, z_centroid, tspar, Ha, boomcoords, stiff_area)
print("No AssertionErrors have arised, and the comp_Iyy function works correctly")
comp_Izz(tskin, beta, tspar, Ha, boomcoords, stiff_area)
print("No AssertionErrors have arised, and the comp_Izz function works correctly")
