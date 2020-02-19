from integration import integrate
import math
import sectionproperties as sp

Am, z_centroid, Iyy, Izz, tskin, boomcoords, boomcoords_hinge, stiff_area = sp.get_geometry()
Ha = 0.173
Ca = 0.484
tspar = 0.0025
hstiff = 0.014
tstiff = 0.0012
wstiff = 0.018
beta = math.atan(8.65 / 39.75)
r = Ha / 2


# SECTION 1
def qb_1(theta):  # 0 to pi/2
    c = -1 / Izz
    boom8_theta = math.asin(boomcoords[7, 1] / r)
    if theta < boom8_theta:
        qb = c * (tskin * r ** 2) * integrate(math.sin, 0, theta)
    else:
        qb = c * ((tskin * r ** 2) * integrate(math.sin, 0, theta) + stiff_area * boomcoords[7, 1])
    return qb


# SECTION 2
def qb_2(s2):  # 0 to 0.0865
    c = -1 / Izz
    qb = c * integrate(lambda x: tspar * x, 0, s2)
    return qb


# SECTION 3
def qb_3(s3):  # 0 to 0.4068
    c = -1 / Izz
    boom9_s3 = -boomcoords_hinge[8, 0] / math.cos(beta)
    boomspace_s3 = (boomcoords_hinge[8, 0] - boomcoords_hinge[9, 0]) / math.cos(beta)
    if s3 < boom9_s3:
        qb = c * tskin * integrate(lambda x: r - x, 0, s3) + qb_1(math.pi/2) + qb_2(0.0865)
    else:
        nbooms = int((s3 - boom9_s3) / boomspace_s3)
        qb_booms = 0
        for i in range(8, 8 + nbooms):
            qb_booms += stiff_area * boomcoords_hinge[i, 1]
        qb_skin = tskin * integrate(lambda x: r - x, 0, s3)
        qb = c * (qb_skin + qb_booms) + qb_1(math.pi/2) + qb_2(0.0865)
    return qb


# SECTION 4
def qb_4(s4):  # 0 to 0.4068
    c = -1 / Izz
    boom1_s4 = (Ca-r+boomcoords_hinge[0, 0]) / math.cos(beta)
    boomspace_s4 = (boomcoords_hinge[0, 0] - boomcoords_hinge[1, 0]) / math.cos(beta)
    if s4 < boom1_s4:
        qb = c * integrate(lambda x: -tskin * math.sin(beta) * x, 0, s4) + qb_3(0.4068)
    else:
        nbooms = int((s4 - boom1_s4) / boomspace_s4)
        qb_booms = 0
        for i in range(0, 0 + nbooms):
            qb_booms += stiff_area * boomcoords_hinge[i, 1]
        qb_skin = integrate(lambda x: -tskin * math.sin(beta) * x, 0, s4)
        qb = c * (qb_skin + qb_booms) + qb_3(0.4068)
    return qb


# SECTION 5
def qb_5(s5):  # 0 to 0.0865
    c = -1 / Izz
    qb = c * integrate(lambda x: -tspar * x, 0, s5)
    return qb


# SECTION 6
def qb_6(theta):
    c = -1 / Izz
    boom6_theta = math.asin(boomcoords[5, 1] / r)
    if theta < boom6_theta:
        qb = c * (tskin * r ** 2) * integrate(math.sin, theta, 0) + qb_5(0.0865) + qb_4(0.4068)
    else:
        qb = c * ((tskin * r ** 2) * integrate(math.sin, theta, 0) + stiff_area * boomcoords[5, 1]) + qb_5(0.0865) + qb_4(0.4068)
    return qb
