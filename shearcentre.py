from integration import integrate
import math
import sectionproperties as sp
import numpy as np

Am, z_centroid, Iyy, Izz, tskin, boomcoords, boomcoords_hinge, stiff_area = sp.get_geometry()
Ha = 0.173
Ca = 0.484
tspar = 0.0025
hstiff = 0.014
tstiff = 0.0012
wstiff = 0.018
beta = math.atan(8.65 / 39.75)
l_topskin = 0.3975 / math.cos(beta)
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
        qb = c * tskin * integrate(lambda x: r - x * math.sin(beta), 0, s3) + qb_1(math.pi / 2) + qb_2(0.0865)
    else:
        nbooms = int((s3 - boom9_s3) / boomspace_s3) + 1
        qb_booms = 0
        for i in range(8, 8 + nbooms):
            qb_booms += stiff_area * boomcoords_hinge[i, 1]
        qb_skin = tskin * integrate(lambda x: r - x * math.sin(beta), 0, s3)
        qb = c * (qb_skin + qb_booms) + qb_1(math.pi / 2) + qb_2(0.0865)
    return qb


# SECTION 4
def qb_4(s4):  # 0 to 0.4068
    c = -1 / Izz
    boom1_s4 = (Ca - r + boomcoords_hinge[0, 0]) / math.cos(beta)
    boomspace_s4 = -(boomcoords_hinge[0, 0] - boomcoords_hinge[1, 0]) / math.cos(beta)
    if s4 < boom1_s4:
        qb = c * integrate(lambda x: -tskin * math.sin(beta) * x, 0, s4) + qb_3(0.4068)
    else:
        nbooms = int((s4 - boom1_s4) / boomspace_s4) + 1
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
def qb_6(theta):  # -pi/2 to 0
    c = -1 / Izz
    boom6_theta = math.asin(boomcoords[5, 1] / r)
    if theta < boom6_theta:
        qb = c * (tskin * r ** 2) * integrate(math.sin, math.pi / 2, theta) + qb_5(0.0865) + qb_4(0.4068)
    else:
        qb = c * ((tskin * r ** 2) * integrate(math.sin, -math.pi / 2, theta) + stiff_area * boomcoords[5, 1]) + qb_5(
            0.0865) + qb_4(0.4068)
    return qb


# SHEAR CENTRE
def shear_centre(n):
    # REDUNDANT SHEAR

    # CELL I
    # SECTION 1 - base
    int_qb1 = 0
    for i in range(1, n + 1):
        int_qb1 += (qb_1(i * (math.pi / (2 * n))) * ((math.pi * r) / (2 * n))) / tskin

    # SECTION 2 - base (cell I)
    int_qb2I = 0
    for i in range(1, n + 1):
        int_qb2I += (qb_2(r - i * (r / n)) * (r / n)) / tspar

    # SECTION 5 - base (cell I)
    int_qb5I = 0
    for i in range(1, n + 1):
        int_qb5I += (qb_5(i * (r / n)) * (r / n)) / tspar

    # SECTION 6 - base
    int_qb6 = 0
    for i in range(1, n + 1):
        int_qb6 += (qb_6((-math.pi / 2) + i * (math.pi / (2 * n))) * ((math.pi * r) / (2 * n))) / tskin
    print('Done computing cell I contribution.')
    qb_intI = int_qb1 - int_qb2I + int_qb5I + int_qb6

    # CELL II
    # SECTION 2 - base (cell II)
    int_qb2II = 0
    for i in range(1, n + 1):
        int_qb2II += (qb_2(i * (r / n)) * (r / n)) / tspar

    # SECTION 3 - base
    int_qb3 = 0
    for i in range(1, n + 1):
        int_qb3 += (qb_3(i * (l_topskin / n)) * (l_topskin / n)) / tskin

    # SECTION 4 - base
    int_qb4 = 0
    for i in range(1, n + 1):
        int_qb4 += (qb_4(i * (l_topskin / n)) * (l_topskin / n)) / tskin

    # SECTION 5 - base (cell II)
    int_qb5II = 0
    for i in range(1, n + 1):
        int_qb5II += (qb_5(r - i * (r / n)) * (r / n)) / tspar
    print('Done computing cell II contribution.')
    qb_intII = int_qb2II + int_qb3 + int_qb4 - int_qb5II

    # SYSTEM SOLVING
    eqs = np.array([
        [(math.pi * r) / tskin + r / tspar, (-2 * r) / tspar],
        [-2 * r / tspar, (2 * r / tspar) + (2 * l_topskin / tskin)]
    ])
    cs = np.array([-qb_intI, -qb_intII])
    qsI, qsII = np.linalg.solve(eqs, cs)
    print('Redundant shear flow in cell I is', qsI, '[N/m], clockwise positive')
    print('Redundant shear flow in cell II is', qsII, '[N/m], clockwise positive')

    def q1(theta):
        q = qb_1(theta) + qsI
        return q

    def q2(s):
        q = qb_2(s) - qsI + qsII
        return q

    def q3(s):
        q = qb_3(s) + qsII
        return q

    def q4(s):
        q = qb_4(s) + qsII
        return q

    def q5(s):
        q = qb_5(s) + qsI - qsII
        return q

    def q6(theta):
        q = qb_6(theta) + qsI
        return q

    # MOMENT CONTRIBUTION
    # SECTION 1
    int_q1 = 0
    for i in range(1, n + 1):
        int_q1 += (q1(i * (math.pi / (2 * n))) * ((math.pi * r) / (2 * n)))
    m1 = int_q1 * r

    # SECTION 3
    int_q3 = 0
    for i in range(1, n + 1):
        int_q3 += (q3(i * (l_topskin / n)) * (l_topskin / n))
    m3 = int_q3 * math.cos(beta) * r

    # SECTION 4
    int_q4 = 0
    for i in range(1, n + 1):
        int_q4 += (q4(i * (l_topskin / n)) * (l_topskin / n))
    m4 = int_q4 * math.cos(beta) * r

    # SECTION 6
    int_q6 = 0
    for i in range(1, n + 1):
        int_q6 += (q6((-math.pi / 2) + i * (math.pi / (2 * n))) * ((math.pi * r) / (2 * n)))
    m6 = int_q6 * r

    # FINAL COMPUTATION
    xi = (m1 + m3 + m4 + m6)

    return qsI, qsII, q1, q2, q3, q4, q5, q6, xi

#qsI, qsII, q1, q2, q3, q4, q5, q6, xi = shear_centre(1000)

#print(xi)