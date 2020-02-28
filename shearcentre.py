from integration import integrate as S
import sectionproperties as SP
from torsionstiffness import torsionalstiffness
import numpy as np

sect = None


def set_sect(section):
    global sect
    sect = section


# BASIC SHEAR FLOWS - VERTICAL SHEAR FORCE
def qb_1(theta):  # SECTION 1: 0 to pi/2
    c = -1 / sect.Izz
    boom8_theta = np.arcsin(sect.boomcoords[7, 1] / sect.r)
    if theta < boom8_theta:
        qb = c * sect.tskin * sect.r ** 2 * S(np.sin, 0, theta)
    else:
        qb = c * (sect.tskin * sect.r ** 2 * S(np.sin, 0, theta) + sect.stiff_area * sect.boomcoords[7, 1])
    return qb


def qb_2(s2):  # SECTION 2: 0 to r
    c = -1 / sect.Izz
    qb = c * S(lambda x: sect.tspar * x, 0, s2)
    return qb


def qb_3(s3):  # SECTION 3: 0 to l_topskin
    c = -1 / sect.Izz
    boom9_s3 = -sect.boomcoords_hinge[8, 0] / np.cos(sect.beta)
    boomspace_s3 = (sect.boomcoords_hinge[8, 0] - sect.boomcoords_hinge[9, 0]) / np.cos(sect.beta)
    if s3 < boom9_s3:
        qb = c * sect.tskin * S(lambda x: sect.r - x * np.sin(sect.beta), 0, s3) + qb_1(np.pi / 2) + qb_2(sect.r)
    else:
        nbooms = int((s3 - boom9_s3) // boomspace_s3) + 1
        qb_booms = 0
        for i in range(8, 8 + nbooms):
            qb_booms += sect.stiff_area * sect.boomcoords_hinge[i, 1]
        qb_skin = sect.tskin * S(lambda x: sect.r - x * np.sin(sect.beta), 0, s3)
        qb = c * (qb_skin + qb_booms) + qb_1(np.pi / 2) + qb_2(sect.r)
    return qb


def qb_4(s4):  # SECTION 4: 0 to l_topskin
    c = -1 / sect.Izz
    boom1_s4 = (sect.Ca - sect.r + sect.boomcoords_hinge[0, 0]) / np.cos(sect.beta)
    boomspace_s4 = -(sect.boomcoords_hinge[0, 0] - sect.boomcoords_hinge[1, 0]) / np.cos(sect.beta)
    if s4 < boom1_s4:
        qb = c * S(lambda x: -sect.tskin * np.sin(sect.beta) * x, 0, s4) + qb_3(sect.l_topskin)
    else:
        nbooms = int((s4 - boom1_s4) / boomspace_s4) + 1
        qb_booms = 0
        for i in range(0, 0 + nbooms):
            qb_booms += sect.stiff_area * sect.boomcoords_hinge[i, 1]
        qb_skin = S(lambda x: -sect.tskin * np.sin(sect.beta) * x, 0, s4)
        qb = c * (qb_skin + qb_booms) + qb_3(sect.l_topskin)
    return qb


def qb_5(s5):  # SECTION 5: 0 to r
    c = -1 / sect.Izz
    qb = c * S(lambda x: -sect.tspar * x, 0, s5)
    return qb


def qb_6(theta):  # SECTION 6: -pi/2 to 0
    c = -1 / sect.Izz
    boom6_theta = np.arcsin(sect.boomcoords[5, 1] / sect.r)
    if theta < boom6_theta:
        qb = c * sect.tskin * sect.r ** 2 * S(np.sin, -np.pi / 2, theta) + qb_5(sect.r) + qb_4(sect.l_topskin)
    else:
        qb = c * (sect.tskin * sect.r ** 2 * S(np.sin, -np.pi / 2, theta) + sect.stiff_area * sect.boomcoords[
            5, 1]) + qb_5(sect.r) + qb_4(sect.l_topskin)
    return qb


# SHEAR FLOWS - HORIZONTAL SHEAR FORCE
def q_horizontal_1(theta):  # SECTION 1: 0 to pi/2
    c = -1 / sect.Iyy
    boom8_theta = np.arcsin(sect.boomcoords[7, 1] / sect.r)
    if theta < boom8_theta:
        q = c * (sect.tskin * sect.r * S(lambda x: sect.r * np.cos(x) - sect.z_centroid, 0,
                                         theta) + 0.5 * sect.stiff_area * sect.boomcoords[6, 0])
    else:
        q = c * (sect.tskin * sect.r * S(lambda x: sect.r * np.cos(x) - sect.z_centroid, 0, theta) + sect.stiff_area *
                 sect.boomcoords[7, 0])
    return q


def q_horizontal_2(s2):  # SECTION 2: 0 to r
    c = -1 / sect.Iyy
    q = c * S(lambda x: sect.tspar * -sect.z_centroid, 0, s2)
    return q


def q_horizontal_3(s3):  # SECTION 3: 0 to l_topskin
    c = -1 / sect.Iyy
    boom9_s3 = -sect.boomcoords_hinge[8, 0] / np.cos(sect.beta)
    boomspace_s3 = (sect.boomcoords_hinge[8, 0] - sect.boomcoords_hinge[9, 0]) / np.cos(sect.beta)
    if s3 < boom9_s3:
        q = -c * sect.tskin * S(lambda x: x * np.cos(sect.beta) + sect.z_centroid, 0, s3) + q_horizontal_1(
            np.pi / 2) + q_horizontal_2(sect.r)
    else:
        nbooms = int((s3 - boom9_s3) // boomspace_s3) + 1
        q_booms = 0
        for i in range(8, 8 + nbooms):
            q_booms += sect.stiff_area * sect.boomcoords[i, 0]
        q_skin = -sect.tskin * S(lambda x: x * np.cos(sect.beta) + sect.z_centroid, 0, s3)
        q = c * (q_skin + q_booms) + q_horizontal_1(np.pi / 2) + q_horizontal_2(sect.r)
    return q


def q_horizontal_4(s4):  # SECTION 4: 0 to l_topskin
    c = -1 / sect.Iyy
    boom1_s4 = (sect.Ca - sect.r + sect.boomcoords_hinge[0, 0]) / np.cos(sect.beta)
    boomspace_s4 = -(sect.boomcoords_hinge[0, 0] - sect.boomcoords_hinge[1, 0]) / np.cos(sect.beta)
    if s4 < boom1_s4:
        q = c * S(lambda x: sect.tskin * (-(sect.l_spar_to_end + sect.z_centroid) + (x * np.cos(sect.beta))), 0,
                  s4) + q_horizontal_3(sect.l_topskin)
    else:
        nbooms = int((s4 - boom1_s4) / boomspace_s4) + 1
        q_booms = 0
        for i in range(0, 0 + nbooms):
            q_booms += sect.stiff_area * sect.boomcoords[i, 0]
        q_skin = S(lambda x: sect.tskin * (-(sect.l_spar_to_end + sect.z_centroid) + (x * np.cos(sect.beta))), 0, s4)
        q = c * (q_skin + q_booms) + q_horizontal_3(sect.l_topskin)
    return q


def q_horizontal_5(s5):  # SECTION 5: 0 to r
    c = -1 / sect.Iyy
    q = c * S(lambda x: sect.tspar * -sect.z_centroid, 0, s5)
    return q


def q_horizontal_6(theta):  # SECTION 6: -pi/2 to 0
    c = -1 / sect.Iyy
    boom6_theta = np.arcsin(sect.boomcoords[5, 1] / sect.r)
    if theta < boom6_theta:
        q = c * sect.tskin * sect.r * S(lambda x: sect.r * np.cos(x) - sect.z_centroid, -np.pi / 2,
                                        theta) + q_horizontal_5(sect.r) + q_horizontal_4(sect.l_topskin)
    elif theta == 0:
        q = c * (sect.tskin * sect.r * S(lambda x: sect.r * np.cos(x) - sect.z_centroid, -np.pi / 2,
                                         theta) + sect.stiff_area * sect.boomcoords[5, 0] + 0.5 * sect.stiff_area *
                 sect.boomcoords[6, 0]) + q_horizontal_5(sect.r) + q_horizontal_4(sect.l_topskin)
    else:
        q = c * (sect.tskin * sect.r * S(lambda x: sect.r * np.cos(x) - sect.z_centroid, -np.pi / 2,
                                         theta) + sect.stiff_area * sect.boomcoords[5, 0]) + q_horizontal_5(
            sect.r) + q_horizontal_4(sect.l_topskin)
    return q


# SHEAR CENTRE CALCULATION

def shear_centre(n):
    global sect
    if sect is None:
        sect = SP.section()

    # REDUNDANT SHEAR

    # CELL I
    # SECTION 1 - base
    int_qb1 = 0
    for i in range(1, n + 1):
        int_qb1 += qb_1(i * np.pi / (2 * n)) * np.pi * sect.r / (2 * n) / sect.tskin

    # SECTION 2 - base (cell I)
    int_qb2I = 0
    for i in range(1, n + 1):
        int_qb2I += qb_2(sect.r - i * sect.r / n) * sect.r / n / sect.tspar

    # SECTION 5 - base (cell I)
    int_qb5I = 0
    for i in range(1, n + 1):
        int_qb5I += qb_5(i * sect.r / n) * sect.r / n / sect.tspar

    # SECTION 6 - base
    int_qb6 = 0
    for i in range(1, n + 1):
        int_qb6 += qb_6(-np.pi / 2 + i * np.pi / (2 * n)) * np.pi * sect.r / (2 * n) / sect.tskin
    print('Done computing cell I contribution.')
    qb_intI = int_qb1 - int_qb2I + int_qb5I + int_qb6

    # CELL II
    # SECTION 2 - base (cell II)
    int_qb2II = 0
    for i in range(1, n + 1):
        int_qb2II += qb_2(i * sect.r / n) * sect.r / n / sect.tspar

    # SECTION 3 - base
    int_qb3 = 0
    for i in range(1, n + 1):
        int_qb3 += qb_3(i * sect.l_topskin / n) * sect.l_topskin / n / sect.tskin

    # SECTION 4 - base
    int_qb4 = 0
    for i in range(1, n + 1):
        int_qb4 += qb_4(i * sect.l_topskin / n) * sect.l_topskin / n / sect.tskin

    # SECTION 5 - base (cell II)
    int_qb5II = 0
    for i in range(1, n + 1):
        int_qb5II += qb_5(sect.r - i * sect.r / n) * sect.r / n / sect.tspar
    print('Done computing cell II contribution.')
    qb_intII = int_qb2II + int_qb3 + int_qb4 - int_qb5II

    # SYSTEM SOLVING
    eqs = np.array([
        [np.pi * sect.r / sect.tskin + sect.r / sect.tspar, -2 * sect.r / sect.tspar],
        [-2 * sect.r / sect.tspar, 2 * sect.r / sect.tspar + 2 * sect.l_topskin / sect.tskin]
    ])
    cs = np.array([-qb_intI, -qb_intII])
    qsI, qsII = np.linalg.solve(eqs, cs)
    print('Redundant shear flow in cell I is', qsI, '[N/m], clockwise positive')
    print('Redundant shear flow in cell II is', qsII, '[N/m], clockwise positive')

    def q1_v(theta):
        q = qb_1(theta) + qsI
        return q

    def q2_v(s):
        q = qb_2(s) - qsI + qsII
        return q

    def q3_v(s):
        q = qb_3(s) + qsII
        return q

    def q4_v(s):
        q = qb_4(s) + qsII
        return q

    def q5_v(s):
        q = qb_5(s) + qsI - qsII
        return q

    def q6_v(theta):
        q = qb_6(theta) + qsI
        return q

    # MOMENT CONTRIBUTION
    # SECTION 1
    int_q1 = 0
    for i in range(1, n + 1):
        int_q1 += q1_v(i * np.pi / (2 * n)) * np.pi * sect.r / (2 * n)
    m1 = int_q1 * sect.r

    # SECTION 3
    int_q3 = 0
    for i in range(1, n + 1):
        int_q3 += q3_v(i * sect.l_topskin / n) * sect.l_topskin / n
    m3 = int_q3 * np.cos(sect.beta) * sect.r

    # SECTION 4
    int_q4 = 0
    for i in range(1, n + 1):
        int_q4 += q4_v(i * sect.l_topskin / n) * sect.l_topskin / n
    m4 = int_q4 * np.cos(sect.beta) * sect.r

    # SECTION 6
    int_q6 = 0
    for i in range(1, n + 1):
        int_q6 += q6_v(-np.pi / 2 + i * np.pi / (2 * n)) * np.pi * sect.r / (2 * n)
    m6 = int_q6 * sect.r

    # FINAL COMPUTATION
    xi = -(m1 + m3 + m4 + m6)

    # TOTAL EQUATIONS

    q0I, q0II, _ = torsionalstiffness(sect)

    def q1(theta, Sz, Sy, T):
        q = q1_v(theta) * Sy + q_horizontal_1(theta) * Sz + T * q0I
        return q

    def q2(s, Sz, Sy, T):
        q = q2_v(s) * Sy + q_horizontal_2(s) * Sz + T * (q0II - q0I)
        return q

    def q3(s, Sz, Sy, T):
        q = q3_v(s) * Sy + q_horizontal_3(s) * Sz + T * q0II
        return q

    def q4(s, Sz, Sy, T):
        q = q4_v(s) * Sy + q_horizontal_4(s) * Sz + T * q0II
        return q

    def q5(s, Sz, Sy, T):
        q = q5_v(s) * Sy + q_horizontal_5(s) * Sz + T * (q0I - q0II)
        return -q

    def q6(theta, Sz, Sy, T):
        q = q6_v(theta) * Sy + q_horizontal_6(theta) * Sz + T * q0I
        return q

    return qsI, qsII, q1, q2, q3, q4, q5, q6, xi

# qsI, qsII, q1, q2, q3, q4, q5, q6, xi = shear_centre(1000)

##Verification
# our_zc = xi - sect.r
# print("our zc:", our_zc)
# ans = -0.09185594953325857
# print("should be:", ans)
# print("error:", our_zc - ans, round((our_zc - ans) * 100 / ans, 4), "%")
