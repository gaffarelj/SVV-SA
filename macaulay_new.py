import sectionproperties as SP
from shearcentre import shear_centre
import math
from interpolation import get_load
import numpy as np
from integration import integrate as S
import matplotlib.pyplot as plt

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
# qsI, qsII, q1, q2, q3, q4, q5, q6, xi = shear_centre(1000)
xi = -0.007513567161803937
J = 0.000187829  # from verification model
G = 28e9
E = 73.1e9

Ha = 0.173
Ca = 0.484
La = 1.691
tspar = 0.0025
hstiff = 0.014
tstiff = 0.0012
wstiff = 0.018
beta = math.atan(8.65 / 39.75)
l_topskin = 0.3975 / math.cos(beta)
r = Ha / 2
x1 = 0.149
x2 = 0.554
x3 = 1.541
xa = 0.272
d1 = 0.00681
d3 = 0.0203
theta = math.radians(26)
P = 37.9e3


def q_aero(x, chord_steps=300, span_steps=150):  # x goes from 0 to span steps-1
    section_coordinates, section_loads = get_load(n_chord=chord_steps, n_span=span_steps, do_plot=False)
    int_x = np.zeros((1, 3))
    sum_chord = 0
    sum_moments = 0
    int_x[0, 0] = section_coordinates[x]
    for j in range(chord_steps):
        z_step = abs(section_loads[x][1][j][0] - section_loads[x][1][j - 1][0])
        sum_chord += section_loads[x][1][j][1] * z_step
        sum_moments += section_loads[x][1][j][1] * z_step * section_loads[x][1][j][0]
    int_x[0, 1] = sum_chord
    int_x[0, 2] = sum_moments / sum_chord
    return int_x


def tau(x, chord_steps=300, span_steps=150):  # about shear centre, clockwise negative
    tau = np.zeros((1, 2))
    x_position = q_aero(x, chord_steps=chord_steps, span_steps=span_steps)
    tau[0, 0] = x_position[0, 0]
    tau[0, 1] = x_position[0, 2] * (x_position[0, 1] + xi)
    return tau


# INTERPOLATION

# Aerodynamic loading
# a0 + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4 + a5 * x^5 [...]
def coefficients_aero(power=10, chord_steps=300, span_steps=150):
    V = np.zeros((span_steps, power + 1))
    f = np.zeros((span_steps, 1))
    for i in range(span_steps):
        current = q_aero(i, chord_steps=chord_steps, span_steps=span_steps)
        x_coord = current[0, 0]
        load = current[0, 1]
        f[i, 0] = load
        for j in range(power + 1):
            V[i, j] = x_coord ** j
    a = np.dot(np.dot(np.linalg.inv(np.dot(V.T, V)), V.T), f)
    SSD = sum((np.subtract(np.dot(V, a), f)) ** 2)
    avg_error = np.mean(np.subtract(np.dot(V, a), f))
    return a


"""""
# Visual inspection of the interpolation

x = np.linspace(0.0, 1.691, 150)
y = np.zeros((150, 1))
coef = coefficients_aero(10)
for i in range(150):
    y_sum = 0
    for j in range(11):
        y_sum += coef[j, 0] * x[i] ** j
    y[i, 0] = y_sum

x2 = np.zeros((150,1))
y2 = np.zeros((150,1))
for j in range(150):
    point = q_aero(j)
    x2[j,0] = point[0, 0]
    y2[j,0] = point[0, 1]

plt.scatter(x2,y2)
plt.plot(x, y)
plt.show()
"""""


# Torque

def coefficients_tau(power=6, chord_steps=300, span_steps=150):
    V = np.zeros((span_steps, power + 1))
    f = np.zeros((span_steps, 1))
    for i in range(span_steps):
        current = tau(i, chord_steps=chord_steps, span_steps=span_steps)
        x_coord = current[0, 0]
        load = current[0, 1]
        f[i, 0] = load
        for j in range(power + 1):
            V[i, j] = x_coord ** j
    a = np.dot(np.dot(np.linalg.inv(np.dot(V.T, V)), V.T), f)
    SSD = sum((np.subtract(np.dot(V, a), f)) ** 2)
    avg_error = np.mean(np.subtract(np.dot(V, a), f))
    return a


"""""
# Visual inspection of the interpolation

x = np.linspace(0.0, 1.691, 150)
y = np.zeros((150, 1))
coef = coefficients_tau(6)
for i in range(150):
    y_sum = 0
    for j in range(7):
        y_sum += coef[j, 0] * x[i] ** j
    y[i, 0] = y_sum

x2 = np.zeros((150,1))
y2 = np.zeros((150,1))
for j in range(150):
    point = tau(j)
    x2[j,0] = point[0, 0]
    y2[j,0] = point[0, 1]

plt.scatter(x2,y2)
plt.plot(x, y)
plt.show()
"""""

# MACAULAY

# Solve for unknowns
# Unknowns are, in order, [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5]

a = coefficients_aero()
b = coefficients_tau()

# Equation 1
eq1 = np.array([0, 0, 0, -(La - x1), -(La - x2), -(La - x3), -math.cos(theta) * (La - (x2 - (xa / 2))), 0, 0, 0, 0, 0])
res1 = np.array([-P * math.cos(theta) * (La - (x2 + (xa / 2)))])

# Equation 2
new_c2 = np.zeros((11, 1))
for i in range(11):
    new_c2[i, 0] = a[i, 0] / (i + 1)
f2 = lambda x: new_c2[0, 0] * x + new_c2[1, 0] * x ** 2 + new_c2[2, 0] * x ** 3 + new_c2[3, 0] * x ** 4 + \
               new_c2[4, 0] * x ** 5 + new_c2[5, 0] * x ** 6 + new_c2[6, 0] * x ** 7 + new_c2[7, 0] * x ** 8 + \
               new_c2[8, 0] * x ** 9 + new_c2[9, 0] * x ** 10 + new_c2[10, 0] * x ** 11
integral2 = S(f2, 0, La)
eq2 = np.array([(La - x1), (La - x2), (La - x3), 0, 0, 0, math.sin(theta) * (La - (x2 - (xa / 2))), 0, 0, 0, 0, 0])
res2 = np.array([P * math.sin(theta) * (La - (x2 + (xa / 2))) + integral2])

# Equation 3
f3 = lambda x: b[0, 0] + b[1, 0] * x + b[2, 0] * x ** 2 + b[3, 0] * x ** 3 + b[4, 0] * x ** 4 + b[5, 0] * x ** 5 + b[
    6, 0] * x ** 6
integral3 = S(f3, 0, La)
eq3 = np.array([-xi * (La - x1), -xi * (La - x2), -xi * (La - x3), 0, 0, 0,
                (math.cos(theta) * (Ha / 2) - xi * math.sin(theta)) * (La - (x2 - (xa / 2))), 0, 0, 0, 0, 0])
res3 = np.array([(math.cos(theta) * (Ha / 2) - xi * math.sin(theta)) * P * (La - (x2 + (xa / 2))) - integral3])

# Equation 4
eq4 = np.array([0, 0, 0, -1, -1, -1, -math.cos(theta), 0, 0, 0, 0, 0])
res4 = np.array([-P * math.cos(theta)])

# Equation 5
f5 = lambda x: a[0, 0] + a[1, 0] * x + a[2, 0] * x ** 2 + a[3, 0] * x ** 3 + a[4, 0] * x ** 4 + a[5, 0] * x ** 5 + a[
    6, 0] * x ** 6 + a[7, 0] * x ** 7 + a[8, 0] * x ** 8 + a[9, 0] * x ** 9 + a[10, 0] * x ** 10
integral5 = S(f5, 0, La)
eq5 = np.array([1, 1, 1, 0, 0, 0, math.sin(theta), 0, 0, 0, 0, 0])
res5 = np.array([P * math.sin(theta) + integral5])

# Quadruple q function
new_cq = np.zeros((11, 1))
for i in range(11):
    new_cq[i, 0] = a[i, 0] / ((i + 1) * (i + 2) * (i + 3))
fq = lambda x: new_cq[0, 0] * x ** 3 + new_cq[1, 0] * x ** 4 + new_cq[2, 0] * x ** 5 + new_cq[3, 0] * x ** 6 + \
               new_cq[4, 0] * x ** 7 + new_cq[5, 0] * x ** 8 + new_cq[6, 0] * x ** 9 + new_cq[7, 0] * x ** 10 + \
               new_cq[8, 0] * x ** 11 + new_cq[9, 0] * x ** 12 + new_cq[10, 0] * x ** 13

# Double tau function
new_ct = np.zeros((7, 1))
for i in range(7):
    new_ct[i, 0] = b[i, 0] / (i + 1)
fd = lambda x: new_ct[0, 0] * x + new_ct[1, 0] * x ** 2 + new_ct[2, 0] * x ** 3 + new_ct[3, 0] * x ** 4 + new_ct[
    4, 0] * x ** 5 + new_ct[5, 0] * x ** 6 + new_ct[6, 0] * x ** 7

# Equation 6
eq6 = np.array([0, 0, 0, 0, 0, 0, 0, x1, 1, 0, 0, xi])
res6 = np.array([d1 * math.cos(theta) - (1 / (E * sect.Izz)) * S(fq, 0, x1) - (1 / (G * J)) * S(fd, 0, x1) * xi])

# Equation 7
eq7 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, x1, 1, 0])
res7 = np.array([-d1 * math.sin(theta)])

# Equation 8
eq8 = np.array([-(1 / (6 * E * sect.Izz)) * (x2 - x1) ** 3 - (xi ** 2 / (G * J)) * (x2 - x1), 0, 0, 0, 0, 0,
                -(1 / (6 * E * sect.Izz)) * math.sin(theta) * (xa / 2) ** 3 + (xi / (G * J)) * (
                        math.cos(theta) * (Ha / 2) - xi * math.sin(theta)) * (xa / 2), x2, 1, 0, 0, xi])
res8 = np.array([-(1 / (E * sect.Izz)) * S(fq, 0, x2) - (1 / (G * J)) * S(fd, 0, x2) * xi])

# Equation 9
eq9 = np.array(
    [0, 0, 0, (1 / (6 * E * sect.Izz)) * (x2 - x1) ** 3, 0, 0, (math.cos(theta) / (6 * E * sect.Izz)) * (xa / 2) ** 3,
     0, 0, x2, 1, 0])
res9 = np.array([0])

# Equation 10
eq10 = np.array([-(1 / (6 * E * sect.Izz)) * (x3 - x1) ** 3 - (xi ** 2 / (G * J)) * (x3 - x1),
                 -(1 / (6 * E * sect.Izz)) * (x3 - x2) ** 3 - (xi ** 2 / (G * J)) * (x3 - x2), 0, 0, 0, 0,
                 -(math.sin(theta) / (6 * E * sect.Izz)) * (x3 - (x2 - (xa / 2))) ** 3 + (xi / (G * J)) * (
                         math.cos(theta) * (Ha / 2) - xi * math.sin(theta)) * (x3 - (x2 - (xa / 2))), x3, 1, 0, 0,
                 xi])
res10 = np.array([d3 * math.cos(theta) - ((P * math.sin(theta)) / (6 * E * sect.Izz)) * (x3 - (x2 + (xa / 2))) ** 3 + (
        1 / (E * sect.Izz)) * S(fq, 0, x3) + (xi / (G * J)) * (
                          math.cos(theta) * (Ha / 2) - xi * math.sin(theta)) * P * (x3 - (x2 + (xa / 2))) - (
                          xi / (G * J)) * S(fd, 0, x3)])

# Equation 11
eq11 = np.array([0, 0, 0, (1 / (6 * E * sect.Iyy)) * (x3 - x1) ** 3, (1 / (6 * E * sect.Iyy)) * (x3 - x2) ** 3, 0,
                 (math.cos(theta) / (6 * E * sect.Iyy)) * (x3 - (x2 - (xa / 2))) ** 3, 0, 0, x3, 1, 0])
res11 = np.array([-d3 * math.sin(theta) + ((P * math.cos(theta)) / (6 * E * sect.Iyy)) * (x3 - (x2 + (xa / 2))) ** 3])

# Equation 12
eq12 = np.array([0, 0, 0, (1 / (6 * E * sect.Iyy)) * (x2 - (xa / 2) - x1) ** 3, 0, 0, 0, 0, 0, (x2 - (xa / 2)), 1, 0])
res12 = np.array([0])

# Solve the system
A = np.array([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12])
y = np.array([res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, res11, res12])
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = np.linalg.lstsq(A, y, rcond=None)[0]


def My(x):
    M = 0
    if x > x1:
        M -= Rz1 * (x - x1)
    if x > x2 - (xa / 2):
        M -= Fa * math.cos(theta) * (x - (x2 - (xa / 2)))
    if x > x2:
        M -= Rz2 * (x - x2)
    if x > x2 + (xa / 2):
        M += P * math.cos(theta) * (x - (x2 + (xa / 2)))
    if x > x3:
        M -= Rz3 * (x - x3)
    return M

def Mz(x):
    M = S(f2, 0, x)
    if x > x1:
        M += Ry1 * (x - x1)
    if x > x2 - (xa / 2):
        M += Fa * math.sin(theta) * (x - (x2 - (xa / 2)))
    if x > x2:
        M += Ry2 * (x - x2)
    if x > x2 + (xa / 2):
        M -= P * math.sin(theta) * (x - (x2 + (xa / 2)))
    if x > x3:
        M += Rz3 * (x - x3)


x = np.linspace(0.0, 1.691, 500)
y = np.zeros((500, 1))
for i in range(500):
    y[i, 0] = Mz(x[i])

plt.plot(x, y)
plt.show()