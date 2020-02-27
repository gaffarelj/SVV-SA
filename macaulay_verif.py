import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import torsionstiffness as TS
import stress as STR
import max_stress as MS
import numpy as np
import math

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
_, _, J = TS.torsionalstiffness(sect)
MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()
G = 28e9
E = 73.1e9
La = 1.691
x1 = 0.149
x2 = 0.554
x3 = 1.541
d1 = 0.00681
d3 = 0.0203
xa = 0.272
theta = np.radians(26)
P = 37.9e3
Iyy = sect.Iyy
Izz = sect.Izz

# CHECKING BOUNDARY CONDITIONS

# Equation 1 - My(L) = 0
MyL = -Rz1 * (La - x1) - Rz2 * (La - x2) - Rz3 * (La - x3) - Fa * math.cos(theta) * (
        La - (x2 - (xa / 2))) + P * math.cos(theta) * (La - (x2 + (xa / 2)))
print("My(L) is", MyL, "- should be 0")

# Equation 4 - Sz(L) = 0
SzL = -Rz1 - Rz2 - Rz3 - Fa * math.cos(theta) + P * math.cos(theta)
print("Sz(L) is", SzL, "- should be 0")

# Equation 7 - w(x1) = -d1 * sin(theta)
wx1 = C3 * x1 + C4
wx1r = -d1 * math.sin(theta)
print("w(x1) is", wx1, "- should be", wx1r)

# Equation 9 - w(x2) = 0
wx2 = (1 / (E * Iyy)) * ((Rz1 / 6) * (x2 - x1) ** 3 + (Fa / 6) * math.cos(theta) * (xa / 2) ** 3) + C3 * x2 + C4
print("W(x2) is", wx2, "- should be 0")

# Equation 11 - w(x3) = -d3 * sin(theta)
wx3 = -(1 / (E * Iyy)) * (-(Rz1 / 6) * (x3 - x1) ** 3 - (Rz2 / 6) * (x3 - x2) ** 3 - (Fa / 6) * math.cos(theta) * (
        x3 - (x2 - (xa / 2))) ** 3 + (P / 6) * math.cos(theta) * (x3 - (x2 + (xa / 2))) ** 3) + C3 * x3 + C4
wx3r = -d3 * math.sin(theta)
print("W(x3) is", wx3, "- should be", wx3r)

# Equation 9 - w(x2) = 0
wx2a = -(1 / (E * Iyy)) * (-(Rz1 / 6) * (x2 - (xa / 2) - x1) ** 3) + C3 * (x2 - (xa / 2)) + C4
print("W(x2 - (xa/2)) is", wx2a, "- should be 0")

# ORDER OF MAGNITUDE & DIRECTION
print("Ry1 is", round(float(Ry1) / 1000, 1), "[kN]")
print("Ry2 is", round(float(Ry2) / 1000, 1), "[kN]")
print("Ry3 is", round(float(Ry3) / 1000, 1), "[kN]")
print("Rz1 is", round(float(Rz1) / 1000, 1), "[kN]")
print("Rz2 is", round(float(Rz2) / 1000, 1), "[kN]")
print("Rz3 is", round(float(Rz3) / 1000, 1), "[kN]")
print("Fa is", round(float(Fa) / 1000, 1), "[kN]")
