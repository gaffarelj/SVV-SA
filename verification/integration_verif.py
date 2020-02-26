import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import numpy as np
from integration import integrate as S

# Polynomial test
def f_a(x):
	return 10 * x**3 - 3 * x**2 - 20
S_a1 = S(f_a, -50, 50, dx=0.1)
S_a2 = S(f_a, -50, 50, dx=0.001)
S_a = -252000
print("polynomial", S_a, (S_a1-S_a)/S_a*100, (S_a2-S_a)/S_a*100)

# Trigonometric test
def f_b(x):
	return 10 * np.sin(x) - 3 * np.cos(x)
S_b1 = S(f_b, -50, 50, dx=0.1)
S_b2 = S(f_b, -50, 50, dx=0.001)
S_b = 1.574249122
print("trigo", S_b, (S_b1-S_b)/S_b*100, (S_b2-S_b)/S_b*100)

# Rational test
def f_c(x):
	return (5 * x **2 + 2 * x + 7) / (x**2 - 5*x - 14)
S_c1 = S(f_c, 8, 20, dx=0.1)
S_c2 = S(f_c, 8, 20, dx=0.001)
S_c = 133.7935566
print("rational", S_c, (S_c1-S_c)/S_c*100, (S_c2-S_c)/S_c*100)

# Exponential test
def f_d(x):
	return 4 ** x
S_d1 = S(f_d, -50, 50, dx=0.1)
S_d2 = S(f_d, -50, 50, dx=0.001)
S_d = 9.144166173e29
print("exp", S_d, (S_d1-S_d)/S_d*100, (S_d2-S_d)/S_d*100)

# Logaritmic test
def f_e(x):
	return 4 * np.log(x) + 3
S_e1 = S(f_e, 1, 50, dx=0.1)
S_e2 = S(f_e, 1, 50, dx=0.001)
S_e = 733.4046011
print("log", S_e, (S_e1-S_e)/S_e*100, (S_e2-S_e)/S_e*100)

# Linear test
def f_f(x):
	return 10 * x - 20
S_f1 = S(f_f, -50, 50, dx=0.1)
S_f2 = S(f_f, -50, 50, dx=0.001)
S_f = -2000
print("linear", S_f, (S_f1-S_f)/S_f*100, (S_f2-S_f)/S_f*100)