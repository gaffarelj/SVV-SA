import sectionproperties as SP
from shearcentre import shear_centre
import math
from interpolation import get_load
import numpy as np
from integration import integrate as S
import matplotlib.pyplot as plt
plt.gcf().subplots_adjust(left=0.25)


xi, J = None, None
La = None
x1, x2, x3, xa = None, None, None, None
d1, d3 = None, None
r, G, E = None, None, None
theta, P = None, None
Izz, Iyy = None, None
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = None, None, None, None, None, None, None, None, None, None, None, None
f2, f3, f5, fq, fd = None, None, None, None, None
load = None


# MACAULAY
def set_vars(xi_i, J_i, r_i, Izz_i, Iyy_i, G_i=28e9, E_i=73.1e9, La_i=1.691, 
             x1_i=0.149, x2_i=0.554, x3_i=1.541, d1_i=0.00681, d3_i=0.0203, 
             xa_i=0.272, theta_i=np.radians(26), P_i=37.9e3, load_i=None):
    global xi, J, La, x1, x2, x3, xa, d1, d3, r, G, E, theta, P, Izz, Iyy, load
    xi, J, La, x1, x2, x3, xa, d1, d3, r, G, E, theta, P, Izz, Iyy, load = \
        xi_i, J_i, La_i, x1_i, x2_i, x3_i, xa_i, d1_i, d3_i, r_i, G_i, E_i, theta_i, P_i, Izz_i, Iyy_i, load_i

# Solve for unknowns
# Unknowns are, in order, [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5]
def system(power=10, power_t=6):
    global Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5
    global f2, f3, f5, fq, fd
    print("Macaulay:")
    a = coefficients_aero()
    b = coefficients_tau()

    print("\t Setting up matrices...")
    eq1 = np.array([0, 0, 0, -(La - x1), -(La - x2), -(La - x3), -math.cos(theta) * (La - (x2 - (xa / 2))), 0, 0, 0, 0, 0])
    res1 = np.array([-P * math.cos(theta) * (La - (x2 + (xa / 2)))])
    
    new_c2 = np.zeros((power+1, 1))
    for i in range(power+1):
        new_c2[i, 0] = a[i, 0] / (i + 1)
    if power == 0:
        f2 = lambda x: new_c2[0, 0] * x
    else:
        f2 = lambda x: new_c2[0, 0] * x + new_c2[1, 0] * x ** 2 + new_c2[2, 0] * x ** 3 + new_c2[3, 0] * x ** 4 + \
                   new_c2[4, 0] * x ** 5 + new_c2[5, 0] * x ** 6 + new_c2[6, 0] * x ** 7 + new_c2[7, 0] * x ** 8 + \
                   new_c2[8, 0] * x ** 9 + new_c2[9, 0] * x ** 10 + new_c2[10, 0] * x ** 11
    integral2 = S(f2, 0, La)
    eq2 = np.array([(La - x1), (La - x2), (La - x3), 0, 0, 0, math.sin(theta) * (La - (x2 - (xa / 2))), 0, 0, 0, 0, 0])
    res2 = np.array([P * math.sin(theta) * (La - (x2 + (xa / 2))) + integral2])
    
    if power_t == 2:
        f3 = lambda x: b[0, 0] + b[1, 0] * x
    else:
        f3 = lambda x: b[0, 0] + b[1, 0] * x + b[2, 0] * x ** 2 + b[3, 0] * x ** 3 + b[4, 0] * x ** 4 + b[5, 0] * x ** 5 + b[
        6, 0] * x ** 6
    integral3 = S(f3, 0, La)
    eq3 = np.array([-xi * (La - x1), -xi * (La - x2), -xi * (La - x3), 0, 0, 0,
                    (math.cos(theta) * r - xi * math.sin(theta)) * (La - (x2 - (xa / 2))), 0, 0, 0, 0, 0])
    res3 = np.array([(math.cos(theta) * r - xi * math.sin(theta)) * P * (La - (x2 + (xa / 2))) - integral3])
    
    eq4 = np.array([0, 0, 0, -1, -1, -1, -math.cos(theta), 0, 0, 0, 0, 0])
    res4 = np.array([-P * math.cos(theta)])
    
    if power == 0:
        f5 = lambda x: a[0, 0]
    else:
        f5 = lambda x: a[0, 0] + a[1, 0] * x + a[2, 0] * x ** 2 + a[3, 0] * x ** 3 + a[4, 0] * x ** 4 + a[5, 0] * x ** 5 + a[
        6, 0] * x ** 6 + a[7, 0] * x ** 7 + a[8, 0] * x ** 8 + a[9, 0] * x ** 9 + a[10, 0] * x ** 10
    integral5 = S(f5, 0, La)
    eq5 = np.array([1, 1, 1, 0, 0, 0, math.sin(theta), 0, 0, 0, 0, 0])
    res5 = np.array([P * math.sin(theta) + integral5])

    # Quadruple q function
    new_cq = np.zeros((power+1, 1))
    for i in range(power+1):
        new_cq[i, 0] = a[i, 0] / ((i + 1) * (i + 2) * (i + 3))
    if power == 0:
        fq = lambda x: new_cq[0, 0] * x ** 3
    else:
        fq = lambda x: new_cq[0, 0] * x ** 3 + new_cq[1, 0] * x ** 4 + new_cq[2, 0] * x ** 5 + new_cq[3, 0] * x ** 6 + \
                   new_cq[4, 0] * x ** 7 + new_cq[5, 0] * x ** 8 + new_cq[6, 0] * x ** 9 + new_cq[7, 0] * x ** 10 + \
                   new_cq[8, 0] * x ** 11 + new_cq[9, 0] * x ** 12 + new_cq[10, 0] * x ** 13

    # Double tau function
    new_ct = np.zeros((power_t+1, 1))
    for i in range(power_t+1):
        new_ct[i, 0] = b[i, 0] / (i + 1)
    if power_t == 2:
        fd = lambda x: new_ct[0, 0] * x + new_ct[1, 0] * x ** 2
    else:
        fd = lambda x: new_ct[0, 0] * x + new_ct[1, 0] * x ** 2 + new_ct[2, 0] * x ** 3 + new_ct[3, 0] * x ** 4 + new_ct[
        4, 0] * x ** 5 + new_ct[5, 0] * x ** 6 + new_ct[6, 0] * x ** 7
    
    eq6 = np.array([0, 0, 0, 0, 0, 0, 0, x1, 1, 0, 0, xi])
    res6 = np.array([d1 * math.cos(theta) - (1 / (E * Izz)) * S(fq, 0, x1) - (1 / (G * J)) * S(fd, 0, x1) * xi])
    
    eq7 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, x1, 1, 0])
    res7 = np.array([-d1 * math.sin(theta)])
    
    eq8 = np.array([-(1 / (6 * E * Izz)) * (x2 - x1) ** 3 - (xi ** 2 / (G * J)) * (x2 - x1), 0, 0, 0, 0, 0,
                    -(1 / (6 * E * Izz)) * math.sin(theta) * (xa / 2) ** 3 + (xi / (G * J)) * (
                            math.cos(theta) * r - xi * math.sin(theta)) * (xa / 2), x2, 1, 0, 0, xi])
    res8 = np.array([-(1 / (E * Izz)) * S(fq, 0, x2) - (1 / (G * J)) * S(fd, 0, x2) * xi])
    
    eq9 = np.array(
        [0, 0, 0, (1 / (6 * E * Iyy)) * (x2 - x1) ** 3, 0, 0, (math.cos(theta) / (6 * E * Iyy)) * (xa / 2) ** 3,
         0, 0, x2, 1, 0])
    res9 = np.array([0])
    
    eq10 = np.array([-(1 / (6 * E * Izz)) * (x3 - x1) ** 3 - (xi ** 2 / (G * J)) * (x3 - x1),
                     -(1 / (6 * E * Izz)) * (x3 - x2) ** 3 - (xi ** 2 / (G * J)) * (x3 - x2), 0, 0, 0, 0,
                     -(math.sin(theta) / (6 * E * Izz)) * (x3 - (x2 - (xa / 2))) ** 3 + (xi / (G * J)) * (
                             math.cos(theta) * r - xi * math.sin(theta)) * (x3 - (x2 - (xa / 2))), x3, 1, 0, 0,
                     xi])
    res10 = np.array([d3 * math.cos(theta) - ((P * math.sin(theta)) / (6 * E * Izz)) * (x3 - (x2 + (xa / 2))) ** 3 + \
        (1 / (E * Izz)) * S(fq, 0, x3) + \
               (xi / (G * J)) * (math.cos(theta) * r - xi * math.sin(theta)) * P * (x3 - (x2 + (xa / 2))) -\
              (xi / (G * J)) * S(fd, 0, x3)])
    
    eq11 = np.array([0, 0, 0, (1 / (6 * E * Iyy)) * (x3 - x1) ** 3, (1 / (6 * E * Iyy)) * (x3 - x2) ** 3, 0,
                     (math.cos(theta) / (6 * E * Iyy)) * (x3 - (x2 - (xa / 2))) ** 3, 0, 0, x3, 1, 0])
    res11 = np.array([-d3 * math.sin(theta) + ((P * math.cos(theta)) / (6 * E * Iyy)) * (x3 - (x2 + (xa / 2))) ** 3])
    
    eq12 = np.array([0, 0, 0, (1 / (6 * E * Iyy)) * (x2 - (xa / 2) - x1) ** 3, 0, 0, 0, 0, 0, (x2 - (xa / 2)), 1, 0])
    res12 = np.array([0])

    # Solve the system
    A = np.array([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12])
    y = np.array([res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, res11, res12])
    print("Macaulay: solving system...")
    Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = np.linalg.lstsq(A, y, rcond=None)[0]
    residuals = y - np.dot(A, [Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5])
    #print(residuals)
    #input()
    #Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = np.linalg.solve(A, y)
    return Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5


# INTERPOLATION

# Aerodynamic loading
# a0 + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4 + a5 * x^5 [...]
def coefficients_aero(power=10, chord_steps=300, span_steps=150):
    print("\t Computing Aero Loading...")
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

# Torque
def coefficients_tau(power=6, chord_steps=300, span_steps=150):
    print("\t Computing Torque...")
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

def q_aero(x, chord_steps=300, span_steps=150):  # x goes from 0 to span steps-1
    global load
    if load is None:
        section_coordinates, section_loads = get_load(n_chord=chord_steps, n_span=span_steps, do_plot=False)
        load = section_coordinates, section_loads
    else:
        section_coordinates, section_loads = load
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
    M = -S(f2, 0, x)
    if x > x1:
        M += Ry1 * (x - x1)
    if x > x2 - (xa / 2):
        M += Fa * math.sin(theta) * (x - (x2 - (xa / 2)))
    if x > x2:
        M += Ry2 * (x - x2)
    if x > x2 + (xa / 2):
        M -= P * math.sin(theta) * (x - (x2 + (xa / 2)))
    if x > x3:
        M += Ry3 * (x - x3)
    return M

def T(x):
    T = S(f3, 0, x)
    if x > x1:
        T += Ry1 * (-xi)
    if x > x2 - (xa / 2):
        T += (math.cos(theta) * r - xi * math.sin(theta)) * Fa
    if x > x2:
        T += Ry2 * (-xi)
    if x > x2 + (xa / 2):
        T -= (math.cos(theta) * r - xi * math.sin(theta)) * P
    if x > x3:
        T += Ry3 * (-xi)
    return -T

def Sz(x):
    F = 0
    if x > x1:
        F -= Rz1
    if x > x2 - (xa / 2):
        F -= math.cos(theta) * Fa
    if x > x2:
        F -= Rz2
    if x > x2 + (xa / 2):
        F += math.cos(theta) * P
    if x > x3:
        F -= Rz3
    return F

def Sy(x):
    F = -S(f5, 0, x)
    if x > x1:
        F += Ry1
    if x > x2 - (xa / 2):
        F += math.sin(theta) * Fa
    if x > x2:
        F += Ry2
    if x > x2 + (xa / 2):
        F -= math.sin(theta) * P
    if x > x3:
        F += Ry3
    return F

def v(x):
    core = -S(fq, 0, x)
    if x > x1:
        core += (Ry1/6) * (x - x1) ** 3
    if x > x2 - (xa / 2):
        core += math.sin(theta) * (Fa/6) * (x - (x2 - (xa / 2))) ** 3
    if x > x2:
        core += (Ry2/6) * (x - x2) ** 3
    if x > x2 + (xa / 2):
        core -= math.sin(theta) * (P/6) * (x - (x2 + (xa / 2))) ** 3
    if x > x3:
        core += (Ry3/6) * (x - x3) ** 3
    v = -(1 / (E * Izz)) * core + C1 * x + C2
    return v

def w(x):
    core = 0
    if x > x1:
        core -= (Rz1/6) * (x - x1) ** 3
    if x > x2 - (xa / 2):
        core -= math.cos(theta) * (Fa/6) * (x - (x2 - (xa / 2))) ** 3
    if x > x2:
        core -= (Rz2/6) * (x - x2) ** 3
    if x > x2 + (xa / 2):
        core += math.cos(theta) * (P/6) * (x - (x2 + (xa / 2))) ** 3
    if x > x3:
        core -= (Rz3/6) * (x - x3) ** 3
    v = -(1 / (E * Iyy)) * core + C3 * x + C4
    return v

def alpha(x): #Check!
    core = S(fd, 0, x)
    if x > x1:
        core += -xi * Ry1 * (x - x1)
    if x > x2 - (xa / 2):
        core += (math.cos(theta) * r - xi * math.sin(theta)) * Fa * (x - (x2 - (xa / 2)))
    if x > x2:
        core -= -xi * Ry2 * (x - x2)
    if x > x2 + (xa / 2):
        core += (math.cos(theta) * r - xi * math.sin(theta)) * P * (x - (x2 + (xa / 2)))
    if x > x3:
        core -= -xi * Ry3 * (x - x3)
    alpha = -(1 / (G * J)) * core + C5
    angle = math.degrees(alpha)
    return -alpha

def plot_result(f, legend, title, show_plot=False, dx=0.005):
    x = np.arange(0.0, La, dx)
    y = [f(xi) for xi in x]
    plt.plot(x, y)
    if show_plot: plt.show()
    plt.ylabel(title)
    plt.xlabel("x [m]")
    plt.savefig(f"plots/macaulay/{legend}.pdf", bbox_inches='tight')
    plt.close()

def do_plots():
    plot_result(alpha, "alpha", "Twist [rad]")
    plot_result(w, "w", "Displacement in the z-direction [m]")
    plot_result(v, "v", "Displacement in the y-direction [m]")
    plot_result(Sy, "Sy", "Shear in the y-direction [N]")
    plot_result(Sz, "Sz", "Shear in the z-direction [N]")
    plot_result(T, "T", "Torque around x [Nm]")
    plot_result(My, "My", "Bending moment around z [Nm]")
    plot_result(Mz, "Mz", "Bending moment around y [Nm]")

    # Plot bendings
    xs = np.arange(0, 1.691, 0.001)
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('x [m]')
    color = 'tab:red'
    ax1.set_ylabel('bending moment around z [Nm]', color=color)
    ax1.plot(xs, [Mz(x) for x in xs], color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('bending moment around y [Nm]', color=color)
    ax2.plot(xs, [My(x) for x in xs], color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.savefig(f"plots/macaulay/bendings_n.pdf", bbox_inches='tight')
    plt.close()