from sectionproperties import get_geometry
from shearcentre import shear_centre
import math
from interpolation import get_load
import numpy as np

Am, z_centroid, Iyy, Izz, tskin, boomcoords, boomcoords_hinge, stiff_area = get_geometry()
# qsI, qsII, q1, q2, q3, q4, q5, q6, xi = shear_centre(1000)

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

load_array = np.array(get_load(n_chord=15, n_span=300, do_plot=True))
#print(load_array)


def q_aero(x, chord_steps=150, span_steps=300):         # x goes from 0 to La
    load_array = np.array(get_load(n_chord=chord_steps, n_span=span_steps, do_plot=False))
    int_x = np.zeros((span_steps, 1))
    x_step = La/span_steps
    z_step = Ca/chord_steps
    for i in range(span_steps):
        sum_chord = 0
        for j in range(chord_steps):
            sum_chord += load_array[i, 2] * z_step
        int_x[i, 0] = sum_chord
    x_index = int(x/x_step)
    return x_index

#print(q_aero(La))
#print(load_array[,149,0])