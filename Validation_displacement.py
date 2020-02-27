import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os.path
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import torsionstiffness as TS
import stress as STR
import interpolation as I

def validation_booms():
    zy1 = [0.1025, 0]
    zy2 = [0.06558, 0.07878]
    zy3 = [-0.01831, 0.09876]
    zy4 = [-0.10635, 0.08081]
    zy5 = [-0.19438, 0.06285]
    zy6 = [-0.28241, 0.04489]
    zy7 = [-0.37045, 0.02694]
    zy8 = [-0.45848, 0.00898]
    zy9 = [-0.45848, -0.00898]
    zy10 = [-0.37045, -0.02694]
    zy11 = [-0.28241, -0.04489]
    zy12 = [-0.19438, -0.06285]
    zy13 = [-0.10635, -0.08081]
    zy14 = [-0.01831, -0.09876]
    zy15 = [0.06558, -0.07878]
    return np.array([zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11, zy12, zy13, zy14, zy15])

sect = SP.section(Nstiffeners=15, Ha=0.205, Ca=0.605, tskin=0.0011, tspar=0.0028,
                    hstiff=0.016, tstiff=0.0012, wstiff=0.019, booms=validation_booms(), remove_booms=True)

SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)

_, _, J = TS.torsionalstiffness(sect)

load = I.get_load(C_a=0.605, l_a=2.661, n_span=150, do_plot=False, fixed_load=5540)

MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy, G_i=27.1e9, E_i=73.1e9, 
            La_i=2.661, x1_i=0.172, x2_i=1.211, x3_i=2.591, d1_i=0.01154, 
            d3_i=0.0184, xa_i=0.35, theta_i=np.radians(28), P_i=97.4e3)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system(power=0)
MC.do_plots()
path = 'Validation/nodes.txt'
file = open(path, "r")
nodes = np.genfromtxt(path, delimiter=",", skip_header=0 )
file.close()

path = 'Validation/elements.txt'
file = open(path, "r")
elements = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()

path = 'Validation/Jambent_displ.csv'
file = open(path, "r")
displ_dat = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()


# find all nodes that lie on the hingeline (0,0) in zy plane
hingeline = np.zeros((np.size(np.unique(nodes[:,1])),5))
original = np.zeros((np.size(np.unique(nodes[:,1])),3))

k=0
for i in range(np.shape(nodes)[0]):
    if nodes[i,2] == 0 and nodes[i,3] == 0:
        # first value in hingeline array are the labels of nodes
        hingeline[k,0] = nodes[i,0]
        # Second value is the original position along the x-axis
        hingeline[k,1] = nodes[i,1]
        k += 1


# get out the displacement per node
# the 5th element in hingeline row is the total magnitude of the displacement

for i in range(np.shape(hingeline)[0]):
    hingeline[i,4] = displ_dat[int(hingeline[i,0]),1]
    hingeline[i,2] = displ_dat[int(hingeline[i,0]),3]
    hingeline[i,3] = displ_dat[int(hingeline[i,0]),4]



x_coords = hingeline[:,1]/1000


# change units to m instead of mm

disp_num = np.zeros((np.shape(x_coords)[0],3))
disp_num[:,0] = x_coords

# compute v = y w = z displacement
for i in range(np.shape(x_coords)[0]):
    disp_num[i,1] = MC.v(x_coords[i])
    disp_num[i,2] = MC.w(x_coords[i])

hingeline = hingeline.transpose()
disp_num = disp_num.transpose()


# compute total Magnitude of displ of the numerical model

mag_num = ((disp_num[1])**2 + (disp_num[2]**2))**0.5

# setup MSE

MSE = 1/np.shape(disp_num)[0] * sum((mag_num-hingeline[4]/1000)**2)
print("MSE = ",MSE)




fig = plt.figure()

ax = Axes3D(fig)
ax.scatter(hingeline[1]/1000,-1*hingeline[2]/1000,-1*hingeline[3]/1000, label="Validation model")
ax.scatter(disp_num[0],-1*disp_num[1],-1*disp_num[2], label="Numerical model")
ax.set_xlabel('x [m]')
ax.set_ylabel('y-displacement [m]')
ax.set_zlabel('z-displacement [m]')
plt.legend()
plt.show()

plt.plot(hingeline[1]/1000,-1*hingeline[2]/1000, label="Validation model", marker="x", linewidth=0)
plt.plot(disp_num[0],-1*disp_num[1], label="Numerical model", marker="x", linewidth=0)
plt.xlabel("x [m]")
plt.ylabel("y-displacement [m]")
plt.legend()
plt.show()

plt.plot(hingeline[1]/1000,-1*hingeline[3]/1000, label="Validation model", marker="x", linewidth=0)
plt.plot(disp_num[0],-1*disp_num[2], label="Numerical model", marker="x", linewidth=0)
plt.xlabel("x [m]")
plt.ylabel("z-displacement [m]")
plt.legend()
plt.show()