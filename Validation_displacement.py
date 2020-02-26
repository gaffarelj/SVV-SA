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

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)

SC.set_sect(sect)
#qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
xi = -0.007513567161803937

_, _, J = TS.tosionalstiffness(sect)

MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()
MC.do_plots()
MC.plot_result(MC.My, "My_b")

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

hingeline = hingeline
disp_num = np.zeros((np.shape(x_coords)[0],3))
disp_num[:,0] = x_coords

# compute v = y w = z displacement
for i in range(np.shape(x_coords)[0]):
    disp_num[i,1] = MC.v(x_coords[i])
    disp_num[i,2] = MC.w(x_coords[i])

hingeline = hingeline.transpose()
disp_num = disp_num.transpose()
# compute total Magnitude of displ
print(hingeline[4])
mag_num = ((disp_num[1])**2 + (disp_num[2]**2))**0.5

# setup MSE


MSE = 1/np.shape(disp_num)[0] * sum((mag_num-hingeline[4]/1000)**2)
print(MSE)




fig = plt.figure()

ax = Axes3D(fig)
ax.scatter(hingeline[1],hingeline[3],hingeline[2])
ax.scatter(disp_num[0],disp_num[1],disp_num[2])
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
#ax.set_xlim3d(0, 2500)
#ax.set_ylim3d(-10,10)
#ax.set_zlim3d(-10, 10)
plt.show()
