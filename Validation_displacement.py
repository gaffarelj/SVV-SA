import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

path = 'data/nodes.txt'
file = open(path, "r")
nodes = np.genfromtxt(path, delimiter=",", skip_header=0 )
file.close()

path = 'data/elements.txt'
file = open(path, "r")
elements = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()

path = 'data/Bending_result_notext_dat.csv'
file = open(path, "r")
data_frame = np.genfromtxt(path, delimiter=";", skip_header=0)
file.close()

path = 'data/Displ_jamstraight.csv'
file = open(path, "r")
displ_dat = np.genfromtxt(path, delimiter=",", skip_header=3)
file.close()


# find all nodes that lie on the hingeline (0,0) in zy plane
hingeline = np.zeros((np.size(np.unique(nodes[:,1])),4))
original = np.zeros((np.size(np.unique(nodes[:,1])),3))

k=0
for i in range(np.shape(nodes)[0]):
    if nodes[i,2] == 0 and nodes[i,3] == 0:
        # first value in hingeline array are the labels of nodes
        hingeline[k,0] = nodes[i,0]
        # Second value is the original position along the x-axis
        hingeline[k,1] = nodes[i,1]
        # creates copy to compare against (Original)
        original[k, 0] = nodes[i,1]
        original[k, 1:] = 0
        k += 1

for i in range(np.shape(hingeline)[0]):
    hingeline[i,1] = hingeline[i,1] + displ_dat[int(hingeline[i,0]),2]
    hingeline[i,2] = displ_dat[int(hingeline[i,0]),3]
    hingeline[i,3] = displ_dat[int(hingeline[i, 0]),4]


hingeline = hingeline.transpose()
original = original.transpose()
# plotting the result

fig = plt.figure()

ax = Axes3D(fig)
ax.scatter(hingeline[1],hingeline[2],hingeline[3])
ax.scatter(original[0],original[1],original[2])
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_xlim3d(0, 2500)
ax.set_ylim3d(-10,10)
ax.set_zlim3d(-10, 10)
plt.show()
