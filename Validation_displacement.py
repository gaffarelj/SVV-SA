import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def scatter3d(x,y,z, cs):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlim3d(0, 2500)
    ax.set_ylim3d(-1250, 1250)
    ax.set_zlim3d(-1000, 1000)
    pl = ax.scatter(x, y, z, c=cs, cmap='coolwarm')
    fig.colorbar(pl)

    plt.show()

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
hingeline = np.zeros((np.size(np.unique(nodes[:,1])),5))
original = np.zeros((np.size(np.unique(nodes[:,1])),3))
k=0
for i in range(np.shape(nodes)[0]):
    if nodes[i,2] == 0 and nodes[i,3] == 0:
        hingeline[k,0] = nodes[i,0]
        hingeline[k,1] = nodes[i,1]
        original[k, 0] = nodes[i, 1]
        original[k, 1:] = 0
        k += 1

for i in range(np.shape(hingeline)[0]):
    hingeline[i,2:] = displ_dat[int(hingeline[i,0]),2:]


hingeline = hingeline.transpose()
original = original.transpose()
#plt.plot(hingeline[2],hingeline[3],hingeline[4])

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = Axes3D(fig)
ax.scatter(hingeline[2],hingeline[3],hingeline[4])
ax.scatter(original[0],original[1],original[2])
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()