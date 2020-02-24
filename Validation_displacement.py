import numpy as np
from matplotlib import pyplot as plt

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



hingeline = np.zeros((np.size(np.unique(nodes[:,1])),2))
k=0
for i in range(np.shape(nodes)[0]):
    if nodes[i,2] == 0 and nodes[i,3] == 0:
        hingeline[k,0] = nodes[i,0]
        hingeline[k,1] = nodes[i,1]
        k += 1

print(np.shape(hingeline))
print(np.shape(np.unique(nodes[:,1])))