''' program build by Mees

Program aims to validate the numerical model.

First the validation model will be loaded, merging the input and output data.
This creates an multi dimensional array containing the location of the nodes with
the corresponding magnitudes of the von misses stress and shear stress.

'''

#load modules
import numpy as np


# load nodes for loadcase 1

path = 'data/Elementstest.txt'
file = open(path, "r")
elements = np.genfromtxt(path, delimiter=",", skip_header=1, comments="*")
file.close()

path = 'data/Nodestest.txt'
file = open(path, "r")
nodes = np.genfromtxt(path, delimiter=",", skip_header=1, comments="*")
file.close()

path = 'data/results.txt'
file = open(path, "r")
output = np.genfromtxt(path, delimiter=",", skip_header=3, comments="*")
file.close()

print(output)


#for i in range(np.shape(elements)[0]):
