import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cmx
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable( cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.show()
'''
def opentext(filename):
    lijst = [filename, '.txt']
    filename = ''.join(lijst)
    document = open(filename, 'r')
    text = document.read()
    document.close()
    return text

def cut_text(text, cutsym):
    listperline = text.split(cutsym)
    arrayperline = np.array(listperline)
    return arrayperline
def plotxyz(array):
    return
text1 = opentext('data/nodes')
nodeline = cut_text(text1, '\n')
nodeline0 = nodeline[0]


allnodes = []
for i in range(len(nodeline)):
    nodeline[i] = ''.join(nodeline[i].split())
    nodelist = cut_text(nodeline[i], ',')
    for i in range(len(nodelist)):
        nodelist[i] = float(nodelist[i])
    allnodes.append(nodelist)

text1 = opentext('data/elements')
nodeline = cut_text(text1, '\n')
nodeline0 = nodeline[0]
allelements = []
for i in range(len(nodeline)):
    nodeline[i] = ''.join(nodeline[i].split())
    nodelist = cut_text(nodeline[i], ',')
    for i in range(len(nodelist)):
        nodelist[i] = float(nodelist[i])
    allelements.append(nodelist)

map(float, allelements)

#print('allnodes', allnodes)
#print('allelements', allelements)
'''
# starting to couple every element with 4 points.

path = 'data/nodes.txt'
file = open(path, "r")
nodes = np.genfromtxt(path, delimiter=",", skip_header=0 )
file.close()

path = 'data/elements.txt'
file = open(path, "r")
elements = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()

path = 'data/Bending_result_dat.csv'
file = open(path, "r")
data_frame = np.genfromtxt(path, delimiter=",", skip_header=15)
file.close()

# define a collecting array for all result data
data = np.zeros((np.size(elements),5))


#data[:,0] = np.sum(nodes[np.where(elements[1])])/4

#print(elements[0])
#currentcoords = elements[0,1:]
#print(currentcoords)
#print(nodes[int(currentcoord)-1])

#average_x = (nodes[int(currentcoords[0])-1,1] + nodes[int(currentcoords[1])-1,1] + nodes[int(currentcoords[2])-1,1] + nodes[int(currentcoords[3])-1,1])/4
#print(average_x)

#nested loop finds the averages (integration points) of all elements and stores them in data array
# at the first 3 positions, the last 2 positions are reserved for the average stress and shear stress from
# top to bottom

for j in range(np.shape(elements)[0]):

    currentcoords = elements[j,1:]

    summationx = 0
    summationy = 0
    summationz = 0

    for i in range(4):
        summationx += nodes[int(currentcoords[i])-1,1]
        summationy += nodes[int(currentcoords[i])-1,2]
        summationz += nodes[int(currentcoords[i])-1,3]

    data[j,0] = summationx/4
    data[j,1] = summationy/4
    data[j,2] = summationz/4
    # add the von misses & shear stresses to the data_frame
    data[j, 3] = (data_frame[j, 2] + data_frame[j, 3]) / 2
    data[j, 4] = (data_frame[j, 4] + data_frame[j, 5]) / 2



data = data.transpose()


scatter3d(data[0],data[1],data[2],data[3])

'''
elementarray = []
print(len(allelements))
for i in range(len(allelements)):
    print('elementnumber', i)

    # first node: find in the elementlist the right index and find its corresponding nodes.
    print('1', allnodes[int(allelements[i][1]) - 1])
    # second node
    print('2', allnodes[int(allelements[i][2]) - 1])
    # third node
    print('3', allnodes[int(allelements[i][3]) - 1])
    # fourth node
    print('4', allnodes[int(allelements[i][4]) - 1])

    # find values for element
    elementlist = [allnodes[int(allelements[i][1]) - 1], allnodes[int(allelements[i][2]) - 1],
                   allnodes[int(allelements[i][3]) - 1], allnodes[int(allelements[i][4]) - 1]]
    print(elementlist)
    elementarray.append(elementlist)
# print (elementarray)


'''