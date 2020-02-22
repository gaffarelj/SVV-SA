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

# starting to couple every element with 4 points.
######################################## loading data ##########################################
path = 'data/nodes.txt'
file = open(path, "r")
nodes = np.genfromtxt(path, delimiter=",", skip_header=0 )
file.close()

path = 'data/elements.txt'
file = open(path, "r")
elements = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()

path = 'data/Jamstraight_result_notext_nosmallvalues_dat2.csv'
file = open(path, "r")
data_frame = np.genfromtxt(path, delimiter=";", skip_header=0)
file.close()

path = 'data/Displ_jamstraight.csv'
file = open(path, "r")
displ_dat = np.genfromtxt(path, delimiter=",", skip_header=3)
file.close()

# optional --> add displacement to all nodes

nodes[:,1:] = nodes[:,1:] + displ_dat[:,2:]

# define a collecting array for all result data
data = np.zeros((np.size(elements),5))


# nested loop finds the averages (integration points) of all elements and stores them in data array
# at the first 3 positions, the last 2 positions are reserved for the average stress and shear stress from
# top to bottom

############################# Making data frame ##############################################
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

######################################### Plotting in 4D #####################################


data = data.transpose()


scatter3d(data[0],data[1],data[2],data[4])


##################################### plotting  slices with unknown spacing ####################
'''
next_section = np.unique(data[:,0])

# select a cross section place to monitor the stresses

discr_miss =[]
discr_shear = []
xloc = []
for numb in range(np.shape(next_section)[0]):
    section_data = np.zeros((62, 6))

    #print(next_section[0:40])
    # spacing of complete 24.6856995‬mm

    # first check for section completeness (62 points)
    i = 0
    for j in range(np.shape(data)[0]):
        if data[j,0] == next_section[numb]:
            i += 1

    # if section complete --> get out data for particular layer from the data array
    if i == 62:
        k = 0
        for j in range(np.shape(data)[0]):
            if data[j,0] == next_section[numb]:
                section_data[k,0] = data[j,1]
                section_data[k,1] = data[j,2]
                section_data[k,2] = data[j,3]
                section_data[k,3] = data[j,4]
                k += 1


        section_data =section_data.transpose()


        #plt.scatter(section_data[1],section_data[0],c=section_data[2])
        #plt.colorbar()
        #plt.show()


    ############################# setup error testing ################################

    # section data only populated by validation model, the two colums left to be filled by sampling the numerical model.


    local_mse_miss = 1/62 * sum((section_data[2]-section_data[4]))**2
    local_mse_shear = 1/62 * sum((section_data[3]-section_data[5]))**2

    # store result in list
    xloc.append(next_section[numb])
    discr_miss.append(local_mse_miss)
    discr_shear.append(local_mse_shear)


# plot list of error

plt.plot(xloc,discr_miss,xloc,discr_shear)
#plt.show()

'''

























'''
# 3d plot the sections that are incomplete --> to check what is their physical
incomplete = []

for k in range(np.shape(next_section)[0]):
    i = 0
    current_section = []
    for j in range(np.shape(data)[0]):
        if data[j, 0] == next_section[k]:
            i += 1
            current_section.append(data[j,:])
    if i != 62:
        incomplete.append(current_section)

incomplete = np.concatenate(incomplete,axis=0)

print(np.shape(incomplete))

'''
'''
section_data = np.zeros((62,3))


# select a cross section place to monitor the stresses
#print(next_section)
numb = 1

# define starting point of complete section
i = 0
found = False
k = 0
while found == False:
    i = 0
    for j in range(np.shape(data)[0]):
        if data[j,0] == next_section[k]:
            i += 1
    k += 1
    if i == 62:
        found = True

start_int = k-1


# find the discepancies between next_section

dummy = np.zeros(np.shape(next_section))
dummy[1:] = next_section[0:-1]

spacing = next_section-dummy

#print(next_section[start_int]+spacing[start_int+1])
#print(next_section[start_int], next_section[start_int+1])

i = 0
for j in range(np.shape(data)[0]):
    if data[j,0] == next_section[numb]:
        i += 1

# if section complete --> get out data for particular layer from the data array
if i != 62:
    for j in range(np.shape(data)[0]):
        if data[j,0] == next_section[numb]:
            data[j,:] = 0
            k += 1


for j in range(np.shape(data)[0]):
    if data[j,0] >= xloc[loc_int] - 0.2*step and data[j,0] <= xloc[loc_int] + 0.2*step:
        section_data[k,0] = data[j,1]
        section_data[k,1] = data[j,2]
        section_data[k,2] = data[j,3]
        k += 1

section_data =section_data.transpose()

plt.scatter(section_data[1],section_data[0],c=section_data[2])
plt.colorbar()
plt.show()

############################################# working code to plot slicewise but with guesswork for spacing #########

'''
