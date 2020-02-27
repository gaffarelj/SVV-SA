import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from io import BytesIO
# get data from numerical model
import macaulay as MC
import stress as STR
import sectionproperties as SP
import shearcentre as SC
import torsionstiffness as TS

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
_, _, J = TS.torsionalstiffness(sect)
MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()

def scatter3d(x,y,z, cs):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlim3d(0, 2500)
    ax.set_ylim3d(-1250, 1250)
    ax.set_zlim3d(-1000, 1000)
    pl = ax.scatter(x, y, z, c=cs, cmap='coolwarm')
    fig.colorbar(pl)

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

path = 'data/Bending_result_notext_dat.csv'
file = open(path, "r")
lines = ' '.join([s.replace(',', '.') for s in file.readlines()])
data_frame = np.genfromtxt(BytesIO(lines.encode('utf-8')), delimiter=';', dtype=np.float32)
# data_frame = np.genfromtxt(path, delimiter=";", skip_header=0, dtype=float)
# data_frame = pd.read_csv(path,delimiter=";")
file.close()

# for item in data_frame:
#     print(item)

path = 'data/Displ_jamstraight.csv'
file = open(path, "r")
displ_dat = np.genfromtxt(path, delimiter=",", skip_header=3)
file.close()
'''
# optional --> add displacement to all nodes

nodes[:,1:] = nodes[:,1:] + displ_dat[:,2:]
'''

# define a collecting array for all result data
data = np.zeros((elements.shape[0],5))

# nested loop finds the averages (integration points) of all elements and stores them in data array
# at the first 3 positions, the last 2 positions are reserved for the average stress and shear stress from
# top to bottom

############################# Making data frame ##############################################
for j in range(elements.shape[0]):

    currentcoords = elements[j,1:]  # the four nodes
    summationx = 0
    summationy = 0
    summationz = 0

    for i in range(4):
        summationx += nodes[int(currentcoords[i])-1,1]
        summationy += nodes[int(currentcoords[i])-1,2]
        summationz += nodes[int(currentcoords[i])-1,3]

    data[j, 0] = summationx/4
    data[j, 1] = summationy/4
    data[j, 2] = summationz/4  # coord of the element
    # add the von misses & shear stresses to the data_frame
    # print(data_frame[j, 2],data_frame[j, 3],(data_frame[j, 2] + data_frame[j, 3]) / 2)
    data[j, 3] = (data_frame[j, 2] + data_frame[j, 3]) / 2
    # print(data_frame[j, 4] + data_frame[j, 5],(data_frame[j, 4] + data_frame[j, 5]) / 2)
    data[j, 4] = (data_frame[j, 4] + data_frame[j, 5]) / 2

######################################### Plotting in 4D #####################################

#
# data = data.transpose()
#
#
# scatter3d(data[0],data[1],data[2],data[3])

##################################### plotting  slices with unknown spacing ####################

next_section = np.unique(data[:, 0])

# select a cross section place to monitor the stresses

discr_miss =[]
discr_shear = []
xloc = []

lst_x = sorted(list(next_section))

# get stress for given x,y,z
numericaldata = {}  # x:[z,y,vm,ss]
for x in lst_x:  # x is mm
    print(x)
    x = x/1000
    s = STR.stress(MC.Mz(x), MC.My(x), MC.Sz(x), MC.Sy(x), MC.T(x), sect, q1, q2, q3, q4, q5, q6, show_plot=False)
    s.section_stress()
    vm = np.array(s.vm_stresses)
    ss = np.array(s.shear_stress)
    b = ss[:, -1]
    b = b.reshape(522,1)
    vmss = np.concatenate((vm, b), axis=1)
    numericaldata[x] = vmss

def get_stresses(x,y,z):  # all inputs are mm
    vmss = numericaldata[x/1000]
    distance = np.sum((vmss[:, 0:2] - np.array([z/1000, y/1000])) ** 2, axis=1)
    indice = np.argmin(distance)
    print(x, y, z, vmss[indice, 2:])
    return vmss[indice, 2:]


for x in lst_x:
    section_data = np.zeros((62, 6))

sections = {}  # x:[y,z,vm,ss]
for row in data:  # row = [x,y,z,vm,ss]
    # print(row)
    if row[0] in sections:
        sections[row[0]].append(row[1:])
    else:
        sections[row[0]] = [row[1:]]

incompletex = []
for x, info in sections.items():  # info=[[y,z,vm,ss],........]
    if len(info) == 62:
        section_data = np.zeros((62, 6))
        for i, item in enumerate(info):
            section_data[i,0:4] = list(item)
            section_data[i,4:] = get_stresses(x, item[0], item[1])
        section_data = section_data.transpose()
        print(section_data[2:4])
        local_mse_miss = 1 / 62 * sum([i**2 for i in (section_data[2] - section_data[4]*10**6)])
        local_mse_shear = 1 / 62 * sum([i**2 for i in (section_data[3] - section_data[5]*10**6)])

        xloc.append(round(x,6))
        # print(local_mse_miss,local_mse_shear)
        discr_miss.append(local_mse_miss)
        discr_shear.append(local_mse_shear)
    else:
        incompletex.append(x)


'''
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


    # local_mse_miss = 1/62 * sum((section_data[2]-section_data[4]))**2
    # local_mse_shear = 1/62 * sum((section_data[3]-section_data[5]))**2

    local_mse_miss = 1 / 62 * sum([i ** 2 for i in (section_data[2] - section_data[4])])
    local_mse_shear = 1 / 62 * sum([i ** 2 for i in (section_data[3] - section_data[5])])

    # store result in list
    xloc.append(round(next_section[numb],6))
    discr_miss.append(local_mse_miss)
    discr_shear.append(local_mse_shear)
'''


# plot list of error
print(sum(xloc),sum(discr_miss),sum(discr_shear))
plt.plot(xloc, discr_miss, xloc, discr_shear)
plt.show()



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
