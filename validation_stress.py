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
# loading data
path = 'data/nodes.txt'
file = open(path, "r")
nodes = np.genfromtxt(path, delimiter=",", skip_header=0 )
file.close()

path = 'data/elements.txt'
file = open(path, "r")
elements = np.genfromtxt(path, delimiter=",", skip_header=0)
file.close()

case = 2
if case == 1:
    # Bending, no loading
    load = I.get_load(C_a=0.605, l_a=2.661, n_span=150, do_plot=False, fixed_load=0)
    P = 0
    fname = "Bending"
    d1, d3 = 0.01154, 0.0184
elif case == 2:
    # Loading, and bending
    load = I.get_load(C_a=0.605, l_a=2.661, n_span=150, do_plot=False, fixed_load=5540)
    P = 97.4e3
    fname = "Jambent"
    d1, d3 = 0.01154, 0.0184
else:
    # Loading, no bending
    load = I.get_load(C_a=0.605, l_a=2.661, n_span=150, do_plot=False, fixed_load=5540)
    P = 97.4e3
    fname = "Jamstraight"
    d1, d3 = 0, 0

    
path = f'data/{fname}_result.csv'
    

sect = SP.section(Nstiffeners=15, Ha=0.205, Ca=0.605, tskin=0.0011, tspar=0.0028,
                    hstiff=0.016, tstiff=0.0012, wstiff=0.019, booms=validation_booms(), remove_booms=True)
SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
_, _, J = TS.torsionalstiffness(sect)

MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy, G_i=28e9, E_i=73.1e9, 
            La_i=2.661, x1_i=0.172, x2_i=1.211, x3_i=2.591, d1_i=d1, 
            d3_i=d3, xa_i=0.35, theta_i=np.radians(28), P_i=P)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system(power=0, power_t=2)

file = open(path, "r")
data_frame = np.genfromtxt(file.readlines(), delimiter=',', dtype=np.float32)
file.close()

path = 'data/Jambent_result.csv'
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

# Making data frame
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
    data[j, 3] = (data_frame[j, 2] + data_frame[j, 3]) / 2
    data[j, 4] = (data_frame[j, 4] + data_frame[j, 5]) / 2

# plotting  slices with unknown spacing
next_section = np.unique(data[:, 0])

# select a cross section place to monitor the stresses

discr_miss =[]
discr_shear = []
xloc = []

lst_x = sorted(list(next_section))

# get stress for given x,y,z
numericaldata = {}  # x:[z,y,vm,ss]
for x in lst_x:  # x is mm
    print("x", x, end="\r")
    x = x/1000
    s = STR.stress(MC.Mz(x), MC.My(x), MC.Sz(x), MC.Sy(x), MC.T(x), sect, q1, q2, q3, q4, q5, q6, show_plot=False)
    s.section_stress()
    
    vm = np.array(s.vm_stresses)
    ss = np.array(s.shear_stress)
    b = ss[:, -1]
    b = b.reshape(620,1)
    vmss = np.concatenate((vm, b), axis=1)
    numericaldata[x] = vmss
print()

def get_stresses(x,y,z):  # all inputs are mm
    vmss = numericaldata[x/1000]
    distance = np.sum((vmss[:, 0:2] - np.array([z/1000, y/1000])) ** 2, axis=1)
    indice = np.argmin(distance)
    #print(x, y, z, vmss[indice, 2:])
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
        local_mae_miss = sum(np.fabs(section_data[2]*10**6/1e9 - section_data[4]/1e9))/len(section_data[2])
        local_mae_shear = sum(np.fabs(section_data[3]*10**6/1e9 - section_data[5]/1e9))/len(section_data[3])

        xloc.append(round(x,6))
        # print(local_mse_miss,local_mse_shear)
        discr_miss.append(local_mae_miss)
        discr_shear.append(local_mae_shear)
    else:
        incompletex.append(x)

# plot list of error
print(np.nanmean(discr_miss),np.nanmean(discr_shear))
plt.scatter([x/1000 for x in xloc], discr_miss, label="von Mises stress", marker="x", s=5)
plt.scatter([x/1000 for x in xloc], discr_shear, label="Shear stress", marker="x", s=5)
plt.legend()
plt.xlabel("x [m]")
plt.ylabel("stress [GPa]")
plt.savefig(f"plots/validation/stress_{case}.pdf", bbox_inches='tight')