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

path = 'data/Jamstraight_result_notext_nosmallvalues_dat2.csv'
file = open(path, "r")
data_frame = np.genfromtxt(path, delimiter=";", skip_header=0)
file.close()

path = 'data/Displ_jamstraight.csv'
file = open(path, "r")
displ_dat = np.genfromtxt(path, delimiter=",", skip_header=3)
file.close()

hingeline = np.zeros()
