import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import  torsionstiffness as TS
import stress as STR
import max_stress as MS
import numpy as np
import math
sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
_, _, J = TS.torsionalstiffness(sect)
MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()

# CHECKING BOUNDARY CONDITIONS
