import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import  torsionstiffness as TS
import stress as STR

sect = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)

SC.set_sect(sect)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)
#xi = -0.007513567161803937

_, _, J = TS.tosionalstiffness(sect)

#MC.set_vars(xi, J, sect.r, sect.Izz, sect.Iyy)
#Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()
#MC.do_plots()
#MC.plot_result(MC.My, "My_b")


x = 0.554
#Mz, My = MC.Mz(x)[0], MC.My(x)[0]
Mz, My = -34645.28897198964, 119630.64930303364

s = STR.stress(Mz, My, sect, q1, q2, q3, q4, q5, q6)
s.section_stress()
#s.plot_stress()
s.plot_shear_flows()