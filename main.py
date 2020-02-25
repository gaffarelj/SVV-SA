import sectionproperties as SP
import shearcentre as SC
import macaulay as MC
import  torsionstiffness as TS

section_prop = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
# The following parameters can be used out of section_prop (example: print(section_prop.z_centroid)
# Iyy, Izz, Am, z_centroid, stiff_area, boomcoords, boomcoords_hinge,
# tskin, tspar, hstiff, tstiff, wstiff, Ha, Ca, r, l_topskin, beta, l_spar_to_end

SC.set_sect(section_prop)
qsI, qsII, q1, q2, q3, q4, q5, q6, xi = SC.shear_centre(1000)

_, _, J = TS.tosionalstiffness()

MC.set_vars(xi, J, section_prop.r, section_prop.Izz, section_prop.Iyy)
Ry1, Ry2, Ry3, Rz1, Rz2, Rz3, Fa, C1, C2, C3, C4, C5 = MC.system()
MC.do_plots()
MC.plot_result(MC.My, "My_b")


