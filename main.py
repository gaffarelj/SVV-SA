import sectionproperties as SP


section_prop = SP.section(Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025, hstiff=0.014, tstiff=0.0012, wstiff=0.018)
# The following parameters can be used out of section_prop (example: print(section_prop.z_centroid)
# Iyy, Izz, Am, z_centroid, stiff_area, boomcoords, boomcoords_hinge,
# tskin, tspar, hstiff, tstiff, wstiff, Ha, Ca, r, l_topskin, beta, l_spar_to_end

print(section_prop.Izz)
print(section_prop.Iyy)
print(section_prop.Am)
print(section_prop.z_centroid)
print(section_prop.tskin)
print(section_prop.stiff_area)
print(section_prop.boomcoords)
print()
print(section_prop.boomcoords_hinge)