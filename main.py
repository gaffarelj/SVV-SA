import sectionproperties as SP


Am, z_centroid, Iyy, Izz, tskin, boomcoords, boomcoords_hinge, stiff_area = SP.get_geometry()
#print(Am, z_centroid, Iyy, Izz, tskin, stiff_area, boomcoords, "\n\n", boomcoords_hinge, stiff_area)
print(Izz, stiff_area)