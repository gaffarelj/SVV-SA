import numpy as np
"""
Section Properties - idealized aileron

Computes the centroid, enclosed area and moments of inertia (Iyy and Izz) for the idealized section.
All units are meters unless otherwise specified.
Stiffener positions are hardcoded (frame of reference is the hinge line), more can be added.
Inputs are below, as given in Aircraft Data, coded for CRJ700.
"""

class section():
    def __init__(self, Nstiffeners=13, Ha=0.173, Ca=0.484, tskin=0.0011, tspar=0.0025,
                    hstiff=0.014, tstiff=0.0012, wstiff=0.018, booms=None):
        self.Nstiffeners = Nstiffeners  # This is actually hard coded by get_boomcoords
        self.Ha = Ha
        self.Ca = Ca
        self.tskin = tskin
        self.tspar = tspar
        self.hstiff = hstiff
        self.tstiff = tstiff
        self.wstiff = wstiff
        self.r = Ha/2
        self.l_spar_to_end = Ca - self.r
        self.beta = np.arctan(self.r / self.l_spar_to_end)
        self.l_topskin = self.l_spar_to_end / np.cos(self.beta)
        # Placeholder vars, to be computed
        self.stiff_area, self.Iyy, self.Izz, self.boomcoords, self.boomcoords_hinge, \
            self.Am, self.z_centroid = 0, 0, 0, [], [], 0, 0
        if booms is None:
            self.get_boomcoords()
        else:
            self.boomcoords = booms
        self.boomcoords_hinge = self.boomcoords
        # stiffener area (in meters^2)
        self.stiff_area = (self.hstiff + self.wstiff) * self.tstiff     # small angle approximation is used
        # enclosed area (in meters^2)
        self.Am = np.pi * self.r ** 2 + self.l_spar_to_end * self.r
        self.comp_centroid()
        # moments of inertia (Iyy and Izz about the centroid, in meters^4)
        self.comp_Iyy()
        self.comp_Izz()

    def get_boomcoords(self):
        # stiffener positions (origin at hinge line, from bottom right going clockwise, in meters)
        # [z, y]
        zy1 = [-0.35671, -0.00888]
        zy2 = [-0.27513, -0.02663]
        zy3 = [-0.19355, -0.04438]
        zy4 = [-0.11197, -0.06213]
        zy5 = [-0.03039, -0.07989]
        zy6 = [0.04924, -0.07112]
        zy7 = [0.0865, 0]
        zy8 = [0.04924, 0.07112]
        zy9 = [-0.03039, 0.07989]
        zy10 = [-0.11197, 0.06213]
        zy11 = [-0.19355, 0.0444]
        zy12 = [-0.27513, 0.02663]
        zy13 = [-0.35671, 0.00888]
        self.boomcoords = np.array([zy1, zy2, zy3, zy4, zy5, zy6, zy7, zy8, zy9, zy10, zy11, zy12, zy13])

    def comp_centroid(self):
        # centroid (z-coordinate from hinge line, in meters)
        sumA_skin = (0.5 * np.pi * self.Ha * self.tskin) + 2 * (self.l_topskin * self.tskin) + (self.tspar * self.Ha)
        sumAz_skin = (0.5 * np.pi * self.Ha * self.tskin) * (self.Ha / np.pi) + \
            2 * (self.l_topskin * self.tskin) * (-self.l_spar_to_end/2)
        sumAz_boom = 0
        for i in range(len(self.boomcoords)):
            sumAz_boom += self.boomcoords[i, 0] * self.stiff_area
        self.z_centroid = (sumAz_skin + sumAz_boom) / (sumA_skin + len(self.boomcoords) * self.stiff_area)
        boomcoords_centroid = np.zeros((len(self.boomcoords), 2))
        for i in range(len(self.boomcoords)):        # shift boom origin from hinge line to centroid
            boomcoords_centroid[i, 0] = self.boomcoords[i, 0] - self.z_centroid
            boomcoords_centroid[i, 1] = self.boomcoords[i, 1]
        self.boomcoords = boomcoords_centroid

    def comp_Iyy(self):
        # Iyy
        Iyy_thin = self.tskin * self.l_topskin ** 3 * (np.cos(self.beta)) ** 2 / 12 + \
            self.tskin * self.l_topskin * (self.l_spar_to_end/2 + self.z_centroid) ** 2
        Iyy_spar = self.Ha * self.tspar ** 3 / 12 + self.tspar * self.Ha * self.z_centroid ** 2
        Iyy_ring = self.tskin * self.r ** 3 * (np.pi ** 2 - 8) / (2 * np.pi) + \
            0.5 * np.pi * self.Ha * self.tskin * (self.Ha / np.pi - self.z_centroid) ** 2
        Iyy_booms = 0
        for i in range(len(self.boomcoords)):
            Iyy_booms += ((self.boomcoords[i, 0]) ** 2) * self.stiff_area
        self.Iyy = 2 * Iyy_thin + Iyy_spar + Iyy_ring + Iyy_booms

    def comp_Izz(self):
        # Izz
        Izz_thin = self.tskin * self.l_topskin ** 3 * np.sin(self.beta) ** 2 / 12 + \
            self.tskin * self.l_topskin * (self.l_topskin/2 * np.sin(self.beta)) ** 2
        Izz_spar = self.tspar * self.Ha ** 3 / 12
        Izz_ring = self.r ** 3 * self.tskin * np.pi / 2
        Izz_booms = 0
        for i in range(len(self.boomcoords)):
            Izz_booms += (self.boomcoords[i, 1]) ** 2 * self.stiff_area
        self.Izz = 2 * Izz_thin + Izz_spar + Izz_ring + Izz_booms