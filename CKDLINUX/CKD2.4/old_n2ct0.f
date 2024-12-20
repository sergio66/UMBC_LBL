C***********************************************************************
C
C  PROGRAM        new_N2CT0  BLOCK DATA
C
C  PURPOSE        NITROGEN CONTINUUM VALUES FOR COLLISION INDUCED 
C                 ABSORPTION BAND
C
C  VERSION        1.0 Scott Hannon November 1993
C                 This is a replacement for N2CT0
C
C  DESCRIPTION    VALUES TAKEN FROM
C                 Applied Optics: Collision-induced Absorption in the
C                 Fundamental Band of N2. V. Menoux, R. Le Doucen,
C                 C. Boulet, A. Roblin, and A.M. Bouchardy
C                 20 January 1993, Vol 32, no3
C                 The data for xxx5 cm-1 has been faked using a
C                 quadratic spline (except near the ends) for the
C                 regions 2100-2300 cm-1 and 2370-2640 cm-1, as the
C                 data in these two regions were only given at a
C                 xx10 cm-1 point spacing.
C                 units of (10^-6 cm-1 amagat-2)
C
C***********************************************************************
C
       BLOCK DATA N2CT0
C-----------------------------------------------------------------------
       REAL*8 VB22,VT22,DV22,NPT22,CN2297(111),CN2273(111),
     + CN2253(111),CN2233(111),CN2213(111),CN2193(111)

       COMMON /CN2C0/ VB22,VT22,DV22,NPT22,CN2297,CN2273,
     + CN2253,CN2233,CN2213,CN2193
C-----------------------------------------------------------------------
C
       DATA VB22,VT22,DV22,NPT22 / 2095.0,2645.0,5.0,111 /
C
       DATA CN2297/ 0.0,
     +  0.015, 0.017, 0.020, 0.024, 0.030, 0.038, 0.050, 0.065,
     +  0.080, 0.093, 0.110, 0.134, 0.160, 0.183, 0.205, 0.231,
     +  0.260, 0.290, 0.320, 0.352, 0.385, 0.419, 0.455, 0.493,
     +  0.535, 0.581, 0.625, 0.664, 0.700, 0.737, 0.775, 0.811,
     +  0.845, 0.875, 0.900, 0.918, 0.940, 0.972, 1.000, 1.020,
     +  1.100, 1.280, 1.450, 1.600, 1.800, 1.930, 2.000, 1.940,
     +  1.820, 1.650, 1.500, 1.430, 1.370, 1.350, 1.345, 1.346,
     +  1.350, 1.351, 1.350, 1.348, 1.340, 1.323, 1.300, 1.276,
     +  1.245, 1.200, 1.150, 1.102, 1.055, 1.005, 0.945, 0.873,
     +  0.800, 0.738, 0.685, 0.636, 0.590, 0.548, 0.505, 0.459,
     +  0.415, 0.376, 0.345, 0.320, 0.295, 0.266, 0.235, 0.208,
     +  0.185, 0.166, 0.150, 0.137, 0.125, 0.112, 0.100, 0.092,
     +  0.085, 0.078, 0.070, 0.060, 0.050, 0.041, 0.035, 0.032,
     +  0.030, 0.025, 0.020, 0.016, 0.015, 0.0/
C
       DATA CN2273/ 0.0,
     +  0.010, 0.013, 0.015, 0.017, 0.020, 0.023, 0.030, 0.043,
     +  0.060, 0.080, 0.100, 0.120, 0.140, 0.159, 0.180, 0.203,
     +  0.230, 0.260, 0.290, 0.318, 0.350, 0.389, 0.430, 0.469,
     +  0.505, 0.541, 0.580, 0.627, 0.680, 0.735, 0.785, 0.825,
     +  0.860, 0.893, 0.925, 0.953, 0.975, 0.994, 1.025, 1.084,
     +  1.185, 1.330, 1.490, 1.675, 1.880, 2.040, 2.100, 2.055,
     +  1.915, 1.755, 1.600, 1.510, 1.450, 1.420, 1.405, 1.400,
     +  1.400, 1.401, 1.400, 1.396, 1.390, 1.382, 1.365, 1.334,
     +  1.290, 1.237, 1.180, 1.122, 1.060, 0.992, 0.920, 0.849,
     +  0.775, 0.697, 0.630, 0.585, 0.550, 0.509, 0.465, 0.423,
     +  0.385, 0.353, 0.325, 0.299, 0.275, 0.254, 0.230, 0.199,
     +  0.170, 0.153, 0.140, 0.125, 0.110, 0.099, 0.090, 0.079,
     +  0.070, 0.064, 0.060, 0.054, 0.045, 0.036, 0.030, 0.030,
     +  0.030, 0.024, 0.015, 0.012, 0.010, 0.0/
C
       DATA CN2253/ 0.0,
     +  0.010, 0.010, 0.010, 0.013, 0.015, 0.016, 0.020, 0.032,
     +  0.050, 0.070, 0.090, 0.107, 0.125, 0.146, 0.170, 0.194,
     +  0.220, 0.248, 0.275, 0.301, 0.330, 0.365, 0.400, 0.429,
     +  0.460, 0.502, 0.550, 0.596, 0.650, 0.717, 0.780, 0.824,
     +  0.860, 0.901, 0.940, 0.967, 0.990, 1.020, 1.060, 1.118,
     +  1.220, 1.370, 1.520, 1.750, 1.975, 2.130, 2.200, 2.165,
     +  2.060, 1.860, 1.665, 1.580, 1.515, 1.490, 1.475, 1.469,
     +  1.470, 1.473, 1.475, 1.476, 1.470, 1.453, 1.425, 1.384,
     +  1.330, 1.265, 1.200, 1.145, 1.080, 0.989, 0.890, 0.805,
     +  0.730, 0.660, 0.600, 0.555, 0.515, 0.471, 0.425, 0.384,
     +  0.350, 0.324, 0.300, 0.271, 0.240, 0.214, 0.190, 0.168,
     +  0.150, 0.140, 0.130, 0.115, 0.100, 0.089, 0.080, 0.069,
     +  0.060, 0.055, 0.050, 0.042, 0.035, 0.032, 0.030, 0.025,
     +  0.020, 0.020, 0.020, 0.018, 0.010, 0.0/
C
       DATA CN2233/ 0.0,
     +  0.010, 0.010, 0.010, 0.013, 0.015, 0.017, 0.020, 0.028,
     +  0.040, 0.057, 0.075, 0.092, 0.110, 0.132, 0.155, 0.176,
     +  0.200, 0.229, 0.260, 0.287, 0.310, 0.331, 0.355, 0.388,
     +  0.430, 0.478, 0.525, 0.569, 0.630, 0.721, 0.800, 0.831,
     +  0.850, 0.895, 0.950, 0.991, 1.025, 1.061, 1.100, 1.149,
     +  1.250, 1.410, 1.560, 1.850, 2.080, 2.285, 2.375, 2.310,
     +  2.165, 2.000, 1.810, 1.700, 1.615, 1.580, 1.570, 1.571,
     +  1.575, 1.574, 1.570, 1.564, 1.550, 1.521, 1.480, 1.432,
     +  1.370, 1.292, 1.220, 1.169, 1.100, 0.979, 0.850, 0.758,
     +  0.690, 0.623, 0.560, 0.509, 0.465, 0.424, 0.385, 0.348,
     +  0.320, 0.305, 0.285, 0.247, 0.210, 0.191, 0.180, 0.161,
     +  0.140, 0.123, 0.110, 0.099, 0.090, 0.080, 0.070, 0.059,
     +  0.050, 0.045, 0.040, 0.032, 0.025, 0.021, 0.020, 0.020,
     +  0.020, 0.016, 0.010, 0.010, 0.010, 0.0/
C
       DATA CN2213/ 0.0,
     +  0.010, 0.010, 0.010, 0.010, 0.010, 0.015, 0.020, 0.023,
     +  0.030, 0.045, 0.060, 0.068, 0.080, 0.104, 0.130, 0.146,
     +  0.160, 0.179, 0.205, 0.235, 0.265, 0.294, 0.325, 0.362,
     +  0.405, 0.453, 0.500, 0.545, 0.610, 0.712, 0.810, 0.866,
     +  0.900, 0.942, 0.990, 1.037, 1.085, 1.139, 1.200, 1.268,
     +  1.355, 1.470, 1.650, 1.960, 2.165, 2.410, 2.540, 2.440,
     +  2.325, 2.140, 1.900, 1.800, 1.750, 1.725, 1.720, 1.720,
     +  1.720, 1.717, 1.710, 1.697, 1.670, 1.625, 1.570, 1.512,
     +  1.440, 1.346, 1.250, 1.170, 1.080, 0.958, 0.830, 0.729,
     +  0.650, 0.583, 0.530, 0.490, 0.450, 0.399, 0.350, 0.319,
     +  0.300, 0.282, 0.260, 0.230, 0.200, 0.178, 0.160, 0.139,
     +  0.120, 0.108, 0.100, 0.092, 0.080, 0.064, 0.050, 0.044,
     +  0.040, 0.033, 0.025, 0.019, 0.015, 0.012, 0.010, 0.010,
     +  0.010, 0.010, 0.010, 0.008, 0.005, 0.0/

       DATA CN2193/ 0.0,
     +  0.005, 0.009, 0.010, 0.010, 0.010, 0.012, 0.015, 0.017,
     +  0.020, 0.024, 0.030, 0.038, 0.050, 0.065, 0.080, 0.095,
     +  0.110, 0.128, 0.150, 0.175, 0.205, 0.239, 0.275, 0.309,
     +  0.350, 0.403, 0.455, 0.502, 0.580, 0.718, 0.860, 0.950,
     +  1.010, 1.071, 1.130, 1.180, 1.230, 1.289, 1.350, 1.409,
     +  1.500, 1.650, 1.850, 2.150, 2.470, 2.640, 2.750, 2.700,
     +  2.580, 2.400, 2.150, 2.025, 1.970, 1.955, 1.960, 1.956,
     +  1.950, 1.950, 1.940, 1.902, 1.850, 1.799, 1.750, 1.694,
     +  1.610, 1.486, 1.350, 1.228, 1.100, 0.945, 0.800, 0.701,
     +  0.630, 0.563, 0.500, 0.445, 0.400, 0.366, 0.335, 0.302,
     +  0.270, 0.243, 0.220, 0.200, 0.180, 0.159, 0.140, 0.125,
     +  0.110, 0.092, 0.075, 0.061, 0.050, 0.041, 0.035, 0.030,
     +  0.025, 0.017, 0.010, 0.010, 0.010, 0.010, 0.010, 0.007,
     +  0.005, 0.005, 0.005, 0.005, 0.002, 0.0/
C
      END
