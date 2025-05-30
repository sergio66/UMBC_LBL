C
C  PROGRAM        brandnew_02ctt  BLOCK DATA
C
C  PURPOSE        OXYGEN CONTINUUM VALUES FOR COLLISION INDUCED 
C                 ABSORPTION BAND
C
C  DESCRIPTION    VALUES TAKEN LBLRTMv5.10
C
C***********************************************************************
cccccc this is for the main InfraRed portion (1340 - 1850 cm-1)
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2  
c see Applied Optics vol 36, pg 563-567 (1997)

      BLOCK DATA bo2f                                                    

      INTEGER NPTS
      REAL*8 V1S,V2S,DVS,
     *          o0001(50),o0051(50),o0101(03), 
     *          ot0001(50),ot0051(50),ot0101(03) 

      COMMON /o2_f  / V1S,V2S,DVS, NPTS,
     *          o0001,o0051,o0101, 
     *          ot0001,ot0051,ot0101 
                                                                         
      DATA V1S,V2S,DVS,NPTS /1340.000,1850.000,   5.000,  103/               
                                                                         
      DATA o0001/                                                        
     *      0.000E+00,  9.744E-09,  2.256E-08,  3.538E-08,  4.820E-08,   
     *      6.100E-08,  7.400E-08,  8.400E-08,  9.600E-08,  1.200E-07,   
     *      1.620E-07,  2.080E-07,  2.460E-07,  2.850E-07,  3.140E-07,   
     *      3.800E-07,  4.440E-07,  5.000E-07,  5.710E-07,  6.730E-07,   
     *      7.680E-07,  8.530E-07,  9.660E-07,  1.100E-06,  1.210E-06,   
     *      1.330E-06,  1.470E-06,  1.590E-06,  1.690E-06,  1.800E-06,   
     *      1.920E-06,  2.040E-06,  2.150E-06,  2.260E-06,  2.370E-06,   
     *      2.510E-06,  2.670E-06,  2.850E-06,  3.070E-06,  3.420E-06,   
     *      3.830E-06,  4.200E-06,  4.450E-06,  4.600E-06,  4.530E-06,   
     *      4.280E-06,  3.960E-06,  3.680E-06,  3.480E-06,  3.350E-06/   
      DATA o0051/                                                        
     *      3.290E-06,  3.250E-06,  3.230E-06,  3.230E-06,  3.210E-06,   
     *      3.190E-06,  3.110E-06,  3.030E-06,  2.910E-06,  2.800E-06,   
     *      2.650E-06,  2.510E-06,  2.320E-06,  2.130E-06,  1.930E-06,   
     *      1.760E-06,  1.590E-06,  1.420E-06,  1.250E-06,  1.110E-06,   
     *      9.900E-07,  8.880E-07,  7.910E-07,  6.780E-07,  5.870E-07,   
     *      5.240E-07,  4.640E-07,  4.030E-07,  3.570E-07,  3.200E-07,   
     *      2.900E-07,  2.670E-07,  2.420E-07,  2.150E-07,  1.820E-07,   
     *      1.600E-07,  1.460E-07,  1.280E-07,  1.030E-07,  8.700E-08,   
     *      8.100E-08,  7.100E-08,  6.400E-08,  5.807E-08,  5.139E-08,   
     *      4.496E-08,  3.854E-08,  3.212E-08,  2.569E-08,  1.927E-08/   
      DATA o0101/                                                        
     *      1.285E-08,  6.423E-09,  0.000E+00/                           
 
      DATA ot0001/                                                        
     *      4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,   
     *      4.670E+02,  4.000E+02,  3.150E+02,  3.790E+02,  3.680E+02,   
     *      4.750E+02,  5.210E+02,  5.310E+02,  5.120E+02,  4.420E+02,   
     *      4.440E+02,  4.300E+02,  3.810E+02,  3.350E+02,  3.240E+02,   
     *      2.960E+02,  2.480E+02,  2.150E+02,  1.930E+02,  1.580E+02,   
     *      1.270E+02,  1.010E+02,  7.100E+01,  3.100E+01, -6.000E+00,   
     *     -2.600E+01, -4.700E+01, -6.300E+01, -7.900E+01, -8.800E+01,   
     *     -8.800E+01, -8.700E+01, -9.000E+01, -9.800E+01, -9.900E+01,   
     *     -1.090E+02, -1.340E+02, -1.600E+02, -1.670E+02, -1.640E+02,   
     *     -1.580E+02, -1.530E+02, -1.510E+02, -1.560E+02, -1.660E+02/   
      DATA ot0051/                                                        
     *     -1.680E+02, -1.730E+02, -1.700E+02, -1.610E+02, -1.450E+02,   
     *     -1.260E+02, -1.080E+02, -8.400E+01, -5.900E+01, -2.900E+01,   
     *      4.000E+00,  4.100E+01,  7.300E+01,  9.700E+01,  1.230E+02,   
     *      1.590E+02,  1.980E+02,  2.200E+02,  2.420E+02,  2.560E+02,   
     *      2.810E+02,  3.110E+02,  3.340E+02,  3.190E+02,  3.130E+02,   
     *      3.210E+02,  3.230E+02,  3.100E+02,  3.150E+02,  3.200E+02,   
     *      3.350E+02,  3.610E+02,  3.780E+02,  3.730E+02,  3.380E+02,   
     *      3.190E+02,  3.460E+02,  3.220E+02,  2.910E+02,  2.900E+02,   
     *      3.500E+02,  3.710E+02,  5.040E+02,  4.000E+02,  4.000E+02,   
     *      4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02/   
      DATA ot0101/                                                        
     *      4.000E+02,  4.000E+02,  4.000E+02/                           
 
      END 
 
C************************************************************************
