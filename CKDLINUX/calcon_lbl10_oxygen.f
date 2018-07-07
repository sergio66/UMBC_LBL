C     ********    O2 OXYGEN COLLISION INDUCED FUNDAMENTAL  ***********  
c
c     version_1 of the Oxygen Collision Induced Fundamental
c
c     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann,
c        and Ch. Boulet,
c        Infrared collision-induced absorption by O2 near 6.4 microns for
c        atmospheric applications: measurements and emprirical modeling, 
c        Appl. Optics, 35, 5911-5917, (1996).

c
C        Only calculate if V2 > 1340. cm-1 and V1 <  1850. cm-1

         if (((V2.gt.1340.0).and.(V1.lt.1850.))) then
c     
            tau_fac = xo2cn *  Wk(7) * 1.e-20 * amagat 
c
c           Wk(7) is the oxygen column amount in units of molec/cm2
c           amagat is in units of amagats (air)
c
c           The temperature correction is done in the subroutine o2_ver_1:
c
            call o2_ver_1 (v1c,v2c,dvc,nptc,c0,tave)
C
c           c0 are the oxygen absorption coefficients at temperature tave 
c              - these absorption coefficients are in units of
c                   [(cm^2/molec) 10^20)]/(cm-1  amagat) 
c              - cm-1 in the denominator arises through the removal
c                   of the radiation field
c              - for this case, an amagat is interpreted as one
c                   loshmidt of air (273K)
c
            DO 80 J = 1, NPTC
               VJ = V1C+DVC* REAL(J-1)
               C(J) = tau_fac * c0(J) 
C
C              Radiation field
C
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
c
 80         CONTINUE

            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)
         endif

C        ********    O2 Collision Induced   ********    
C
C        O2 continuum formulated by Mate et al. over the spectral region
C        7550-8486 cm-1:  "Absolute Intensities for the O2 1.27 micron
C        continuum absorption", B. Mate, C. Lugez, G.T. Fraser, and
C        W.J. Lafferty, J. Geophys. Res., 104, 30,585-30,590, 1999. 
c
c        The units of these continua coefficients are  1 / (amagat_O2*amag
c
c        Also, refer to the paper "Observed  Atmospheric
C        Collision Induced Absorption in Near Infrared Oxygen Bands",
C        Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
C        Journal of Geophysical Research (1997).
C
C        Only calculate if V2 > 7536. cm-1 and V1 <  8500. cm-1
c
         if (((V2.gt.7536.0).and.(V1.lt.8500.))) then
c
            a_o2  = 1./0.446
            a_n2  = 0.3/0.446
            a_h2o = 1.

            tau_fac = xn2cn * (Wk(7)/xlosmt) * amagat * 
     &           (a_o2*x_vmr_o2+a_n2*x_vmr_n2+a_h2o*x_vmr_h2o)

c
            CALL O2INF1 (V1C,V2C,DVC,NPTC,C0)                     
c
            DO 92 J = 1, NPTC                                             
               C(J) = tau_fac * C0(J)
               VJ = V1C+DVC* REAL(J-1)                                    
C                                                                      
C              Radiation field
C                                                                      
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)     
                 
 92         CONTINUE                                                      
c 
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      
c
         endif
C
C        O2 continuum formulated by Mlawer et al. over the spectral region
C        9100-11000 cm-1. Refer to the paper "Observed  Atmospheric
C        Collision Induced Absorption in Near Infrared Oxygen Bands",
C        Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
C        Journal of Geophysical Research (1997).
C
C        Only calculate if V2 > 9100. cm-1 and V1 <  11000. cm-1
c
         if ((V2.gt.9100.0).and.(V1.lt.11000.)) then
c
            CALL O2INF2 (V1C,V2C,DVC,NPTC,C0)                      
            WO2 = xo2cn * (WK(7)*1.e-20) * RHOAVE
            ADJWO2 = (WK(7)/WTOT) * (1./0.209) * WO2
c
            DO 93 J = 1, NPTC                                             
               C(J) = C0(J)*ADJWO2
               VJ = V1C+DVC* REAL(J-1)                                    
C                                                                      
C              Radiation field
C                                                                      
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   
c
 93         CONTINUE                                                      
c
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      
c
         endif
C
C        O2 continuum formulated by Greenblatt et al. over the spectral re
C        8797-29870 cm-1:  "Absorption Coefficients of Oxygen Between 
c        330 and 1140 nm, G.D. Green blatt, J.J. Orlando, J.B. Burkholder,
c        and A.R. Ravishabkara,  J. Geophys. Res., 95, 18577-18582, 1990. 
c
c        The units conversion to (cm^2/molec)/atm(o2)  has been done in 
c        subroutine o2_vis
C
C        Only calculate if V2 > 15000. cm-1 and V1 <  29870. cm-1
c
         if (((V2.gt.15000.0).and.(V1.lt.29870.))) then
c
            WO2 = WK(7) * 1.e-20 * ((pave/1013.)*(273./tave)) * xo2cn
            CHIO2 =  WK(7)/WTOT 
            ADJWO2 = chio2 * WO2
c
            CALL O2_vis (V1C,V2C,DVC,NPTC,C0)                     
c
            DO 94 J = 1, NPTC                                             
               C(J) = C0(J)*ADJWO2
               VJ = V1C+DVC* REAL(J-1)                                    
C                                                                      
C              Radiation field
C                                                                      
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   
c
 94         CONTINUE                                                      
c 
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      
c
         endif
c
C        Only calculate if V2 > 36000. cm-1

         if (V2.gt.36000.0) then

            WO2 = WK(7) * 1.e-20 * xo2cn
c
            CALL O2HERZ (V1C,V2C,DVC,NPTC,C0,TAVE,PAVE)                   
            DO 90 J = 1, NPTC                                             
               C(J) = C0(J)*WO2
               VJ = V1C+DVC* REAL(J-1)                                    
C                                                                         
C              Radiation field                                            
C                                                                         
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   
 90         CONTINUE                                                      
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      

         endif

      RETURN                                                              

C                                                                         
      END                                                                 
C
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      subroutine o2_ver_1 (v1c,v2c,dvc,nptc,c,T)
c
      IMPLICIT REAL*8 (v)                                                 

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)

      COMMON /o2_f  / V1S,V2S,DVS,NPTS,xo2(103),xo2t(103)

      dimension c(*)

c
c     Oxygen Collision Induced Fundamental

c     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann,
c                                                         and Ch. Boulet
c        Infrared collision-induced absorption by O2 near 6.4 microns for
c        atmospheric applications: measurements and emprirical modeling, 
c         Appl. Optics, 35, 5911-5917, (1996).

      DATA T_0/ 296./, xlosmt/ 2.68675e+19/
c
      xktfac = (1./T_0)-(1./T)
c     
c     correct formulation for consistency with LBLRTM:
c
      factor = (1.e+20 /xlosmt) 
c
c     A factor of 0.21, the mixing ration of oxygen, in the Thibault et al
c     formulation is not included here.  This factor is in the column amou
C                           
      DVC = DVS             
      V1C = V1ABS-DVC       
      V2C = V2ABS+DVC       
C                           
      IF (V1C.LT.V1S) then
         I1 = -1 
      else
         I1 = (V1C-V1S)/DVS + 0.01
      end if
C                                    
      V1C = V1S + DVS*REAL(I1-1)        
      I2 = (V2C-V1S)/DVS + 0.01            
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C + DVS*REAL(NPTC-1)       
c
      do 10 j=1,nptc
         i = i1+j
         C(J) = 0.                                                        
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            
         VJ = V1C+DVC* REAL(J-1)                                          
c     the radiation field is removed with 1/vj
c
         c(j) = factor * xo2(i)* exp(xo2t(i)*xktfac) / vj
c
 10   end do
c
 920  format (f10.2,1p,e12.2,0p,f10.2,1p2e12.2)
c
      return
      end

      BLOCK DATA bo2f                                                   
                                                                        
      IMPLICIT REAL*8 (V)                                               
                                                                        
      COMMON /o2_f  / V1S,V2S,DVS,NPTS,
     *          o0001(50),o0051(50),o0101(03),
     *          ot0001(50),ot0051(50),ot0101(03)
                                                                        
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

C     --------------------------------------------------------------
C
      SUBROUTINE O2INF1 (V1C,V2C,DVC,NPTC,C)                        
C                                                                        
      IMPLICIT REAL*8           (V)                                     
C                                                                        
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)               
      DIMENSION C(*)                                                     

      COMMON /o2inf1_mate/ V1S,V2S,DVS,NPTS,xo2inf1(483)

C                                                                        
C        O2 continuum formulated by Mate et al. over the spectral region
C        7550-8486 cm-1:  "Absolute Intensities for the O2 1.27 micron
C        continuum absorption", B. Mate, C. Lugez, G.T. Fraser, and
C        W.J. Lafferty, J. Geophys. Res., 104, 30,585-30,590, 1999. 
c
c        The units of these continua coefficients are  1 / (amagat_O2*amag
c
c        Also, refer to the paper "Observed  Atmospheric
C        Collision Induced Absorption in Near Infrared Oxygen Bands",
C        Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
C        Journal of Geophysical Research (1997).
c   ***********

      DVC = DVS                                                          
C                                                                        
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
C                                                                        
      IF (V1C.LT.V1S) then
         I1 = -1 
      else
         I1 = (V1C-V1S)/DVS + 0.01
      end if
C                                    
      V1C = V1S + DVS*REAL(I1-1)        
      I2 = (V2C-V1S)/DVS + 0.01            
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C + DVS*REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   
         I = I1+(J-1)                                                     
         C(J) = 0.                                                        
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            
         vj = v1c + dvc* REAL(j-1)
         C(J) = xo2inf1(I)/vj
   10 CONTINUE                                                            
C                                                                         
      RETURN                                                             
C                                                                        
      END                                                                

C     --------------------------------------------------------------
C
      BLOCK DATA bo2inf1
                                                                        
      IMPLICIT REAL*8 (V)                                               
                                                                        
      COMMON /o2inf1_mate/ V1,V2,DV,NPT,                                  
     *          o0001(50),o0051(50),o0101(50),o0151(50),o0201(50),      
     *          o0251(50),o0301(50),o0351(50),o0401(50),o0451(33)      
                                                                        
      DATA V1,V2,DV,NPT /7536.000,8500.000,   2.000,  483/              
                                                                        
      DATA o0001/                                                       
     *      0.000E+00,  4.355E-11,  8.709E-11,  1.742E-10,  3.484E-10,  
     *      6.968E-10,  1.394E-09,  2.787E-09,  3.561E-09,  3.314E-09,  
     *      3.368E-09,  3.435E-09,  2.855E-09,  3.244E-09,  3.447E-09,  
     *      3.891E-09,  4.355E-09,  3.709E-09,  4.265E-09,  4.772E-09,  
     *      4.541E-09,  4.557E-09,  4.915E-09,  4.688E-09,  5.282E-09,  
     *      5.755E-09,  5.096E-09,  5.027E-09,  4.860E-09,  4.724E-09,  
     *      5.048E-09,  5.248E-09,  5.473E-09,  4.852E-09,  5.362E-09,  
     *      6.157E-09,  6.150E-09,  6.347E-09,  6.388E-09,  6.213E-09,  
     *      6.521E-09,  8.470E-09,  8.236E-09,  8.269E-09,  8.776E-09,  
     *      9.122E-09,  9.189E-09,  9.778E-09,  8.433E-09,  9.964E-09/  
      DATA o0051/                                                       
     *      9.827E-09,  1.064E-08,  1.063E-08,  1.031E-08,  1.098E-08,  
     *      1.156E-08,  1.295E-08,  1.326E-08,  1.467E-08,  1.427E-08,  
     *      1.452E-08,  1.456E-08,  1.554E-08,  1.605E-08,  1.659E-08,  
     *      1.754E-08,  1.757E-08,  1.876E-08,  1.903E-08,  1.876E-08,  
     *      1.869E-08,  2.036E-08,  2.203E-08,  2.221E-08,  2.284E-08,  
     *      2.288E-08,  2.394E-08,  2.509E-08,  2.663E-08,  2.720E-08,  
     *      2.839E-08,  2.923E-08,  2.893E-08,  2.949E-08,  2.962E-08,  
     *      3.057E-08,  3.056E-08,  3.364E-08,  3.563E-08,  3.743E-08,  
     *      3.813E-08,  3.946E-08,  4.082E-08,  4.201E-08,  4.297E-08,  
     *      4.528E-08,  4.587E-08,  4.704E-08,  4.962E-08,  5.115E-08/  
      DATA o0101/                                                       
     *      5.341E-08,  5.365E-08,  5.557E-08,  5.891E-08,  6.084E-08,  
     *      6.270E-08,  6.448E-08,  6.622E-08,  6.939E-08,  7.233E-08,  
     *      7.498E-08,  7.749E-08,  8.027E-08,  8.387E-08,  8.605E-08,  
     *      8.888E-08,  9.277E-08,  9.523E-08,  9.880E-08,  1.037E-07,  
     *      1.076E-07,  1.114E-07,  1.151E-07,  1.203E-07,  1.246E-07,  
     *      1.285E-07,  1.345E-07,  1.408E-07,  1.465E-07,  1.519E-07,  
     *      1.578E-07,  1.628E-07,  1.685E-07,  1.760E-07,  1.847E-07,  
     *      1.929E-07,  2.002E-07,  2.070E-07,  2.177E-07,  2.262E-07,  
     *      2.365E-07,  2.482E-07,  2.587E-07,  2.655E-07,  2.789E-07,  
     *      2.925E-07,  3.023E-07,  3.153E-07,  3.296E-07,  3.409E-07/  
      DATA o0151/                                                       
     *      3.532E-07,  3.680E-07,  3.859E-07,  3.951E-07,  4.074E-07,  
     *      4.210E-07,  4.381E-07,  4.588E-07,  4.792E-07,  4.958E-07,  
     *      5.104E-07,  5.271E-07,  5.501E-07,  5.674E-07,  5.913E-07,  
     *      6.243E-07,  6.471E-07,  6.622E-07,  6.831E-07,  6.987E-07,  
     *      7.159E-07,  7.412E-07,  7.698E-07,  7.599E-07,  7.600E-07,  
     *      7.918E-07,  8.026E-07,  8.051E-07,  8.049E-07,  7.914E-07,  
     *      7.968E-07,  7.945E-07,  7.861E-07,  7.864E-07,  7.741E-07,  
     *      7.675E-07,  7.592E-07,  7.400E-07,  7.362E-07,  7.285E-07,  
     *      7.173E-07,  6.966E-07,  6.744E-07,  6.597E-07,  6.413E-07,  
     *      6.265E-07,  6.110E-07,  5.929E-07,  5.717E-07,  5.592E-07/  
      DATA o0201/                                                       
     *      5.411E-07,  5.235E-07,  5.061E-07,  4.845E-07,  4.732E-07,  
     *      4.593E-07,  4.467E-07,  4.328E-07,  4.161E-07,  4.035E-07,  
     *      3.922E-07,  3.820E-07,  3.707E-07,  3.585E-07,  3.475E-07,  
     *      3.407E-07,  3.317E-07,  3.226E-07,  3.134E-07,  3.016E-07,  
     *      2.969E-07,  2.894E-07,  2.814E-07,  2.749E-07,  2.657E-07,  
     *      2.610E-07,  2.536E-07,  2.467E-07,  2.394E-07,  2.337E-07,  
     *      2.302E-07,  2.241E-07,  2.191E-07,  2.140E-07,  2.093E-07,  
     *      2.052E-07,  1.998E-07,  1.963E-07,  1.920E-07,  1.862E-07,  
     *      1.834E-07,  1.795E-07,  1.745E-07,  1.723E-07,  1.686E-07,  
     *      1.658E-07,  1.629E-07,  1.595E-07,  1.558E-07,  1.523E-07/  
      DATA o0251/                                                       
     *      1.498E-07,  1.466E-07,  1.452E-07,  1.431E-07,  1.408E-07,  
     *      1.381E-07,  1.362E-07,  1.320E-07,  1.298E-07,  1.262E-07,  
     *      1.247E-07,  1.234E-07,  1.221E-07,  1.197E-07,  1.176E-07,  
     *      1.142E-07,  1.121E-07,  1.099E-07,  1.081E-07,  1.073E-07,  
     *      1.061E-07,  1.041E-07,  1.019E-07,  9.969E-08,  9.727E-08,  
     *      9.642E-08,  9.487E-08,  9.318E-08,  9.116E-08,  9.046E-08,  
     *      8.827E-08,  8.689E-08,  8.433E-08,  8.324E-08,  8.204E-08,  
     *      8.036E-08,  7.951E-08,  7.804E-08,  7.524E-08,  7.392E-08,  
     *      7.227E-08,  7.176E-08,  6.975E-08,  6.914E-08,  6.859E-08,  
     *      6.664E-08,  6.506E-08,  6.368E-08,  6.262E-08,  6.026E-08/  
      DATA o0301/                                                       
     *      6.002E-08,  5.866E-08,  5.867E-08,  5.641E-08,  5.589E-08,  
     *      5.499E-08,  5.309E-08,  5.188E-08,  5.139E-08,  4.991E-08,  
     *      4.951E-08,  4.833E-08,  4.640E-08,  4.524E-08,  4.479E-08,  
     *      4.304E-08,  4.228E-08,  4.251E-08,  4.130E-08,  3.984E-08,  
     *      3.894E-08,  3.815E-08,  3.732E-08,  3.664E-08,  3.512E-08,  
     *      3.463E-08,  3.503E-08,  3.218E-08,  3.253E-08,  3.107E-08,  
     *      2.964E-08,  2.920E-08,  2.888E-08,  2.981E-08,  2.830E-08,  
     *      2.750E-08,  2.580E-08,  2.528E-08,  2.444E-08,  2.378E-08,  
     *      2.413E-08,  2.234E-08,  2.316E-08,  2.199E-08,  2.088E-08,  
     *      1.998E-08,  1.920E-08,  1.942E-08,  1.859E-08,  1.954E-08/  
      DATA o0351/                                                       
     *      1.955E-08,  1.749E-08,  1.720E-08,  1.702E-08,  1.521E-08,  
     *      1.589E-08,  1.469E-08,  1.471E-08,  1.543E-08,  1.433E-08,  
     *      1.298E-08,  1.274E-08,  1.226E-08,  1.204E-08,  1.201E-08,  
     *      1.298E-08,  1.220E-08,  1.220E-08,  1.096E-08,  1.080E-08,  
     *      9.868E-09,  9.701E-09,  1.130E-08,  9.874E-09,  9.754E-09,  
     *      9.651E-09,  9.725E-09,  8.413E-09,  7.705E-09,  7.846E-09,  
     *      8.037E-09,  9.163E-09,  8.098E-09,  8.160E-09,  7.511E-09,  
     *      7.011E-09,  6.281E-09,  6.502E-09,  7.323E-09,  7.569E-09,  
     *      5.941E-09,  5.867E-09,  5.676E-09,  4.840E-09,  5.063E-09,  
     *      5.207E-09,  4.917E-09,  5.033E-09,  5.356E-09,  3.795E-09/  
      DATA o0401/                                                       
     *      4.983E-09,  4.600E-09,  3.635E-09,  3.099E-09,  2.502E-09,  
     *      3.823E-09,  3.464E-09,  4.332E-09,  3.612E-09,  3.682E-09,  
     *      3.709E-09,  3.043E-09,  3.593E-09,  3.995E-09,  4.460E-09,  
     *      3.583E-09,  3.290E-09,  3.132E-09,  2.812E-09,  3.109E-09,  
     *      3.874E-09,  3.802E-09,  4.024E-09,  3.901E-09,  2.370E-09,  
     *      1.821E-09,  2.519E-09,  4.701E-09,  3.855E-09,  4.685E-09,  
     *      5.170E-09,  4.387E-09,  4.148E-09,  4.043E-09,  3.545E-09,  
     *      3.392E-09,  3.609E-09,  4.635E-09,  3.467E-09,  2.558E-09,  
     *      3.389E-09,  2.672E-09,  2.468E-09,  1.989E-09,  2.816E-09,  
     *      4.023E-09,  2.664E-09,  2.219E-09,  3.169E-09,  1.654E-09/  
      DATA o0451/                                                       
     *      3.189E-09,  2.535E-09,  2.618E-09,  3.265E-09,  2.138E-09,  
     *      1.822E-09,  2.920E-09,  2.002E-09,  1.300E-09,  3.764E-09,  
     *      3.212E-09,  3.222E-09,  2.961E-09,  2.108E-09,  1.708E-09,  
     *      2.636E-09,  2.937E-09,  2.939E-09,  2.732E-09,  2.218E-09,  
     *      1.046E-09,  6.419E-10,  1.842E-09,  1.112E-09,  1.265E-09,  
     *      4.087E-09,  2.044E-09,  1.022E-09,  5.109E-10,  2.554E-10,  
     *      1.277E-10,  6.386E-11,  0.000E+00/                          

      END


C     --------------------------------------------------------------
C
      SUBROUTINE O2INF2 (V1C,V2C,DVC,NPTC,C)                         
C                                                                        
      IMPLICIT REAL*8           (V)                                     
C                                                                        
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)               
      DIMENSION C(*)                                                     
C                                                                        
      DATA V1_osc /9375./, HW1 /58.96/, V2_osc /9439./, HW2 /45.04/
      DATA S1 /1.166E-04/, S2 /3.086E-05/
C                                                                        
      V1S = 9100.                                                        
      v2s = 11000.
      DVS = 2.                                                          
      DVC = DVS                                                          
C                                                                        
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
C                                                                        
      NPTC = (v2c-v1c)/dvc + 3.01
      V2C = V1C+DVc* REAL(NPTC-1)                                        
c
      DO 10 J = 1, NPTC                                                  
         C(J) = 0.                                                       
         VJ = V1C+DVC* REAL(J-1)                                         
         IF ((Vj.gt.v1s) .and. (Vj.lt.v2s)) then
            DV1 = Vj - V1_osc
            DV2 = Vj - V2_osc
            IF (DV1 .LT. 0.0) THEN
               DAMP1 = EXP (DV1 / 176.1)
            ELSE
               DAMP1 = 1.0
            ENDIF
            IF (DV2 .LT. 0.0) THEN
               DAMP2 = EXP (DV2 / 176.1)
            ELSE
               DAMP2 = 1.0
            ENDIF
            O2INF = 0.31831 * (((S1 * DAMP1 / HW1)/(1. + (DV1/HW1)**2))
     *           + ((S2 * DAMP2 / HW2)/(1. + (DV2/HW2)**2))) * 1.054
            C(J) = O2INF/VJ  
         endif
   10 CONTINUE                                                           
C                                                                        
      RETURN                                                             
C                                                                        
      END                                                                

C     --------------------------------------------------------------
C
      SUBROUTINE O2_vis (V1C,V2C,DVC,NPTC,C)                        
C                                                                        
      IMPLICIT REAL*8           (V)                                     
C                                                                        
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)               
      DIMENSION C(*)                                                     
      COMMON /o2_o2_vis/ V1s,V2s,DVs,NPTs, s(1488)                        

      DATA XLOSMT / 2.68675E+19 /   
C
C        O2 continuum formulated by Greenblatt et al. over the spectral re
C        8797-29870 cm-1:  "Absorption Coefficients of Oxygen Between 
c        330 and 1140 nm, G.D. Green blatt, J.J. Orlando, J.B. Burkholder,
c        and A.R. Ravishabkara,  J. Geophys. Res., 95, 18577-18582, 1990. 
c
c        The units conversion  is to (cm^2/molec)/atm(o2)
c
c      these are the conditions reported in the paper by Greenblatt et al.
c     the spectrum of Fig. 1.
c
c     conditions:  55 atm.; 296 K; 89.5 cm path
c
      factor = 1./((xlosmt*1.e-20*(55.*273./296.)**2)*89.5)
c
      DVC = DVS                                                          
C                                                                        
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
C                                                                        
      IF (V1C.LT.V1S) then
         I1 = -1 
      else
         I1 = (V1C-V1S)/DVS + 0.01
      end if
C                                    
      V1C = V1S + DVS*REAL(I1-1)        
      I2 = (V2C-V1S)/DVS + 0.01            
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C + DVS*REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   
         I = I1+(J-1)                                                     
         C(J) = 0.                                                        
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            
         vj = v1c + dvc* REAL(j-1)
c
         C(J) = factor*S(I)/vj  

   10 CONTINUE                                                            
C                                                                         
      RETURN                                                             
C                                                                        
      END                                                                
C
C     --------------------------------------------------------------
C
      BLOCK DATA bo2in_vis
                                                                        
      IMPLICIT REAL*8 (V)                                               
                                                                        
      COMMON /o2_o2_vis/ V1,V2,DV,NPT,                                    
     *  o2vis0001(50),o2vis0051(50),o2vis0101(50),o2vis0151(50),
     *  o2vis0201(50),o2vis0251(50),o2vis0301(50),o2vis0351(50),
     *  o2vis0401(50),o2vis0451(50),o2vis0501(50),o2vis0551(50),
     *  o2vis0601(50),o2vis0651(50),o2vis0701(50),o2vis0751(50),
     *  o2vis0801(50),o2vis0851(50),o2vis0901(50),o2vis0951(50),
     *  o2vis1001(50),o2vis1051(50),o2vis1101(50),o2vis1151(50),
     *  o2vis1201(50),o2vis1251(50),o2vis1301(50),o2vis1351(50),
     *  o2vis1401(50),o2vis1451(38)
                                                                        
      DATA V1,V2,DV,NPT /15000.0, 29870.0, 10.0,  1488/              
                                                                        
      DATA o2vis0001/
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   0.00E+00,
     *      0.00E+00,   0.00E+00,   6.06E-04,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.49E-03,   3.00E-03,   3.00E-03,   3.00E-03,   4.00E-03,
     *      4.00E-03,   5.00E-03,   5.00E-03,   6.00E-03,   7.00E-03/
      DATA o2vis0051/
     *      8.00E-03,   9.00E-03,   1.00E-02,   1.10E-02,   1.25E-02,
     *      1.46E-02,   1.60E-02,   1.80E-02,   2.00E-02,   2.23E-02,
     *      2.50E-02,   2.69E-02,   3.00E-02,   3.30E-02,   3.63E-02,
     *      4.01E-02,   4.42E-02,   4.67E-02,   5.14E-02,   5.55E-02,
     *      5.96E-02,   6.43E-02,   6.94E-02,   7.37E-02,   7.88E-02,
     *      8.38E-02,   8.86E-02,   9.37E-02,   9.89E-02,   1.03E-01,
     *      1.07E-01,   1.10E-01,   1.14E-01,   1.16E-01,   1.18E-01,
     *      1.19E-01,   1.20E-01,   1.21E-01,   1.20E-01,   1.20E-01,
     *      1.19E-01,   1.17E-01,   1.16E-01,   1.13E-01,   1.10E-01,
     *      1.07E-01,   1.03E-01,   9.97E-02,   9.58E-02,   9.15E-02/
      DATA o2vis0101/
     *      8.80E-02,   8.41E-02,   7.94E-02,   7.53E-02,   7.17E-02,
     *      6.83E-02,   6.43E-02,   6.08E-02,   5.69E-02,   5.31E-02,
     *      5.02E-02,   4.77E-02,   4.40E-02,   4.23E-02,   3.94E-02,
     *      3.70E-02,   3.51E-02,   3.30E-02,   3.10E-02,   2.90E-02,
     *      2.79E-02,   2.60E-02,   2.50E-02,   2.32E-02,   2.20E-02,
     *      2.10E-02,   2.00E-02,   1.90E-02,   1.80E-02,   1.70E-02,
     *      1.65E-02,   1.50E-02,   1.40E-02,   1.30E-02,   1.30E-02,
     *      1.20E-02,   1.10E-02,   1.10E-02,   1.00E-02,   1.00E-02,
     *      9.00E-03,   9.00E-03,   9.00E-03,   8.00E-03,   8.00E-03,
     *      7.01E-03,   7.00E-03,   7.00E-03,   6.98E-03,   6.00E-03/
      DATA o2vis0151/
     *      5.80E-03,   5.00E-03,   5.00E-03,   5.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   5.00E-03,
     *      5.00E-03,   6.00E-03,   6.00E-03,   7.00E-03,   7.41E-03,
     *      8.15E-03,   9.00E-03,   1.01E-02,   1.10E-02,   1.20E-02/
      DATA o2vis0201/
     *      1.40E-02,   1.50E-02,   1.70E-02,   1.85E-02,   1.97E-02,
     *      2.24E-02,   2.47E-02,   2.74E-02,   3.06E-02,   3.36E-02,
     *      3.70E-02,   4.05E-02,   4.49E-02,   4.93E-02,   5.47E-02,
     *      6.01E-02,   6.52E-02,   7.23E-02,   7.89E-02,   8.80E-02,
     *      9.61E-02,   1.05E-01,   1.17E-01,   1.26E-01,   1.39E-01,
     *      1.49E-01,   1.60E-01,   1.68E-01,   1.74E-01,   1.79E-01,
     *      1.82E-01,   1.84E-01,   1.85E-01,   1.84E-01,   1.83E-01,
     *      1.81E-01,   1.80E-01,   1.77E-01,   1.74E-01,   1.71E-01,
     *      1.68E-01,   1.64E-01,   1.60E-01,   1.55E-01,   1.51E-01,
     *      1.46E-01,   1.40E-01,   1.36E-01,   1.30E-01,   1.25E-01/
      DATA o2vis0251/
     *      1.20E-01,   1.14E-01,   1.09E-01,   1.05E-01,   9.93E-02,
     *      9.30E-02,   8.88E-02,   8.38E-02,   7.94E-02,   7.51E-02,
     *      7.08E-02,   6.66E-02,   6.32E-02,   6.01E-02,   5.55E-02,
     *      5.24E-02,   4.93E-02,   4.63E-02,   4.41E-02,   4.15E-02,
     *      3.90E-02,   3.63E-02,   3.50E-02,   3.26E-02,   3.05E-02,
     *      2.94E-02,   2.73E-02,   2.62E-02,   2.46E-02,   2.36E-02,
     *      2.25E-02,   2.10E-02,   2.00E-02,   1.90E-02,   1.80E-02,
     *      1.76E-02,   1.70E-02,   1.60E-02,   1.50E-02,   1.49E-02,
     *      1.40E-02,   1.30E-02,   1.30E-02,   1.22E-02,   1.20E-02,
     *      1.20E-02,   1.10E-02,   1.10E-02,   1.10E-02,   1.00E-02/
      DATA o2vis0301/
     *      1.00E-02,   1.00E-02,   1.00E-02,   9.16E-03,   9.00E-03,
     *      9.00E-03,   9.00E-03,   9.00E-03,   8.49E-03,   8.00E-03,
     *      8.00E-03,   8.00E-03,   8.00E-03,   8.00E-03,   8.00E-03,
     *      8.00E-03,   7.00E-03,   8.00E-03,   7.00E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   8.00E-03,   8.00E-03,   8.00E-03/
      DATA o2vis0351/
     *      8.00E-03,   8.00E-03,   8.00E-03,   9.00E-03,   9.00E-03,
     *      9.00E-03,   9.07E-03,   1.00E-02,   1.00E-02,   1.00E-02,
     *      1.10E-02,   1.10E-02,   1.20E-02,   1.22E-02,   1.30E-02,
     *      1.31E-02,   1.40E-02,   1.50E-02,   1.60E-02,   1.70E-02,
     *      1.82E-02,   2.00E-02,   2.01E-02,   2.10E-02,   2.20E-02,
     *      2.28E-02,   2.30E-02,   2.30E-02,   2.30E-02,   2.30E-02,
     *      2.30E-02,   2.30E-02,   2.30E-02,   2.20E-02,   2.20E-02,
     *      2.20E-02,   2.10E-02,   2.10E-02,   2.00E-02,   2.00E-02,
     *      1.90E-02,   1.90E-02,   1.82E-02,   1.80E-02,   1.74E-02,
     *      1.70E-02,   1.63E-02,   1.60E-02,   1.50E-02,   1.49E-02/
      DATA o2vis0401/
     *      1.40E-02,   1.37E-02,   1.30E-02,   1.30E-02,   1.21E-02,
     *      1.20E-02,   1.13E-02,   1.09E-02,   1.00E-02,   9.34E-03,
     *      9.00E-03,   8.43E-03,   8.00E-03,   7.39E-03,   7.00E-03,
     *      6.00E-03,   6.00E-03,   5.74E-03,   5.00E-03,   5.00E-03,
     *      5.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      3.17E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03/
      DATA o2vis0451/
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      1.04E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
      DATA o2vis0501/
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.41E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   1.98E-03,   1.46E-03,
     *      1.05E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
      DATA o2vis0551/
     *      1.00E-03,   1.00E-03,   1.71E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   3.00E-03,   3.00E-03,
     *      3.82E-03,   4.00E-03,   4.17E-03,   5.00E-03,   6.00E-03,
     *      7.00E-03,   7.73E-03,   8.07E-03,   9.70E-03,   1.17E-02,
     *      1.31E-02,   1.47E-02,   1.64E-02,   1.81E-02,   2.07E-02,
     *      2.37E-02,   2.70E-02,   2.97E-02,   3.27E-02,   3.70E-02,
     *      4.13E-02,   4.49E-02,   4.89E-02,   5.38E-02,   5.98E-02,
     *      6.45E-02,   6.94E-02,   7.41E-02,   8.01E-02,   8.51E-02,
     *      9.00E-02,   9.49E-02,   9.88E-02,   1.01E-01,   1.04E-01,
     *      1.07E-01,   1.07E-01,   1.06E-01,   1.03E-01,   1.00E-01/
      DATA o2vis0601/
     *      9.66E-02,   8.93E-02,   8.35E-02,   7.92E-02,   7.33E-02,
     *      6.84E-02,   6.40E-02,   5.91E-02,   5.57E-02,   5.26E-02,
     *      5.03E-02,   4.75E-02,   4.48E-02,   4.26E-02,   4.07E-02,
     *      3.83E-02,   3.69E-02,   3.47E-02,   3.24E-02,   3.11E-02,
     *      2.85E-02,   2.69E-02,   2.55E-02,   2.42E-02,   2.21E-02,
     *      2.09E-02,   1.93E-02,   1.77E-02,   1.62E-02,   1.60E-02,
     *      1.44E-02,   1.36E-02,   1.30E-02,   1.16E-02,   1.10E-02,
     *      1.00E-02,   1.00E-02,   9.00E-03,   8.27E-03,   8.00E-03,
     *      7.45E-03,   7.00E-03,   7.00E-03,   6.18E-03,   6.00E-03,
     *      6.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03/
      DATA o2vis0651/
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   2.07E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   1.28E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03/
      DATA o2vis0701/
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.57E-03,   5.00E-03,
     *      5.00E-03,   5.64E-03,   6.00E-03,   6.67E-03,   7.00E-03,
     *      7.35E-03,   8.00E-03,   8.36E-03,   9.00E-03,   9.00E-03,
     *      1.00E-02,   1.00E-02,   1.00E-02,   1.00E-02,   1.00E-02,
     *      1.00E-02,   1.00E-02,   9.65E-03,   9.00E-03,   9.00E-03,
     *      8.00E-03,   8.00E-03,   7.69E-03,   7.00E-03,   7.00E-03/
      DATA o2vis0751/
     *      6.44E-03,   6.00E-03,   6.00E-03,   6.00E-03,   5.00E-03,
     *      5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   3.98E-03,
     *      3.01E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   2.54E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03/
      DATA o2vis0801/
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      1.33E-03,   1.89E-03,   1.07E-03,   1.06E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
      DATA o2vis0851/
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
      DATA o2vis0901/
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   5.50E-04,
     *      0.00E+00,   0.00E+00,   1.00E-03,   1.00E-03,   7.51E-04,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00/
      DATA o2vis0951/
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   1.34E-04,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   7.65E-05,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
      DATA o2vis1001/
     *      1.00E-03,   1.20E-04,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00/
      DATA o2vis1051/
     *      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,
     *      0.00E+00,   6.09E-04,   3.47E-04,   6.97E-04,   2.60E-04,
     *      7.81E-04,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.68E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.76E-03,   3.00E-03,   3.00E-03,   3.00E-03/
      DATA o2vis1101/
     *      3.80E-03,   4.00E-03,   4.82E-03,   5.00E-03,   5.84E-03,
     *      6.00E-03,   6.85E-03,   7.85E-03,   8.86E-03,   9.86E-03,
     *      1.09E-02,   1.19E-02,   1.29E-02,   1.47E-02,   1.59E-02,
     *      1.77E-02,   1.97E-02,   2.09E-02,   2.27E-02,   2.47E-02,
     *      2.67E-02,   2.87E-02,   3.07E-02,   3.26E-02,   3.38E-02,
     *      3.56E-02,   3.68E-02,   3.86E-02,   3.90E-02,   3.98E-02,
     *      4.07E-02,   4.10E-02,   4.10E-02,   4.03E-02,   3.93E-02,
     *      3.83E-02,   3.73E-02,   3.64E-02,   3.48E-02,   3.34E-02,
     *      3.18E-02,   2.99E-02,   2.85E-02,   2.70E-02,   2.50E-02,
     *      2.31E-02,   2.11E-02,   1.92E-02,   1.76E-02,   1.63E-02/
      DATA o2vis1151/
     *      1.47E-02,   1.34E-02,   1.17E-02,   1.07E-02,   9.78E-03,
     *      8.81E-03,   7.84E-03,   6.88E-03,   6.00E-03,   5.94E-03,
     *      5.00E-03,   5.00E-03,   4.05E-03,   4.00E-03,   3.13E-03,
     *      3.00E-03,   3.00E-03,   2.24E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   1.54E-03,
     *      1.41E-03,   1.64E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
      DATA o2vis1201/
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,
     *      1.00E-03,   1.15E-03,   2.00E-03,   2.00E-03,   2.00E-03,
     *      2.00E-03,   2.00E-03,   2.00E-03,   2.56E-03,   3.00E-03,
     *      3.00E-03,   3.30E-03,   4.00E-03,   4.00E-03,   4.04E-03,
     *      4.95E-03,   5.85E-03,   6.00E-03,   6.67E-03,   7.58E-03,
     *      8.48E-03,   9.39E-03,   1.03E-02,   1.14E-02,   1.31E-02/
      DATA o2vis1251/
     *      1.40E-02,   1.58E-02,   1.76E-02,   1.94E-02,   2.12E-02,
     *      2.30E-02,   2.56E-02,   2.89E-02,   3.16E-02,   3.44E-02,
     *      3.80E-02,   4.16E-02,   4.52E-02,   4.87E-02,   5.23E-02,
     *      5.59E-02,   5.91E-02,   6.20E-02,   6.53E-02,   6.71E-02,
     *      6.89E-02,   6.98E-02,   7.07E-02,   7.10E-02,   7.10E-02,
     *      7.06E-02,   6.97E-02,   6.89E-02,   6.80E-02,   6.71E-02,
     *      6.54E-02,   6.43E-02,   6.29E-02,   6.11E-02,   5.94E-02,
     *      5.74E-02,   5.48E-02,   5.31E-02,   5.05E-02,   4.86E-02,
     *      4.62E-02,   4.41E-02,   4.23E-02,   4.03E-02,   3.78E-02,
     *      3.61E-02,   3.43E-02,   3.26E-02,   3.08E-02,   2.91E-02/
      DATA o2vis1301/
     *      2.73E-02,   2.58E-02,   2.49E-02,   2.31E-02,   2.22E-02,
     *      2.07E-02,   1.95E-02,   1.86E-02,   1.77E-02,   1.69E-02,
     *      1.60E-02,   1.51E-02,   1.43E-02,   1.40E-02,   1.35E-02,
     *      1.27E-02,   1.18E-02,   1.10E-02,   1.10E-02,   1.02E-02,
     *      1.00E-02,   1.00E-02,   9.67E-03,   8.81E-03,   8.05E-03,
     *      8.90E-03,   8.24E-03,   8.00E-03,   7.53E-03,   7.00E-03,
     *      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   6.42E-03,
     *      6.00E-03,   6.00E-03,   6.00E-03,   6.00E-03,   5.18E-03,
     *      5.00E-03,   5.00E-03,   5.00E-03,   4.80E-03,   4.04E-03,
     *      4.89E-03,   4.27E-03,   4.00E-03,   4.00E-03,   4.00E-03/
      DATA o2vis1351/
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   3.20E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,
     *      3.00E-03,   3.75E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.69E-03,   5.00E-03,   5.00E-03,
     *      5.15E-03,   5.97E-03,   6.00E-03,   6.61E-03,   7.43E-03,
     *      8.00E-03,   8.06E-03,   8.88E-03,   9.70E-03,   1.05E-02,
     *      1.13E-02,   1.21E-02,   1.30E-02,   1.38E-02,   1.52E-02/
      DATA o2vis1401/
     *      1.64E-02,   1.72E-02,   1.80E-02,   1.88E-02,   1.96E-02,
     *      2.04E-02,   2.10E-02,   2.10E-02,   2.10E-02,   2.10E-02,
     *      2.10E-02,   2.10E-02,   2.10E-02,   2.10E-02,   2.10E-02,
     *      2.05E-02,   2.00E-02,   1.99E-02,   1.91E-02,   1.90E-02,
     *      1.85E-02,   1.80E-02,   1.79E-02,   1.71E-02,   1.63E-02,
     *      1.55E-02,   1.47E-02,   1.40E-02,   1.40E-02,   1.33E-02,
     *      1.25E-02,   1.20E-02,   1.19E-02,   1.11E-02,   1.03E-02,
     *      1.00E-02,   9.75E-03,   9.00E-03,   9.00E-03,   8.37E-03,
     *      8.00E-03,   8.00E-03,   8.00E-03,   7.22E-03,   7.00E-03,
     *      7.00E-03,   6.86E-03,   6.07E-03,   6.00E-03,   6.00E-03/
      DATA o2vis1451/
     *      6.00E-03,   5.93E-03,   5.15E-03,   5.00E-03,   5.00E-03,
     *      5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,
     *      5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,
     *      5.00E-03,   5.00E-03,   5.00E-03,   4.68E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,
     *      1.00E-03,   2.00E-04,   0./
c
      end

C     --------------------------------------------------------------
C
      SUBROUTINE O2HERZ (V1C,V2C,DVC,NPTC,C,T,P)                          
C                                                                         
      IMPLICIT REAL*8           (V)                                     ! 
C                                                                         
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                
      DIMENSION C(*)                                                      
C                                                                         
      V1S = 36000.                                                        
      DVS = 10.                                                           
      DVC = DVS                                                           
C                                                                         
      V1C = V1ABS-DVC                                                     
      V2C = V2ABS+DVC                                                     
C                                                                         
      IF (V1C.LT.V1S) then
         I1 = -1 
      else
         I1 = (V1C-V1S)/DVS + 0.01
      end if
C                                    
      V1C = V1S + DVS*REAL(I1-1)        
      I2 = (V2C-V1S)/DVS + 0.01            
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C + DVS*REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   
         I = I1+(J-1)                                                     
         C(J) = 0.                                                        
         IF (I.LT.1) GO TO 10                                             
         VJ = V1C+DVC* REAL(J-1)                                          
         CALL HERTDA (HERZ,VJ)                                            
         CALL HERPRS (HERZ,T,P)                                           
         C(J) = HERZ/VJ                                                   
   10 CONTINUE                                                            
C                                                                         
      RETURN                                                              
C                                                                         
      END                                                                 
C
C     --------------------------------------------------------------
C
      SUBROUTINE HERTDA (HERZ,V)                                          
C                                                                         
      IMPLICIT REAL*8           (V)                                     ! 
C                                                                         
C     HERZBERG O2 ABSORPTION                                              
C     HALL,1987 PRIVATE COMMUNICATION, BASED ON:                          
C                                                                         
C     REF. JOHNSTON, ET AL., JGR,89,11661-11665,1984                      
C          NICOLET, 1987 (RECENT STUDIES IN ATOMIC                        
C                         & MOLECULAR PROCESSES,                          
C                         PLENUM PUBLISHING CORP, NY 1987)                
C                                                                         
C     AND YOSHINO, ET AL., 1988 (PREPRINT OF "IMPROVED ABSORPTION         
C          CROSS SECTIONS OF OXYGEN IN THE WAVELENGTH REGION 205-240NM    
C          OF THE HERZBERG CONTINUUM")                                    
C                                                                         
C     **** NOTE:  CROSS SECTION AT 0 PRESSURE  ***                        
C     THE PRESSURE DEPENDENT TERM IS IN SUBROUTINE HERPRS                 
C                                                                         
CC    COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                       
C                                                                         
      HERZ = 0.0                                                          
      IF (V.LE.36000.00) RETURN                                           
C                                                                         
C     EXTRAPOLATE SMOOTHLY THROUGH THE HERZBERG BAND REGION               
C     NOTE: HERZBERG BANDS ARE NOT CORRECTLY INCLUDED                     
C                                                                         
      CORR = 0.                                                           
      IF (V.LE.40000.) CORR = ((40000.-V)/4000.)*7.917E-27                
C                                                                         
C     UNITS ARE (CM2)                                                     
C                                                                         
C     HALL'S NEW HERZBERG  (LEAST SQRS FIT, LN(P))                        
C                                                                         
C     YRATIO=2048.7/WL(I)  ****IN ANGSTOMS****                            
C           =.20487/WN(I)     IN MICRONS                                  
C           =WCM(I)/48811.0   IN CM-1                                     
C                                                                         
      YRATIO = V/48811.0                                                  
c     HERZ = 6.884E-24*(YRATIO)*EXP(-69.738*( LOG(YRATIO))**2)-CORR       
c     factor of 1.e-20 removed; put in front factor
      HERZ = 6.884E-04*(YRATIO)*EXP(-69.738*( LOG(YRATIO))**2)-CORR       
C                                                                         
      RETURN                                                              
C                                                                         
      END                                                                 
C
C     --------------------------------------------------------------
C
      SUBROUTINE HERPRS (HERZ,T,P)                                        
C                                                                         
C     CORRECT THE HERZBERG CONTINUUM CROSS SECTION FOR PRESSURE           
C     DEPENDENCE; BASED ON SHARDANAND, JQRST, 18, 525-530, 1977.          
C                 FOR UN2| BROADENING                                     
C                 AND YOSHINO ET AL 1988 FOR UO2| BROADENING              
C                                                                         
C     PO2= PARTIAL PRESSURE OF O2                                         
C     PN2= PARTIAL PRESSURE OF N2; BN2=.45*BO2                            
C                                                                         
C     DATA BO2 / 1.72E-3 /                                                
C
C     Changed in Herzberg continuum pressure, 
C     Reference:
C     "Atmospheric Propagation in the UV, Visible, IR and MM-wave
C     Region and Related Systems Aspects".
C     G.P. Anderson,F.X. Kneizys, E.P. Shettle, L.W. Abreu,
C     J.H. Chetwynd, R.E. Huffman, and L.A. Hall; Conference
C     Proceedings No. 454 of the Advisory Group for Aerospace
C     Research & Development; 1990.
C
      DATA BO2 / 1.81E-3 /
      DATA PO / 1013. /,TO / 273.16 /                                     
C                                                                         
C     NOTE:  THE HERZBERG CONTINUUM OBEYS BEER'S LAW                      
C            OPTICAL DEPTH(TOTAL)=SUM OVER LAYER O.D.(I)                  
C                                                                         
C     BO2= RATIO OF SIGMA(O2-O2)/(SIGMA(O2)) * 760(TORR)*.2095            
C     BN2=.45*BO2= RATIO OF SIGMA(O2-N2)/(SIGMA(O2)) * 760(TORR)*.78      
C                                                                         
C     BO2*760*(.2095+.45*.78) = .73 , AS BELOW                            
C
C     Changed Herzberg continuum pressure (see above reference)
C
C     BO2*760*(.2095+.45*.78) = .83 , AS BELOW
C                                                                         
C                                                                         
      HERZ = HERZ*(1.+.83*(P/PO)*(TO/T))                                  
C                                                                         
      RETURN                                                              
C                                                                         
      END                                                                 


