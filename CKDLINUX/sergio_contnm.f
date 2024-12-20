c copied from LBLRTMv6.12 and modified so that it runs from run7.m
c sergio de souza-machado, Dec 2002

c main mods are to input 
      SUBROUTINE CONTNM(JRAD)                                             F00010
C                                                                         F00020
      IMPLICIT REAL*8           (V)                                     ! F00030
C                                                                         F00040
C     SUBROUTINE CONTNM CONTAINS THE CONTINUUM DATA                       F00050
C     WHICH IS INTERPOLATED INTO THE ARRAY ABSRB                          F00060
C                                                                         F00070
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F00080
      COMMON /XCONT/ V1C,V2C,DVC,NPTC,C(6000)
C                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
C                                                                         F00120
      COMMON /CVRCNT/ HVRCNT
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       F00130
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   F00140
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    F00150
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         F00170
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        F00180
     *              NLTEFL,LNFIL4,LNGTH4                                  F00190

      common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
c
      DIMENSION C0(5050),C1(5050),C2(5050)
      DIMENSION SH2OT0(5050),SH2OT1(5050),FH2O(5050),      
     *          CN2T0(5050),FCO2(5050),CT1(5050),CT2(5050) 
      DIMENSION CCH0(5150),CCH1(5150),CCH2(5150)
C
      REAL ABSBSV(5050)
C                                                                         F00230
      CHARACTER*15 HVRCNT
c
      equivalence (fscdid(4), iaersl)
c
      EQUIVALENCE (C0,SH2OT0,CN2T0,FCO2) , (C1,SH2OT1,CT1),               F00240
     *            (C2,FH2O,CT2)                                           F00250
C                                                                         F00260
      DATA P0 / 1013. /,T0 / 296. /                                       F00270
      DATA XLOSMT / 2.68675E+19 /                                         F00280
c     
c     These are self-continuum modification factors from 700-1200 cm-1
c
      DIMENSION XFAC(0:50)
c
      DATA (XFAC(I),I=0,50)/
     1    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
     2    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
     3    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
     4    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
     5    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
     6    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
     7    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
     8    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
     9    1.00000,1.00000,1.00000/
C                                                                         F00290
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
C     Continuum calculation flags:
C     ---------------------------
C     ICNTNM Value      Self     Foreign    Rayleigh     Others
C           0            no        no          no          no
C           1            yes       yes         yes         yes
C           2            no        yes         yes         yes
C           3            yes       no          yes         yes
C           4            no        no          yes         yes
C           5            yes       yes         no          yes
C           6   READ IN XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, 
C               and XRAYL in Record 1.2a
C
C     ASSIGN CVS  VERSION NUMBER TO MODULE 
c
      HVRCNT = '$Revision: 5.19 $'
C
      RHOAVE = (PAVE/P0)*(T0/TAVE)                                        F00300
      XKT = TAVE/RADCN2                                                   F00310
      WTOT = WBROAD                                                       F00320
      DO 10 M = 1, NMOL                                                   F00330
         WTOT = WTOT+WK(M)                                                F00340
   10 CONTINUE                                                            F00350
      WTOT = WTOT*1.E-20                                                  F00360
      PATM = PAVE/P0                                                      F00370
C                                                                         F00380
C=======================================================================

C               ********    WATER VAPOR   ********                        F00390
C                                                                         F00400
C=======================================================================
C                             SELF



C        Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
c
         if ((V2.gt.-20.0).and.(V1.lt.20000.)) then
c
            CALL SL296 (V1C,V2C,DVC,NPTC,SH2OT0)                          F00410
            CALL SL260 (V1C,V2C,DVC,NPTC,SH2OT1)                          F00420
C                                                                         F00440
            W1 = WK(1)*1.0E-20                                            F00450
C                                                                         F00460
C           The factor of 1.e-20 is handled this way to avoid underflows  F00470
C                                                                         F00480
            PH2O = PATM*(W1/WTOT)                                         F00490
            RH2O = PH2O*(T0/TAVE)                                         F00500
            XKT = TAVE/RADCN2                                             F00530
            TFAC = (TAVE-T0)/(260.-T0)                                    F00540

C
C--------------------------------------------------------------------
C     *****    Continuum Correction Patches    ********
C--------------------------------------------------------------------
C                             SELF
            V0S1 = 0.
            HWSQ1 = 100.**2
            BETAS1 = 1.E-04
c           FACTRS1 = 0.3                      ! CKD2.2 value
            FACTRS1 = 0.688
C
            V0S2 = 1050.
            HWSQ2 = 200.**2
            FACTRS2 = -0.2333
C
            V0S3 = 1310.
            HWSQ3 = 120.**2
            BETAS3 = 5.E-06
            FACTRS3 = -0.15
C--------------------------------------------------------------------
C
C
c           Loop calculating self continuum optical depth

            DO 20 J = 1, NPTC                                             F00560
               VJ = V1C+DVC* REAL(J-1)                                    F00570
               SH2O = 0.                                                  F00580
               IF (SH2OT0(J).GT.0.) THEN                                  F00590
                  SH2O = SH2OT0(J)*(SH2OT1(J)/SH2OT0(J))**TFAC            F00600
C     
                  SFAC = 1.
                  IF (VJ.GE.700. .AND.  VJ.LE.1200.) THEN 
                     JFAC = (VJ-700.)/10. + 0.00001
                     SFAC = XFAC(JFAC)
                  ENDIF
C     
C              ---------------------------------------------------------
C              Correction to self continuum (1 SEPT 85); factor of    
C                 0.78 at 1000 and  .......
C                                                                         F00630
                  VS2 = (VJ-V0S1)**2
                  VS4 = VS2*VS2
                  SFAC = SFAC * 
     *                 (1.+FACTRS1*(HWSQ1/(VJ**2+(BETAS1*VS4)+HWSQ1)))  
c
                  VS2 = (VJ-V0S2)**2
                  SFAC = SFAC *
     *                 (1.+FACTRS2*(HWSQ2/(VS2+HWSQ2)))
c
                  VS2 = (VJ-V0S3)**2
                  VS4 = VS2*VS2
                  SFAC = SFAC *
     *                 (1.+FACTRS3*(HWSQ3/(VS2+(BETAS3*VS4)+HWSQ3))) 
C                                                                         
                  SH2O = SFAC * SH2O
c
C              ---------------------------------------------------------

               ENDIF
c
               C(J) = W1*(SH2O*RH2O)*XSELF
C                                                                         F00720
C              ---------------------------------------------------------
C              Radiation field                                            F00730
C                                                                         F00740
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   F00750
C              ---------------------------------------------------------

 20         CONTINUE                                                      F00760
C
c           Interpolate to total optical depth grid

            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F00770

         endif
C                                                                         F00780
C=======================================================================
C                             FOREIGN

C
         PFRGN = PATM-PH2O
         RFRGN = PFRGN*(T0/TAVE)

C--------------------------------------------------------------------

         V0F1 = 350.
         HWSQF1 = 200.**2
         BETAF1 = 5.e-09 
         FACTRF1 = -0.7
C
         V0F1a = 630.
         HWSQF1a = 65.**2
         BETAF1a = 2.e-08 
         FACTRF1a = +0.75
C
         V0F2 =1130.
         HWSQF2 = 330.**2
         BETAF2 = 8.E-11
         FACTRF2 = -0.97
C
         V0F3 = 1975.
         HWSQF3 = 250.**2
         BETAF3 = 5.E-06
         FACTRF3 = -0.65
C
C        ------------------------------------------------------------

C        Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
c
         if ((V2.gt.-20.0).and.(V1.lt.20000.)) then

            CALL FRN296 (V1C,V2C,DVC,NPTC,FH2O)         
C                                                               
            DO 24 J = 1, NPTC                                         
               VJ = V1C+DVC* REAL(J-1)                                    F00570
C
C              CORRECTION TO FOREIGN CONTINUUM
C
               VF2 = (VJ-V0F1)**2
               VF6 = VF2 * VF2 * VF2
               FSCAL = (1.+FACTRF1*(HWSQF1/(VF2+(BETAF1*VF6)+HWSQF1)))
C
               VF2 = (VJ-V0F1a)**2
               VF6 = VF2 * VF2 * VF2
               FSCAL = FSCAL* 
     *              (1.+FACTRF1a*(HWSQF1a/(VF2+(BETAF1a*VF6)+HWSQF1a)))
C
               VF2 = (VJ-V0F2)**2
               VF6 = VF2 * VF2 * VF2
               FSCAL = FSCAL* 
     *              (1.+FACTRF2*(HWSQF2/(VF2+(BETAF2*VF6)+HWSQF2)))
C
               VF2 = (VJ-V0F3)**2
               VF4 = VF2*VF2
               FSCAL = FSCAL* 
     *              (1.+FACTRF3*(HWSQF3/(VF2+BETAF3*VF4+HWSQF3)))
C     
               FH2O(J)=FH2O(J)*FSCAL
C     
               C(J) = W1*(FH2O(J)*RFRGN)*XFRGN
C                                          
C              ---------------------------------------------------------
C              Radiation field                                                  
C                                                                    
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)               
C              ---------------------------------------------------------
C
 24         CONTINUE                                                  
C
C           ------------------------------------------------------------
c
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)
C           ------------------------------------------------------------
C                                                                         F00780


C                                                                   
         endif

C=======================================================================


C     ********    CARBON DIOXIDE   ********                               F00790
C                                                                         F00800
C                                                                         F00810
C        Only calculate if V2 > -20. cm-1 and V1 <  10000. cm-1
c
         if ((V2.gt.-20.0).and.(V1.lt.10000.)) then
c
            WCO2 = WK(2)*RHOAVE*1.0E-20                                   F00820
C                                                                         F00830
            CALL FRNCO2 (V1C,V2C,DVC,NPTC,FCO2)                           F00840
            DO 30 J = 1, NPTC                                             F00850
               VJ = V1C+DVC* REAL(J-1)                                    F00860
               C(J) = FCO2(J)*WCO2*XCO2C                                  F00870

c****2.4.+++  The co2 continuum has been increased by a factor of 7. in the 
c                   nu2 band

               if (vj.gt.500 .and. vj.lt.900)  then
                  c(j) = 7.*c(j)
               endif


c****2.4.+++


C                                                                         F00880
C              Radiation field                                            F00890
C                                                                         F00900
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   F00910
 30         CONTINUE                                                      F00920
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F00930

         endif

C                                                                         F00940
C     ********    DIFFUSE OZONE  ********                                 F01180
C                                                                         F01190
C     Smoothing coefficients from 8920.0-9165.0 cm-1 and from
C     24570.0-24665.0 cm-1.  Data covers 9170.0-24565.0 cm-1
C     region.
C
         IF (V2.GT.8920.0.AND.V1.LE.24665.0) THEN                         F01200
            WO3 = WK(3)*1.0E-20                                           F01210
            CALL XO3CHP (V1C,V2C,DVC,NPTO3,CCH0,CCH1,CCH2)                F01220
C                                                                         F01230
            DT=TAVE-273.15
            DO 50 J = 1, NPTO3                                            F01240
               CCH0(J)=(CCH0(J)+(CCH1(J)+CCH2(J)*DT)*DT)*WO3*XO3CN
               VJ = V1C+DVC* REAL(J-1)                                    F01260
               IF (JRAD.EQ.1) CCH0(J) = CCH0(J)*RADFN(VJ,XKT)             F01270
 50         CONTINUE                                                      F01280
            CALL XINT (V1C,V2C,DVC,CCH0,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)   F01290
         ENDIF
C                                                                         F01300
         IF (V2.GT.27370..AND.V1.LT.40800.) THEN                          F01310
            WO3 = WK(3)*1.E-20                                            F01320
            TC = TAVE-273.15                                              F01330
            CALL O3HHT0 (V1C,V2C,DVC,NPTO3,C0)                            F01340
            CALL O3HHT1 (V1T1,V2T1,DVT1,NPT1,CT1)                         F01350
            CALL O3HHT2 (V1T2,V2T2,DVT2,NPT2,CT2)                         F01360
C                                                                         F01370
            DO 60 J = 1, NPTO3                                            F01380
               C(J) = C0(J)*WO3*XO3CN                                     F01390
C                                                                         F01400
               VJ = V1C+DVC* REAL(J-1)                                    F01410
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   F01420
               C(J) = C(J)*(1.+CT1(J)*TC+CT2(J)*TC*TC)                    F01430
 60         CONTINUE                                                      F01440
C
C           Save non-Hartley Huggins optical depth contribution to 
C           prevent double counting for wavenumber region beyond 
C           40800 cm-1.
C
            IF ((VJ.GT.40815.).AND.(V2.GT.40800)) THEN
               IFIX = (40800.-V1ABS)/DVABS+1.001
               DO 62 I=IFIX,NPTABS
                  ABSBSV(I) = ABSRB(I)
 62            CONTINUE
            ENDIF
C
C           Combine Hartley Huggins with previous optical depths
C
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F01450
C
C           If V2 > 40800 cm-1, replace points with previously
C           saved values (non-Hartley Huggins contribution)
C
            IF ((VJ.GT.40815.).AND.(V2.GT.40800)) THEN
               DO 64 I=IFIX,NPTABS
                  ABSRB(I) = ABSBSV(I)
 64            CONTINUE
            ENDIF
         ENDIF
C
C        If V2 > 40800 cm-1, add UV Hartley Huggins contribution
C
         IF (V2.GT.40800..AND.V1.LT.54000.) THEN
            WO3 = WK(3)                                                   F01470
            CALL O3HHUV (V1C,V2C,DVC,NPTO3,C0)                            F01480
C                                                                         F01490
            DO 70 J = 1, NPTO3                                            F01500
               C(J) = C0(J)*WO3*XO3CN                                     F01510
               VJ = V1C+DVC* REAL(J-1)                                    F01520
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   F01530
 70         CONTINUE                                                      F01540
C
C           Save non-Hartley Huggins UV optical depth contribution to
C           prevent double counting for wavenumber region before 
C           40800 cm-1.
C
            IF (V1.LT.40800) THEN
               IFIX = (40800.-V1ABS)/DVABS+1.001
               DO 72 I=1,IFIX-1
                  ABSBSV(I) = ABSRB(I)
 72            CONTINUE
            ENDIF
C
C           Combine UV Hartley Huggins with previous optical depths
C
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F01550
C
C           If V1 < 40800 cm-1, replace points with previously
C           saved values (non-Hartley Huggins UV contribution)
C
            IF (V1.LT.40800) THEN
               DO 74 I=1,IFIX-1
                  ABSRB(I) = ABSBSV(I)
 74            CONTINUE
            ENDIF
C                                                                         F01560
         ENDIF                                                            F01570
C                                                                         F01580


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
            
            rhofac = (Pave/P0)*(273./Tave)
c     
            tau_fac = Wk(7) * 1.e-20 * rhofac 
c
c           Wk(7) is the oxygen column amount in units of molec/cm2
c           rhofac is in units of amagats (air)
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
               C(J) = tau_fac * c0(J) * XO2CN
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
c        The units of these continua coefficients are  1 / (amagat_O2*amagat_air)
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
            WO2 = (WK(7)/xlosmt) * ((pave/1013.)*(273./tave))
            CHIO2 = (WK(7)*1.E-20)/WTOT 
            CHIN2 = (WBROAD*1.E-20)/WTOT 
            ADJFAC = (CHIO2+0.3*CHIN2)/0.446
            ADJWO2 = ADJFAC * WO2
c
            CALL O2INF1 (V1C,V2C,DVC,NPTC,C0)                     
c
            DO 92 J = 1, NPTC                                                
               C(J) = C0(J)*ADJWO2*XO2CN
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
            WO2 = RHOAVE*WK(7)*1.e-20
            CHIO2 = (WK(7)*1.E-20)/WTOT 
            ADJFAC = CHIO2/0.209
            ADJWO2 = ADJFAC * WO2
c
            DO 93 J = 1, NPTC                                                
               C(J) = C0(J)*ADJWO2*XO2CN
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
C        O2 continuum formulated by Greenblatt et al. over the spectral region
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
            WO2 = (WK(7)*1.e-20) * ((pave/1013.)*(273./tave))
            CHIO2 = (WK(7)*1.E-20)/WTOT 
            ADJWO2 = chio2 * WO2
c
            CALL O2_vis (V1C,V2C,DVC,NPTC,C0)                     
c
            DO 94 J = 1, NPTC                                                
               C(J) = C0(J)*ADJWO2*XO2CN
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

            CALL O2HERZ (V1C,V2C,DVC,NPTC,C0,TAVE,PAVE)                   F01820
            DO 90 J = 1, NPTC                                             F01830
               C(J) = C0(J)*WK(7)*XO2CN                                   F01840
               VJ = V1C+DVC* REAL(J-1)                                    F01850
C                                                                         F01860
C              Radiation field                                            F01870
C                                                                         F01880
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   F01890
 90         CONTINUE                                                      F01900
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F01910

         endif
c
C     *********************  NITROGEN CONTINUA  ********************
c
         IF (NMOL.GE.22) THEN                                             F01000
            WN2 = WK(22)                                                  F01010
         ELSE                                                             F01020
            WN2 = WBROAD                                                  F01030
         ENDIF                                                            F01040
C                                                                         F00460
C                                                                         F00940
C     ******** NITROGEN COLLISION INDUCED PURE ROTATION BAND  ********
C
C        Model used:
C         Borysow, A, and L. Frommhold, "Collision-induced
C            rototranslational absorption spectra of N2-N2
C            pairs for temperatures from 50 to 300 K", The
C            Astrophysical Journal, 311, 1043-1057, 1986.
C
C                                                                         F00960
C        THIS NITROGEN CONTINUUM IS IN UNITS OF 1./(CM AMAGAT^2)
C
C        Only calculate if V2 > -10. cm-1 and V1 <  350. cm-1
c
         if ((V2.gt.-10.0).and.(V1.lt.350.)) then
c
C           The following puts WXN2 units in 1./(CM AMAGAT)
C
            WXN2 = WN2/XLOSMT
C                                                                         F00480
C           RHOFAC units are AMAGATS
C
            RHOFAC = ((WN2*1.e-20)/WTOT)*(PAVE/P0)*(273./TAVE)
C     
            TFAC = (TAVE-T0)/(220.-T0)

            CALL N2R296 (V1C,V2C,DVC,NPTC,C0)
            CALL N2R220 (V1C,V2C,DVC,NPTC,C1)
C
            DO 40 J = 1, NPTC                                             F01080
               VJ = V1C+DVC* REAL(J-1)                                    F01090
C
               C(J) = 0.
               IF (C0(J).GT.0. .AND. C1(J).GT.0.)
     *            C(J) = (WXN2*RHOFAC*C0(J)*(C1(J)/C0(J))**TFAC)*XN2CN
C                                                                         F01110
C              Radiation field                                            F01120
C                                                                         F01130
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)                   F01140
 40         CONTINUE                                                      F01150
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F01160

         endif
C                                                                         F01170
C                                                                         F00940
C        ********    NITROGEN COLLISION INDUCED FUNDAMENTAL ********      F00950
C                                                                         F00960
c        version_1 of the Nitrogen Collision Induced Fundamental
c
c        Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and 
c        J._M. Hartmann, Infrared collision-induced absorption by 
c        N2 near 4.3 microns for atmospheric applications: 
c        Measurements and emprirical modeling, Appl. Optics, 35, 
c        5911-5917, (1996).
c
C        Only calculate if V2 > 2085. cm-1 and V1 <  2670. cm-1
C
         if ((V2.gt.2085.0).and.(V1.lt.2670.)) then

            rhofac = (Pave/P0)*(273./Tave)
c     
            tau_fac = Wn2 * 1.e-20 * rhofac 
c
c           Wn2 is in units of molec/cm2
c           rhofac is in units of amagats (air)
c
c           The temperature correction is done in subroutine n2_ver_1:
c
            call n2_ver_1 (v1c,v2c,dvc,nptc,c0,tave)
C
c           c0 are the nitrogen absorption coefficients at 
c           temperature tave 
c              - these absorption coefficients are in units of
c                   [(cm^2/molec) 10^20)]/(cm-1  amagat) 
c              - cm-1 in the denominator arises through the removal
c                   of the radiation field
c              - for this case, an amagat is interpreted as one
c                   loshmidt of air (273K)
c
            DO 45 J = 1, NPTC
               VJ = V1C+DVC* REAL(J-1)
               C(J) = tau_fac * c0(J) * XN2CN
C
C              Radiation field
C
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
 45         CONTINUE
C
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)
c
         endif
c
C     ********** Rayleigh Scattering calculation **********
c
c     The effects of Rayleigh scattering are also included in module ! sac  06/06/02
c     lbllow.f with the aerosol/cloud properties.  In the case that 
c     that lbllow.f is selected (iaersl .ne. 0), the decision has been
c     made to include this effect with the scattering processes,
c     even though it is molecular in nature.  Otherwise the effects
c     of Rayleigh scattering are included here.
c
c     The formulation, adopted from MODTRAN_3.5 (using approximation
c     of Shettle et al., (Appl Opt, 2873-4, 1980) with depolarization
c     = 0.0279, output in km-1 for T=273K & P=1 ATM) is as follows:
C
c     The rayleigh extinction coefficient (scattering out of the direct
c     beam), ray_ext, can be defined as
c
c         ray_ext = (vrayleigh**4/(9.38076E18-1.08426E09*vrayleigh**2))
c     *        *wmol_tot*conv_cm2mol
c
c     where vrayleigh is the wavenumber value, wmol_tot is the total
c     molecular amount in the layer, and conv_cm2mol is the conversion
c     factor derived by multiplying air density (2.68675E19 mol/cm3)
c     at 273 K with the number of km per cm (1.e-5 km/cm).
c
c     For numerical purposes, the layer amount of all molecules is
c     calculated above as WTOT, which has been scaled by 1.e-20. We
c     have taken 1.e20 out of the air density in the denominator
c     as well. In addition, a factor of 10,000 (cm-1) has been 
c     divided out of vrayleigh. Finally, the radiation field is
c     excluded, so xvrayleigh**4 is replaced by xvrayleigh**3. When
c     JRAD=1, the radiation field is put in by multiplying the
c     absorption by xvrayleigh.
c
c     Rayleigh scattering in the direct beam is only calculated for
c     model runs > 3100 cm-1.
c
         If (iaersl .eq.0 .and. v2.ge.3100.) then
c
c        Thus the current formulation is

            conv_cm2mol = 1./(2.68675e-1*1.e5)
         
            do 95 i=1,nptabs
               vrayleigh = v1abs+(i-1)*dvabs
               xvrayleigh = vrayleigh/1.e4
          ray_ext = (xvrayleigh**3/(9.38076E2-10.8426*xvrayleigh**2))
     *           *wtot*conv_cm2mol*XRAYL

C           Radiation field

               IF (JRAD.EQ.1) ray_ext = ray_ext*xvrayleigh

            absrb(i) = absrb(i)+ray_ext
 95      continue
c
      endif
C                                                                         F01920
 100     continue
      RETURN                                                              F01930
C                                                                         F01940
 900  FORMAT (/,'0    *********************************************',/,
     *          '     *      BYPASS O2 CONTINUUM TO HERZBERG      *',/,
     *          '     *       AS A RESULT of TAVE > 350. K        *',/,
     *          '     *********************************************',/)
C
      END                                                                 F01950
C
C     --------------------------------------------------------------
C
C
      SUBROUTINE PRCNTM                                                   A10270
C                                                                         A10280
C     THIS SUBROUTINE PRINTS THE CONTINUUM INFORMATION TO FILE IPR        A10290
C                                                                         A10300
      CHARACTER*51 CINFO1(2,9),CINFO2(2,11),CINFO3(2,9),CINFO4(2,12)
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A10310
     *              NLNGTH,KFILE,KPANEL,LINFIL,NDFLE,IAFIL,IEXFIL,        A10320
     *              NLTEFL,LNFIL4,LNGTH4                                  A10330
      COMMON /CNTPR/ CINFO1,CINFO2,CINFO3,CINFO4
C                                                                         A10340
      WRITE (IPR,910) ((CINFO1(I,J),I=1,2),J=1,9)
      WRITE (IPR,910) ((CINFO2(I,J),I=1,2),J=1,11)
      WRITE (IPR,910) ((CINFO3(I,J),I=1,2),J=1,9)
      WRITE (IPR,910) ((CINFO4(I,J),I=1,2),J=1,12)
C                                                                         A10360
      RETURN                                                              A10370
C                                                                         A10380
 910  FORMAT (2A51)
C                                                                         A10580
      END                                                                 A10590
C
C     --------------------------------------------------------------
      BLOCK DATA CNTINF
C
C     Continuum information for output to TAPE6 in SUBROUTINE PRCNTM
C
      CHARACTER*51 CINFO1(2,9),CINFO2(2,11),CINFO3(2,9),CINFO4(2,12)
      COMMON /CNTPR/ CINFO1,CINFO2,CINFO3,CINFO4
C
      DATA CINFO1/
c           123456789-123456789-123456789-123456789-123456789-1
     1     '                                                   ',
     2     '                                                   ',
     3     '                                                   ',
     4     '                                                   ',
     5     '0  *****  CONTINUA ckd_2.4.2                       ',
     6     '                                                   ',
     7     '                                                   ',
     8     '                                                   ',
     9     '                     H2O   SELF  (T)      0 - 20000',
     *     ' CM-1    ckd_2.4                                   ',
     1     '                           AIR   (T)      0 - 20000',
     2     ' CM-1    ckd_2.4                                   ',
     3     '                     CO2   AIR            0 - 20000',
     4     ' CM-1    co2 nu2 increased * 7         (July 2002) ',
     5     '                     N2    SELF           0 -   350',
     6     ' CM-1    BORYSOW FROMMHOLD                         ',
     7     '                           AIR         2085 -  2670',
     8     ' CM-1                                 (March 1998) ' /
C
      DATA CINFO2/
     *     '                     O2    AIR   (T)   1340 -  1850',
     *     ' CM-1                                 (March 1998) ',
     *     '                           O2/N2       7550 -  8486',
     *     ' CM-1                              (February 2000) ',
     *     '                           AIR         9100 - 11000',
     *     ' CM-1                                (August 1999) ',
     *     '                           O2         15000 - 29870',
     *     ' CM-1                                  (May  2000) ',
     *     '                           O2/N2      36000 -  >>>>',
     *     ' CM-1    HERZBERG                                  ',
     *     '                     O3    AIR         9170 - 24565',
     *     ' CM-1    CHAPPUIS / WULF                           ',
     *     '                                 (T)  27370 - 40800',
     *     ' CM-1    HARTLEY HUGGINS                           ',
     *     '                                      40800 - 54000',
     *     ' CM-1    HARTLEY HUGGINS                           ',
     *     6*'                                                 ' /
C
      DATA CINFO3/
     *     '  H2O SELF HAS BEEN REDUCED IN THE 800-1200 CM-1 RE',
     *     'GION                                 (01 SEPT 1985)',
     *     '  03       TEMPERATURE DEPENDENCE HAS BEEN CORRECTE',
     *     'D                                     (01 MAY 1987)',
     *     '  02       (1390-1760) HAS BEEN REDUCED (FACTOR = 0',
     *     '.78)                                (07 MARCH 1990)',
     *     '  H2O SELF HAS BEEN REDUCED IN THE 1100-1500 CM-1 R',
     *     'EGION                               (01 APRIL 1993)',
     *     '  H2O FOREIGN HAS BEEN REDUCED AT ~1300 CM-1 AND IN',
     *     ' ALL THE WINDOW REGIONS             (01 APRIL 1993)',
     *     '  H2O SELF HAS BEEN MODIFIED IN THE 700-1500 CM-1 R',
     *     'EGION                                 (01 MAY 1994)',
     *     '  H2O FOREIGN HAS BEEN MODIFIED IN THE ENTIRE 1200-',
     *     '2200 CM-1 SPECTRAL RANGE              (01 MAY 1994)',
     *     '  H2O SELF HAS BEEN INCREASED 30% IN THE MICROWAVE ',
     *     'REGION                                (09 FEB 1996)',
     *     '  N2 COLLISION INDUCED PURE ROTATION BAND ADDED    ',
     *     '                                      (09 FEB 1996)'/

      DATA CINFO4/
     *     '  O3 CHAPPUIS CHANGED TO VALUES FROM MODTRAN3      ',
     *     '                                      (09 FEB 1996)',
     *     '  INTERPOLATION EFFECTS BETWEEN HARTLEY HUGGINS DAT',
     *     'A AROUND 40800 CM-1                  (18 SEPT 1996)',
     *     '  O2 COLLISION INDUCED FUNDAMENTAL BAND (1340-1850 ',
     *     'CM-1) CHANGED TO THIBAULT ET AL., 1996 (MARCH 1998)',
     *     '  N2 COLLISION INDUCED FUNDAMENTAL BAND (2085-2670 ',
     *     'CM-1) CHANGED TO LAFFERTY ET AL., 1996 (MARCH 1998)',
     *     '  H2O FOREIGN HAS BEEN MODIFIED IN THE 0-800 CM-1 A',
     *     'ND 1200-2200 CM-1 RANGE BASED ON      (04 JUN 1999)',
     *     '           AERI-ER FROM SHEBA, TOBIN ET AL., 1998  ',
     *     '                                                   ',
     *     '  H2O SELF HAS BEEN INCREASED IN THE 0-200 CM-1 REG',
     *     'ION                                   (04 JUN 1999)',
     *     '  O2 COLLISION INDUCED BAND HAS BEEN ADDED  (9100 -',
     *     ' 11000 CM-1); MLAWER ET AL. 1998      (16 AUG 1999)',
     *     '  O2 COLLISION INDUCED BAND HAS BEEN CHANGED (7555 ',
     *     '- 8486 CM-1); MATE ET AL. 1999        (10 FEB 2000)',
     *     '  O2 COLLISION INDUCED BANDS HAVE BEEN ADDED (15000',
     *     ' - 29870 CM-1); GREENBLATT ET AL. 1990 (8 MAY 2000)',
     *     '  The nu2 CO2 increased by a factor of 7: based on ',
     *     'U  Wisc. AFWEX and AERI_xr at ARM NSA   (July 2002)',
     *     '  -------------------------------------------------',
     *     '------------------------------------------         '/
C
      END
C
C
C
C
      SUBROUTINE SL296 (V1C,V2C,DVC,NPTC,C)                               F02220
C                                                                         F02230
      IMPLICIT REAL*8           (V)                                     ! F02240
C                                                                         F02250
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F02260
      COMMON /SH2O/ V1S,V2S,DVS,NPTS,S(2003)                              F02270
      DIMENSION C(*)                                                      F02280
C                                                                         F02290
      DVC = DVS                                                           F02300
      V1C = V1ABS-DVC                                                     F02310
      V2C = V2ABS+DVC                                                     F02320
C                                                                         F02330
      I1 = (V1C-V1S)/DVS                                                  F02340
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F02410
         I = I1+J                                                         F02420
         C(J) = 0.                                                        F02430
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F02440
         C(J) = S(I)                                                      F02450
   10 CONTINUE                                                            F02460
C                                                                         F02470
      RETURN                                                              F02480
C                                                                         F02490
      END                                                                 F02500
C
C     --------------------------------------------------------------
C
      BLOCK DATA BS296                                                    F02510
C                                                                         F02520
      IMPLICIT REAL*8           (V)                                     ! F02530
C                                                                         F02540
C               06/28/82                                                  F02550
C               UNITS OF (CM**3/MOL) * 1.E-20                             F02560
C                                                                         F02570
      COMMON /SH2O/ V1,V2,DV,NPT,                                         F02580
     *              S0000( 2),S0001(50),S0051(50),S0101(50),S0151(50),    F02590
     *              S0201(50),S0251(50),S0301(50),S0351(50),S0401(50),    F02600
     *              S0451(50),S0501(50),S0551(50),S0601(50),S0651(50),    F02610
     *              S0701(50),S0751(50),S0801(50),S0851(50),S0901(50),    F02620
     *              S0951(50),S1001(50),S1051(50),S1101(50),S1151(50),    F02630
     *              S1201(50),S1251(50),S1301(50),S1351(50),S1401(50),    F02640
     *              S1451(50),S1501(50),S1551(50),S1601(50),S1651(50),    F02650
     *              S1701(50),S1751(50),S1801(50),S1851(50),S1901(50),    F02660
     *              S1951(50),S2001(1)                                    F02670
C                                                                         F02680
       DATA V1,V2,DV,NPT /  -20.0, 20000.0, 10.0, 2003/                   F02690
C                                                                         F02700
      DATA S0000/                                                         F02710
     *     1.1109E-01 ,1.0573E-01/                                        F02720
      DATA S0001/                                                         F02730
     *     1.0162E-01, 1.0573E-01, 1.1109E-01, 1.2574E-01, 1.3499E-01,    F02740
     *     1.4327E-01, 1.5065E-01, 1.5164E-01, 1.5022E-01, 1.3677E-01,    F02750
     *     1.3115E-01, 1.2253E-01, 1.1271E-01, 1.0070E-01, 8.7495E-02,    F02760
     *     8.0118E-02, 6.9940E-02, 6.2034E-02, 5.6051E-02, 4.7663E-02,    F02770
     *     4.2450E-02, 3.6690E-02, 3.3441E-02, 3.0711E-02, 2.5205E-02,    F02780
     *     2.2113E-02, 1.8880E-02, 1.6653E-02, 1.4626E-02, 1.2065E-02,    F02790
     *     1.0709E-02, 9.1783E-03, 7.7274E-03, 6.7302E-03, 5.6164E-03,    F02800
     *     4.9089E-03, 4.1497E-03, 3.5823E-03, 3.1124E-03, 2.6414E-03,    F02810
     *     2.3167E-03, 2.0156E-03, 1.7829E-03, 1.5666E-03, 1.3928E-03,    F02820
     *     1.2338E-03, 1.0932E-03, 9.7939E-04, 8.8241E-04, 7.9173E-04/    F02830
      DATA S0051/                                                         F02840
     *     7.1296E-04, 6.4179E-04, 5.8031E-04, 5.2647E-04, 4.7762E-04,    F02850
     *     4.3349E-04, 3.9355E-04, 3.5887E-04, 3.2723E-04, 2.9919E-04,    F02860
     *     2.7363E-04, 2.5013E-04, 2.2876E-04, 2.0924E-04, 1.9193E-04,    F02870
     *     1.7618E-04, 1.6188E-04, 1.4891E-04, 1.3717E-04, 1.2647E-04,    F02880
     *     1.1671E-04, 1.0786E-04, 9.9785E-05, 9.2350E-05, 8.5539E-05,    F02890
     *     7.9377E-05, 7.3781E-05, 6.8677E-05, 6.3993E-05, 5.9705E-05,    F02900
     *     5.5788E-05, 5.2196E-05, 4.8899E-05, 4.5865E-05, 4.3079E-05,    F02910
     *     4.0526E-05, 3.8182E-05, 3.6025E-05, 3.4038E-05, 3.2203E-05,    F02920
     *     3.0511E-05, 2.8949E-05, 2.7505E-05, 2.6170E-05, 2.4933E-05,    F02930
     *     2.3786E-05, 2.2722E-05, 2.1736E-05, 2.0819E-05, 1.9968E-05/    F02940
      DATA S0101/                                                         F02950
     *     1.9178E-05, 1.8442E-05, 1.7760E-05, 1.7127E-05, 1.6541E-05,    F02960
     *     1.5997E-05, 1.5495E-05, 1.5034E-05, 1.4614E-05, 1.4230E-05,    F02970
     *     1.3883E-05, 1.3578E-05, 1.3304E-05, 1.3069E-05, 1.2876E-05,    F02980
     *     1.2732E-05, 1.2626E-05, 1.2556E-05, 1.2544E-05, 1.2604E-05,    F02990
     *     1.2719E-05, 1.2883E-05, 1.3164E-05, 1.3581E-05, 1.4187E-05,    F03000
     *     1.4866E-05, 1.5669E-05, 1.6717E-05, 1.8148E-05, 2.0268E-05,    F03010
     *     2.2456E-05, 2.5582E-05, 2.9183E-05, 3.3612E-05, 3.9996E-05,    F03020
     *     4.6829E-05, 5.5055E-05, 6.5897E-05, 7.5360E-05, 8.7213E-05,    F03030
     *     1.0046E-04, 1.1496E-04, 1.2943E-04, 1.5049E-04, 1.6973E-04,    F03040
     *     1.8711E-04, 2.0286E-04, 2.2823E-04, 2.6780E-04, 2.8766E-04/    F03050
      DATA S0151/                                                         F03060
     *     3.1164E-04, 3.3640E-04, 3.6884E-04, 3.9159E-04, 3.8712E-04,    F03070
     *     3.7433E-04, 3.4503E-04, 3.1003E-04, 2.8027E-04, 2.5253E-04,    F03080
     *     2.3408E-04, 2.2836E-04, 2.4442E-04, 2.7521E-04, 2.9048E-04,    F03090
     *     3.0489E-04, 3.2646E-04, 3.3880E-04, 3.3492E-04, 3.0987E-04,    F03100
     *     2.9482E-04, 2.8711E-04, 2.6068E-04, 2.2683E-04, 1.9996E-04,    F03110
     *     1.7788E-04, 1.6101E-04, 1.3911E-04, 1.2013E-04, 1.0544E-04,    F03120
     *     9.4224E-05, 8.1256E-05, 7.3667E-05, 6.2233E-05, 5.5906E-05,    F03130
     *     5.1619E-05, 4.5140E-05, 4.0273E-05, 3.3268E-05, 3.0258E-05,    F03140
     *     2.6440E-05, 2.3103E-05, 2.0749E-05, 1.8258E-05, 1.6459E-05,    F03150
     *     1.4097E-05, 1.2052E-05, 1.0759E-05, 9.1400E-06, 8.1432E-06/    F03160
      DATA S0201/                                                         F03170
     *     7.1460E-06, 6.4006E-06, 5.6995E-06, 4.9372E-06, 4.4455E-06,    F03180
     *     3.9033E-06, 3.4740E-06, 3.1269E-06, 2.8059E-06, 2.5558E-06,    F03190
     *     2.2919E-06, 2.0846E-06, 1.8983E-06, 1.7329E-06, 1.5929E-06,    F03200
     *     1.4631E-06, 1.3513E-06, 1.2461E-06, 1.1519E-06, 1.0682E-06,    F03210
     *     9.9256E-07, 9.2505E-07, 8.6367E-07, 8.0857E-07, 7.5674E-07,    F03220
     *     7.0934E-07, 6.6580E-07, 6.2580E-07, 5.8853E-07, 5.5333E-07,    F03230
     *     5.2143E-07, 4.9169E-07, 4.6431E-07, 4.3898E-07, 4.1564E-07,    F03240
     *     3.9405E-07, 3.7403E-07, 3.5544E-07, 3.3819E-07, 3.2212E-07,    F03250
     *     3.0714E-07, 2.9313E-07, 2.8003E-07, 2.6777E-07, 2.5628E-07,    F03260
     *     2.4551E-07, 2.3540E-07, 2.2591E-07, 2.1701E-07, 2.0866E-07/    F03270
      DATA S0251/                                                         F03280
     *     2.0082E-07, 1.9349E-07, 1.8665E-07, 1.8027E-07, 1.7439E-07,    F03290
     *     1.6894E-07, 1.6400E-07, 1.5953E-07, 1.5557E-07, 1.5195E-07,    F03300
     *     1.4888E-07, 1.4603E-07, 1.4337E-07, 1.4093E-07, 1.3828E-07,    F03310
     *     1.3569E-07, 1.3270E-07, 1.2984E-07, 1.2714E-07, 1.2541E-07,    F03320
     *     1.2399E-07, 1.2102E-07, 1.1878E-07, 1.1728E-07, 1.1644E-07,    F03330
     *     1.1491E-07, 1.1305E-07, 1.1235E-07, 1.1228E-07, 1.1224E-07,    F03340
     *     1.1191E-07, 1.1151E-07, 1.1098E-07, 1.1068E-07, 1.1109E-07,    F03350
     *     1.1213E-07, 1.1431E-07, 1.1826E-07, 1.2322E-07, 1.3025E-07,    F03360
     *     1.4066E-07, 1.5657E-07, 1.7214E-07, 1.9449E-07, 2.2662E-07,    F03370
     *     2.6953E-07, 3.1723E-07, 3.7028E-07, 4.4482E-07, 5.3852E-07/    F03380
      DATA S0301/                                                         F03390
     *     6.2639E-07, 7.2175E-07, 7.7626E-07, 8.7248E-07, 9.6759E-07,    F03400
     *     1.0102E-06, 1.0620E-06, 1.1201E-06, 1.2107E-06, 1.2998E-06,    F03410
     *     1.3130E-06, 1.2856E-06, 1.2350E-06, 1.1489E-06, 1.0819E-06,    F03420
     *     1.0120E-06, 9.4795E-07, 9.2858E-07, 9.8060E-07, 1.0999E-06,    F03430
     *     1.1967E-06, 1.2672E-06, 1.3418E-06, 1.3864E-06, 1.4330E-06,    F03440
     *     1.4592E-06, 1.4598E-06, 1.4774E-06, 1.4726E-06, 1.4820E-06,    F03450
     *     1.5050E-06, 1.4984E-06, 1.5181E-06, 1.5888E-06, 1.6850E-06,    F03460
     *     1.7690E-06, 1.9277E-06, 2.1107E-06, 2.3068E-06, 2.5347E-06,    F03470
     *     2.8069E-06, 3.1345E-06, 3.5822E-06, 3.9051E-06, 4.3422E-06,    F03480
     *     4.8704E-06, 5.5351E-06, 6.3454E-06, 7.2690E-06, 8.2974E-06/    F03490
      DATA S0351/                                                         F03500
     *     9.7609E-06, 1.1237E-05, 1.3187E-05, 1.5548E-05, 1.8784E-05,    F03510
     *     2.1694E-05, 2.5487E-05, 3.0092E-05, 3.5385E-05, 4.2764E-05,    F03520
     *     4.9313E-05, 5.5800E-05, 6.2968E-05, 7.1060E-05, 7.7699E-05,    F03530
     *     8.7216E-05, 8.9335E-05, 9.2151E-05, 9.2779E-05, 9.4643E-05,    F03540
     *     9.7978E-05, 1.0008E-04, 1.0702E-04, 1.1026E-04, 1.0828E-04,    F03550
     *     1.0550E-04, 1.0432E-04, 1.0428E-04, 9.8980E-05, 9.4992E-05,    F03560
     *     9.5159E-05, 1.0058E-04, 1.0738E-04, 1.1550E-04, 1.1229E-04,    F03570
     *     1.0596E-04, 1.0062E-04, 9.1742E-05, 8.4492E-05, 6.8099E-05,    F03580
     *     5.6295E-05, 4.6502E-05, 3.8071E-05, 3.0721E-05, 2.3297E-05,    F03590
     *     1.8688E-05, 1.4830E-05, 1.2049E-05, 9.6754E-06, 7.9192E-06/    F03600
      DATA S0401/                                                         F03610
     *     6.6673E-06, 5.6468E-06, 4.8904E-06, 4.2289E-06, 3.6880E-06,    F03620
     *     3.2396E-06, 2.8525E-06, 2.5363E-06, 2.2431E-06, 1.9949E-06,    F03630
     *     1.7931E-06, 1.6164E-06, 1.4431E-06, 1.2997E-06, 1.1559E-06,    F03640
     *     1.0404E-06, 9.4300E-07, 8.4597E-07, 7.6133E-07, 6.8623E-07,    F03650
     *     6.2137E-07, 5.6345E-07, 5.1076E-07, 4.6246E-07, 4.1906E-07,    F03660
     *     3.8063E-07, 3.4610E-07, 3.1554E-07, 2.8795E-07, 2.6252E-07,    F03670
     *     2.3967E-07, 2.1901E-07, 2.0052E-07, 1.8384E-07, 1.6847E-07,    F03680
     *     1.5459E-07, 1.4204E-07, 1.3068E-07, 1.2036E-07, 1.1095E-07,    F03690
     *     1.0237E-07, 9.4592E-08, 8.7530E-08, 8.1121E-08, 7.5282E-08,    F03700
     *     6.9985E-08, 6.5189E-08, 6.0874E-08, 5.6989E-08, 5.3530E-08/    F03710
      DATA S0451/                                                         F03720
     *     5.0418E-08, 4.7745E-08, 4.5367E-08, 4.3253E-08, 4.1309E-08,    F03730
     *     3.9695E-08, 3.8094E-08, 3.6482E-08, 3.4897E-08, 3.3500E-08,    F03740
     *     3.2302E-08, 3.0854E-08, 2.9698E-08, 2.8567E-08, 2.7600E-08,    F03750
     *     2.6746E-08, 2.5982E-08, 2.5510E-08, 2.5121E-08, 2.4922E-08,    F03760
     *     2.4909E-08, 2.5013E-08, 2.5216E-08, 2.5589E-08, 2.6049E-08,    F03770
     *     2.6451E-08, 2.6978E-08, 2.7687E-08, 2.8600E-08, 2.9643E-08,    F03780
     *     3.0701E-08, 3.2058E-08, 3.3695E-08, 3.5558E-08, 3.7634E-08,    F03790
     *     3.9875E-08, 4.2458E-08, 4.5480E-08, 4.8858E-08, 5.2599E-08,    F03800
     *     5.7030E-08, 6.2067E-08, 6.7911E-08, 7.4579E-08, 8.1902E-08,    F03810
     *     8.9978E-08, 9.9870E-08, 1.1102E-07, 1.2343E-07, 1.3732E-07/    F03820
      DATA S0501/                                                         F03830
     *     1.5394E-07, 1.7318E-07, 1.9383E-07, 2.1819E-07, 2.4666E-07,    F03840
     *     2.8109E-07, 3.2236E-07, 3.7760E-07, 4.4417E-07, 5.2422E-07,    F03850
     *     6.1941E-07, 7.4897E-07, 9.2041E-07, 1.1574E-06, 1.4126E-06,    F03860
     *     1.7197E-06, 2.1399E-06, 2.6266E-06, 3.3424E-06, 3.8418E-06,    F03870
     *     4.5140E-06, 5.0653E-06, 5.8485E-06, 6.5856E-06, 6.8937E-06,    F03880
     *     6.9121E-06, 6.9005E-06, 6.9861E-06, 6.8200E-06, 6.6089E-06,    F03890
     *     6.5809E-06, 7.3496E-06, 8.0311E-06, 8.3186E-06, 8.4260E-06,    F03900
     *     9.0644E-06, 9.4965E-06, 9.4909E-06, 9.0160E-06, 9.1494E-06,    F03910
     *     9.3629E-06, 9.5944E-06, 9.5459E-06, 8.9919E-06, 8.6040E-06,    F03920
     *     7.8613E-06, 7.1567E-06, 6.2677E-06, 5.1899E-06, 4.4188E-06/    F03930
      DATA S0551/                                                         F03940
     *     3.7167E-06, 3.0636E-06, 2.5573E-06, 2.0317E-06, 1.6371E-06,    F03950
     *     1.3257E-06, 1.0928E-06, 8.9986E-07, 7.4653E-07, 6.1111E-07,    F03960
     *     5.1395E-07, 4.3500E-07, 3.7584E-07, 3.2633E-07, 2.8413E-07,    F03970
     *     2.4723E-07, 2.1709E-07, 1.9294E-07, 1.7258E-07, 1.5492E-07,    F03980
     *     1.3820E-07, 1.2389E-07, 1.1189E-07, 1.0046E-07, 9.0832E-08,    F03990
     *     8.2764E-08, 7.4191E-08, 6.7085E-08, 6.0708E-08, 5.4963E-08,    F04000
     *     4.9851E-08, 4.5044E-08, 4.0916E-08, 3.7220E-08, 3.3678E-08,    F04010
     *     3.0663E-08, 2.7979E-08, 2.5495E-08, 2.3286E-08, 2.1233E-08,    F04020
     *     1.9409E-08, 1.7770E-08, 1.6260E-08, 1.4885E-08, 1.3674E-08,    F04030
     *     1.2543E-08, 1.1551E-08, 1.0655E-08, 9.8585E-09, 9.1398E-09/    F04040
      DATA S0601/                                                         F04050
     *     8.4806E-09, 7.8899E-09, 7.3547E-09, 6.8670E-09, 6.4131E-09,    F04060
     *     5.9930E-09, 5.6096E-09, 5.2592E-09, 4.9352E-09, 4.6354E-09,    F04070
     *     4.3722E-09, 4.1250E-09, 3.9081E-09, 3.7118E-09, 3.5372E-09,    F04080
     *     3.3862E-09, 3.2499E-09, 3.1324E-09, 3.0313E-09, 2.9438E-09,    F04090
     *     2.8686E-09, 2.8050E-09, 2.7545E-09, 2.7149E-09, 2.6907E-09,    F04100
     *     2.6724E-09, 2.6649E-09, 2.6642E-09, 2.6725E-09, 2.6871E-09,    F04110
     *     2.7056E-09, 2.7357E-09, 2.7781E-09, 2.8358E-09, 2.9067E-09,    F04120
     *     2.9952E-09, 3.1020E-09, 3.2253E-09, 3.3647E-09, 3.5232E-09,    F04130
     *     3.7037E-09, 3.9076E-09, 4.1385E-09, 4.3927E-09, 4.6861E-09,    F04140
     *     5.0238E-09, 5.4027E-09, 5.8303E-09, 6.3208E-09, 6.8878E-09/    F04150
      DATA S0651/                                                         F04160
     *     7.5419E-09, 8.3130E-09, 9.1952E-09, 1.0228E-08, 1.1386E-08,    F04170
     *     1.2792E-08, 1.4521E-08, 1.6437E-08, 1.8674E-08, 2.1160E-08,    F04180
     *     2.4506E-08, 2.8113E-08, 3.2636E-08, 3.7355E-08, 4.2234E-08,    F04190
     *     4.9282E-08, 5.7358E-08, 6.6743E-08, 7.8821E-08, 9.4264E-08,    F04200
     *     1.1542E-07, 1.3684E-07, 1.6337E-07, 2.0056E-07, 2.3252E-07,    F04210
     *     2.6127E-07, 2.9211E-07, 3.3804E-07, 3.7397E-07, 3.8205E-07,    F04220
     *     3.8810E-07, 3.9499E-07, 3.9508E-07, 3.7652E-07, 3.5859E-07,    F04230
     *     3.6198E-07, 3.7871E-07, 4.0925E-07, 4.2717E-07, 4.8241E-07,    F04240
     *     5.2008E-07, 5.6530E-07, 5.9531E-07, 6.1994E-07, 6.5080E-07,    F04250
     *     6.6355E-07, 6.9193E-07, 6.9930E-07, 7.3058E-07, 7.4678E-07/    F04260
      DATA S0701/                                                         F04270
     *     7.9193E-07, 8.3627E-07, 9.1267E-07, 1.0021E-06, 1.1218E-06,    F04280
     *     1.2899E-06, 1.4447E-06, 1.7268E-06, 2.0025E-06, 2.3139E-06,    F04290
     *     2.5599E-06, 2.8920E-06, 3.3059E-06, 3.5425E-06, 3.9522E-06,    F04300
     *     4.0551E-06, 4.2818E-06, 4.2892E-06, 4.4210E-06, 4.5614E-06,    F04310
     *     4.6739E-06, 4.9482E-06, 5.1118E-06, 5.0986E-06, 4.9417E-06,    F04320
     *     4.9022E-06, 4.8449E-06, 4.8694E-06, 4.8111E-06, 4.9378E-06,    F04330
     *     5.3231E-06, 5.7362E-06, 6.2350E-06, 6.0951E-06, 5.7281E-06,    F04340
     *     5.4585E-06, 4.9032E-06, 4.3009E-06, 3.4776E-06, 2.8108E-06,    F04350
     *     2.2993E-06, 1.7999E-06, 1.3870E-06, 1.0750E-06, 8.5191E-07,    F04360
     *     6.7951E-07, 5.5336E-07, 4.6439E-07, 4.0243E-07, 3.5368E-07/    F04370
      DATA S0751/                                                         F04380
     *     3.1427E-07, 2.7775E-07, 2.4486E-07, 2.1788E-07, 1.9249E-07,    F04390
     *     1.7162E-07, 1.5115E-07, 1.3478E-07, 1.2236E-07, 1.1139E-07,    F04400
     *     1.0092E-07, 9.0795E-08, 8.2214E-08, 7.4691E-08, 6.7486E-08,    F04410
     *     6.0414E-08, 5.4584E-08, 4.8754E-08, 4.3501E-08, 3.8767E-08,    F04420
     *     3.4363E-08, 3.0703E-08, 2.7562E-08, 2.4831E-08, 2.2241E-08,    F04430
     *     1.9939E-08, 1.8049E-08, 1.6368E-08, 1.4863E-08, 1.3460E-08,    F04440
     *     1.2212E-08, 1.1155E-08, 1.0185E-08, 9.3417E-09, 8.5671E-09,    F04450
     *     7.8292E-09, 7.1749E-09, 6.5856E-09, 6.0588E-09, 5.5835E-09,    F04460
     *     5.1350E-09, 4.7395E-09, 4.3771E-09, 4.0476E-09, 3.7560E-09,    F04470
     *     3.4861E-09, 3.2427E-09, 3.0240E-09, 2.8278E-09, 2.6531E-09/    F04480
      DATA S0801/                                                         F04490
     *     2.4937E-09, 2.3511E-09, 2.2245E-09, 2.1133E-09, 2.0159E-09,    F04500
     *     1.9330E-09, 1.8669E-09, 1.8152E-09, 1.7852E-09, 1.7752E-09,    F04510
     *     1.7823E-09, 1.8194E-09, 1.8866E-09, 1.9759E-09, 2.0736E-09,    F04520
     *     2.2083E-09, 2.3587E-09, 2.4984E-09, 2.6333E-09, 2.8160E-09,    F04530
     *     3.0759E-09, 3.3720E-09, 3.6457E-09, 4.0668E-09, 4.4541E-09,    F04540
     *     4.7976E-09, 5.0908E-09, 5.4811E-09, 6.1394E-09, 6.3669E-09,    F04550
     *     6.5714E-09, 6.8384E-09, 7.1918E-09, 7.3741E-09, 7.2079E-09,    F04560
     *     7.2172E-09, 7.2572E-09, 7.3912E-09, 7.6188E-09, 8.3291E-09,    F04570
     *     8.7885E-09, 9.2412E-09, 1.0021E-08, 1.0752E-08, 1.1546E-08,    F04580
     *     1.1607E-08, 1.1949E-08, 1.2346E-08, 1.2516E-08, 1.2826E-08/    F04590
      DATA S0851/                                                         F04600
     *     1.3053E-08, 1.3556E-08, 1.4221E-08, 1.5201E-08, 1.6661E-08,    F04610
     *     1.8385E-08, 2.0585E-08, 2.3674E-08, 2.7928E-08, 3.3901E-08,    F04620
     *     4.1017E-08, 4.9595E-08, 6.0432E-08, 7.6304E-08, 9.0764E-08,    F04630
     *     1.0798E-07, 1.2442E-07, 1.4404E-07, 1.6331E-07, 1.8339E-07,    F04640
     *     2.0445E-07, 2.2288E-07, 2.3083E-07, 2.3196E-07, 2.3919E-07,    F04650
     *     2.3339E-07, 2.3502E-07, 2.3444E-07, 2.6395E-07, 2.9928E-07,    F04660
     *     3.0025E-07, 3.0496E-07, 3.1777E-07, 3.4198E-07, 3.4739E-07,    F04670
     *     3.2696E-07, 3.4100E-07, 3.5405E-07, 3.7774E-07, 3.8285E-07,    F04680
     *     3.6797E-07, 3.5800E-07, 3.2283E-07, 2.9361E-07, 2.4881E-07,    F04690
     *     2.0599E-07, 1.7121E-07, 1.3641E-07, 1.1111E-07, 8.9413E-08/    F04700
      DATA S0901/                                                         F04710
     *     7.3455E-08, 6.2078E-08, 5.2538E-08, 4.5325E-08, 3.9005E-08,    F04720
     *     3.4772E-08, 3.1203E-08, 2.8132E-08, 2.5250E-08, 2.2371E-08,    F04730
     *     2.0131E-08, 1.7992E-08, 1.6076E-08, 1.4222E-08, 1.2490E-08,    F04740
     *     1.1401E-08, 1.0249E-08, 9.2279E-09, 8.5654E-09, 7.6227E-09,    F04750
     *     6.9648E-09, 6.2466E-09, 5.7252E-09, 5.3800E-09, 4.6960E-09,    F04760
     *     4.2194E-09, 3.7746E-09, 3.3813E-09, 3.0656E-09, 2.6885E-09,    F04770
     *     2.4311E-09, 2.1572E-09, 1.8892E-09, 1.7038E-09, 1.4914E-09,    F04780
     *     1.3277E-09, 1.1694E-09, 1.0391E-09, 9.2779E-10, 8.3123E-10,    F04790
     *     7.4968E-10, 6.8385E-10, 6.2915E-10, 5.7784E-10, 5.2838E-10,    F04800
     *     4.8382E-10, 4.4543E-10, 4.1155E-10, 3.7158E-10, 3.3731E-10/    F04810
      DATA S0951/                                                         F04820
     *     3.0969E-10, 2.8535E-10, 2.6416E-10, 2.4583E-10, 2.2878E-10,    F04830
     *     2.1379E-10, 2.0073E-10, 1.8907E-10, 1.7866E-10, 1.6936E-10,    F04840
     *     1.6119E-10, 1.5424E-10, 1.4847E-10, 1.4401E-10, 1.4068E-10,    F04850
     *     1.3937E-10, 1.3943E-10, 1.4281E-10, 1.4766E-10, 1.5701E-10,    F04860
     *     1.7079E-10, 1.8691E-10, 2.0081E-10, 2.1740E-10, 2.4847E-10,    F04870
     *     2.6463E-10, 2.7087E-10, 2.7313E-10, 2.8352E-10, 2.9511E-10,    F04880
     *     2.8058E-10, 2.7227E-10, 2.7356E-10, 2.8012E-10, 2.8034E-10,    F04890
     *     2.9031E-10, 3.1030E-10, 3.3745E-10, 3.8152E-10, 4.0622E-10,    F04900
     *     4.2673E-10, 4.3879E-10, 4.5488E-10, 4.7179E-10, 4.6140E-10,    F04910
     *     4.6339E-10, 4.6716E-10, 4.7024E-10, 4.7931E-10, 4.8503E-10/    F04920
      DATA S1001/                                                         F04930
     *     4.9589E-10, 4.9499E-10, 5.0363E-10, 5.3184E-10, 5.6451E-10,    F04940
     *     6.0932E-10, 6.6469E-10, 7.4076E-10, 8.3605E-10, 9.4898E-10,    F04950
     *     1.0935E-09, 1.2593E-09, 1.4913E-09, 1.8099E-09, 2.1842E-09,    F04960
     *     2.7284E-09, 3.2159E-09, 3.7426E-09, 4.5226E-09, 5.3512E-09,    F04970
     *     6.1787E-09, 6.8237E-09, 7.9421E-09, 9.0002E-09, 9.6841E-09,    F04980
     *     9.9558E-09, 1.0232E-08, 1.0591E-08, 1.0657E-08, 1.0441E-08,    F04990
     *     1.0719E-08, 1.1526E-08, 1.2962E-08, 1.4336E-08, 1.6150E-08,    F05000
     *     1.8417E-08, 2.0725E-08, 2.3426E-08, 2.5619E-08, 2.7828E-08,    F05010
     *     3.0563E-08, 3.3438E-08, 3.6317E-08, 4.0400E-08, 4.4556E-08,    F05020
     *     5.0397E-08, 5.3315E-08, 5.9185E-08, 6.5311E-08, 6.9188E-08/    F05030
      DATA S1051/                                                         F05040
     *     7.7728E-08, 7.9789E-08, 8.6598E-08, 8.7768E-08, 9.1773E-08,    F05050
     *     9.7533E-08, 1.0007E-07, 1.0650E-07, 1.0992E-07, 1.0864E-07,    F05060
     *     1.0494E-07, 1.0303E-07, 1.0031E-07, 1.0436E-07, 1.0537E-07,    F05070
     *     1.1184E-07, 1.2364E-07, 1.3651E-07, 1.4881E-07, 1.4723E-07,    F05080
     *     1.4118E-07, 1.3371E-07, 1.1902E-07, 1.0007E-07, 7.9628E-08,    F05090
     *     6.4362E-08, 5.0243E-08, 3.8133E-08, 2.9400E-08, 2.3443E-08,    F05100
     *     1.9319E-08, 1.6196E-08, 1.4221E-08, 1.2817E-08, 1.1863E-08,    F05110
     *     1.1383E-08, 1.1221E-08, 1.1574E-08, 1.1661E-08, 1.2157E-08,    F05120
     *     1.2883E-08, 1.3295E-08, 1.4243E-08, 1.4240E-08, 1.4614E-08,    F05130
     *     1.4529E-08, 1.4685E-08, 1.4974E-08, 1.4790E-08, 1.4890E-08/    F05140
      DATA S1101/                                                         F05150
     *     1.4704E-08, 1.4142E-08, 1.3374E-08, 1.2746E-08, 1.2172E-08,    F05160
     *     1.2336E-08, 1.2546E-08, 1.3065E-08, 1.4090E-08, 1.5215E-08,    F05170
     *     1.6540E-08, 1.6144E-08, 1.5282E-08, 1.4358E-08, 1.2849E-08,    F05180
     *     1.0998E-08, 8.6956E-09, 7.0881E-09, 5.5767E-09, 4.2792E-09,    F05190
     *     3.2233E-09, 2.5020E-09, 1.9985E-09, 1.5834E-09, 1.3015E-09,    F05200
     *     1.0948E-09, 9.4141E-10, 8.1465E-10, 7.1517E-10, 6.2906E-10,    F05210
     *     5.5756E-10, 4.9805E-10, 4.3961E-10, 3.9181E-10, 3.5227E-10,    F05220
     *     3.1670E-10, 2.8667E-10, 2.5745E-10, 2.3212E-10, 2.0948E-10,    F05230
     *     1.8970E-10, 1.7239E-10, 1.5659E-10, 1.4301E-10, 1.3104E-10,    F05240
     *     1.2031E-10, 1.1095E-10, 1.0262E-10, 9.5130E-11, 8.8595E-11/    F05250
      DATA S1151/                                                         F05260
     *     8.2842E-11, 7.7727E-11, 7.3199E-11, 6.9286E-11, 6.5994E-11,    F05270
     *     6.3316E-11, 6.1244E-11, 5.9669E-11, 5.8843E-11, 5.8832E-11,    F05280
     *     5.9547E-11, 6.1635E-11, 6.4926E-11, 7.0745E-11, 7.8802E-11,    F05290
     *     8.6724E-11, 1.0052E-10, 1.1575E-10, 1.3626E-10, 1.5126E-10,    F05300
     *     1.6751E-10, 1.9239E-10, 2.1748E-10, 2.2654E-10, 2.2902E-10,    F05310
     *     2.3240E-10, 2.4081E-10, 2.3930E-10, 2.2378E-10, 2.2476E-10,    F05320
     *     2.2791E-10, 2.4047E-10, 2.5305E-10, 2.8073E-10, 3.1741E-10,    F05330
     *     3.6592E-10, 4.1495E-10, 4.6565E-10, 5.0990E-10, 5.5607E-10,    F05340
     *     6.1928E-10, 6.6779E-10, 7.3350E-10, 8.1434E-10, 8.9635E-10,    F05350
     *     9.9678E-10, 1.1256E-09, 1.2999E-09, 1.4888E-09, 1.7642E-09/    F05360
      DATA S1201/                                                         F05370
     *     1.9606E-09, 2.2066E-09, 2.4601E-09, 2.7218E-09, 3.0375E-09,    F05380
     *     3.1591E-09, 3.2852E-09, 3.2464E-09, 3.3046E-09, 3.2710E-09,    F05390
     *     3.2601E-09, 3.3398E-09, 3.7446E-09, 4.0795E-09, 4.0284E-09,    F05400
     *     4.0584E-09, 4.1677E-09, 4.5358E-09, 4.4097E-09, 4.2744E-09,    F05410
     *     4.5449E-09, 4.8147E-09, 5.2656E-09, 5.2476E-09, 5.0275E-09,    F05420
     *     4.7968E-09, 4.3654E-09, 3.9530E-09, 3.2447E-09, 2.6489E-09,    F05430
     *     2.1795E-09, 1.7880E-09, 1.4309E-09, 1.1256E-09, 9.1903E-10,    F05440
     *     7.6533E-10, 6.3989E-10, 5.5496E-10, 4.9581E-10, 4.5722E-10,    F05450
     *     4.3898E-10, 4.3505E-10, 4.3671E-10, 4.5329E-10, 4.6827E-10,    F05460
     *     4.9394E-10, 5.1122E-10, 5.1649E-10, 5.0965E-10, 4.9551E-10/    F05470
      DATA S1251/                                                         F05480
     *     4.8928E-10, 4.7947E-10, 4.7989E-10, 4.9071E-10, 4.8867E-10,    F05490
     *     4.7260E-10, 4.5756E-10, 4.5400E-10, 4.5993E-10, 4.4042E-10,    F05500
     *     4.3309E-10, 4.4182E-10, 4.6735E-10, 5.0378E-10, 5.2204E-10,    F05510
     *     5.0166E-10, 4.6799E-10, 4.3119E-10, 3.8803E-10, 3.3291E-10,    F05520
     *     2.6289E-10, 2.1029E-10, 1.7011E-10, 1.3345E-10, 1.0224E-10,    F05530
     *     7.8207E-11, 6.2451E-11, 5.0481E-11, 4.1507E-11, 3.5419E-11,    F05540
     *     3.0582E-11, 2.6900E-11, 2.3778E-11, 2.1343E-11, 1.9182E-11,    F05550
     *     1.7162E-11, 1.5391E-11, 1.3877E-11, 1.2619E-11, 1.1450E-11,    F05560
     *     1.0461E-11, 9.6578E-12, 8.9579E-12, 8.3463E-12, 7.8127E-12,    F05570
     *     7.3322E-12, 6.9414E-12, 6.6037E-12, 6.3285E-12, 6.1095E-12/    F05580
      DATA S1301/                                                         F05590
     *     5.9387E-12, 5.8118E-12, 5.7260E-12, 5.6794E-12, 5.6711E-12,    F05600
     *     5.7003E-12, 5.7670E-12, 5.8717E-12, 6.0151E-12, 6.1984E-12,    F05610
     *     6.4232E-12, 6.6918E-12, 7.0065E-12, 7.3705E-12, 7.7873E-12,    F05620
     *     8.2612E-12, 8.7972E-12, 9.4009E-12, 1.0079E-11, 1.0840E-11,    F05630
     *     1.1692E-11, 1.2648E-11, 1.3723E-11, 1.4935E-11, 1.6313E-11,    F05640
     *     1.7905E-11, 1.9740E-11, 2.1898E-11, 2.4419E-11, 2.7426E-11,    F05650
     *     3.0869E-11, 3.4235E-11, 3.7841E-11, 4.1929E-11, 4.6776E-11,    F05660
     *     5.2123E-11, 5.8497E-11, 6.5294E-11, 7.4038E-11, 8.4793E-11,    F05670
     *     9.6453E-11, 1.1223E-10, 1.2786E-10, 1.4882E-10, 1.7799E-10,    F05680
     *     2.0766E-10, 2.4523E-10, 2.8591E-10, 3.3386E-10, 4.0531E-10/    F05690
      DATA S1351/                                                         F05700
     *     4.7663E-10, 5.4858E-10, 6.3377E-10, 7.1688E-10, 8.4184E-10,    F05710
     *     9.5144E-10, 1.0481E-09, 1.1356E-09, 1.2339E-09, 1.3396E-09,    F05720
     *     1.4375E-09, 1.5831E-09, 1.7323E-09, 1.9671E-09, 2.2976E-09,    F05730
     *     2.6679E-09, 3.0777E-09, 3.4321E-09, 3.8192E-09, 4.2711E-09,    F05740
     *     4.4903E-09, 4.8931E-09, 5.2253E-09, 5.4040E-09, 5.6387E-09,    F05750
     *     5.6704E-09, 6.0345E-09, 6.1079E-09, 6.2576E-09, 6.4039E-09,    F05760
     *     6.3776E-09, 6.1878E-09, 5.8616E-09, 5.7036E-09, 5.5840E-09,    F05770
     *     5.6905E-09, 5.8931E-09, 6.2478E-09, 6.8291E-09, 7.4528E-09,    F05780
     *     7.6078E-09, 7.3898E-09, 6.7573E-09, 5.9827E-09, 5.0927E-09,    F05790
     *     4.0099E-09, 3.1933E-09, 2.4296E-09, 1.8485E-09, 1.4595E-09/    F05800
      DATA S1401/                                                         F05810
     *     1.2017E-09, 1.0164E-09, 8.7433E-10, 7.7108E-10, 7.0049E-10,    F05820
     *     6.5291E-10, 6.1477E-10, 5.9254E-10, 5.8150E-10, 5.7591E-10,    F05830
     *     5.8490E-10, 5.8587E-10, 5.9636E-10, 6.2408E-10, 6.5479E-10,    F05840
     *     7.0480E-10, 7.2313E-10, 7.5524E-10, 8.0863E-10, 8.3386E-10,    F05850
     *     9.2342E-10, 9.6754E-10, 1.0293E-09, 1.0895E-09, 1.1330E-09,    F05860
     *     1.2210E-09, 1.2413E-09, 1.2613E-09, 1.2671E-09, 1.2225E-09,    F05870
     *     1.1609E-09, 1.0991E-09, 1.0600E-09, 1.0570E-09, 1.0818E-09,    F05880
     *     1.1421E-09, 1.2270E-09, 1.3370E-09, 1.4742E-09, 1.4946E-09,    F05890
     *     1.4322E-09, 1.3210E-09, 1.1749E-09, 1.0051E-09, 7.8387E-10,    F05900
     *     6.1844E-10, 4.6288E-10, 3.4164E-10, 2.5412E-10, 1.9857E-10/    F05910
      DATA S1451/                                                         F05920
     *     1.5876E-10, 1.2966E-10, 1.0920E-10, 9.4811E-11, 8.3733E-11,    F05930
     *     7.3906E-11, 6.7259E-11, 6.1146E-11, 5.7119E-11, 5.3546E-11,    F05940
     *     4.8625E-11, 4.4749E-11, 4.1089E-11, 3.7825E-11, 3.4465E-11,    F05950
     *     3.1018E-11, 2.8109E-11, 2.5610E-11, 2.2859E-11, 2.0490E-11,    F05960
     *     1.8133E-11, 1.5835E-11, 1.3949E-11, 1.2295E-11, 1.0799E-11,    F05970
     *     9.6544E-12, 8.7597E-12, 7.9990E-12, 7.3973E-12, 6.9035E-12,    F05980
     *     6.4935E-12, 6.1195E-12, 5.8235E-12, 5.5928E-12, 5.4191E-12,    F05990
     *     5.2993E-12, 5.2338E-12, 5.2272E-12, 5.2923E-12, 5.4252E-12,    F06000
     *     5.6523E-12, 5.9433E-12, 6.3197E-12, 6.9016E-12, 7.5016E-12,    F06010
     *     8.2885E-12, 9.4050E-12, 1.0605E-11, 1.2257E-11, 1.3622E-11/    F06020
      DATA S1501/                                                         F06030
     *     1.5353E-11, 1.7543E-11, 1.9809E-11, 2.2197E-11, 2.4065E-11,    F06040
     *     2.6777E-11, 2.9751E-11, 3.2543E-11, 3.5536E-11, 3.9942E-11,    F06050
     *     4.6283E-11, 5.4556E-11, 6.5490E-11, 7.6803E-11, 9.0053E-11,    F06060
     *     1.0852E-10, 1.2946E-10, 1.4916E-10, 1.7748E-10, 2.0073E-10,    F06070
     *     2.2485E-10, 2.5114E-10, 2.7715E-10, 3.1319E-10, 3.3305E-10,    F06080
     *     3.5059E-10, 3.5746E-10, 3.6311E-10, 3.7344E-10, 3.6574E-10,    F06090
     *     3.7539E-10, 3.9434E-10, 4.3510E-10, 4.3340E-10, 4.2588E-10,    F06100
     *     4.3977E-10, 4.6062E-10, 4.7687E-10, 4.6457E-10, 4.8578E-10,    F06110
     *     5.2344E-10, 5.6752E-10, 5.8702E-10, 5.6603E-10, 5.3784E-10,    F06120
     *     4.9181E-10, 4.3272E-10, 3.5681E-10, 2.8814E-10, 2.3320E-10/    F06130
      DATA S1551/                                                         F06140
     *     1.8631E-10, 1.4587E-10, 1.1782E-10, 9.8132E-11, 8.2528E-11,    F06150
     *     6.9174E-11, 6.1056E-11, 5.3459E-11, 4.7116E-11, 4.1878E-11,    F06160
     *     3.8125E-11, 3.6347E-11, 3.5071E-11, 3.3897E-11, 3.3541E-11,    F06170
     *     3.3563E-11, 3.5469E-11, 3.8111E-11, 3.8675E-11, 4.1333E-11,    F06180
     *     4.3475E-11, 4.6476E-11, 4.9761E-11, 5.1380E-11, 5.4135E-11,    F06190
     *     5.3802E-11, 5.5158E-11, 5.6864E-11, 5.9311E-11, 6.3827E-11,    F06200
     *     6.7893E-11, 6.8230E-11, 6.6694E-11, 6.6018E-11, 6.4863E-11,    F06210
     *     6.5893E-11, 6.3813E-11, 6.4741E-11, 6.8630E-11, 7.0255E-11,    F06220
     *     7.0667E-11, 6.8810E-11, 6.4104E-11, 5.8136E-11, 4.7242E-11,    F06230
     *     3.7625E-11, 3.1742E-11, 2.5581E-11, 1.8824E-11, 1.3303E-11/    F06240
      DATA S1601/                                                         F06250
     *     9.6919E-12, 7.5353E-12, 6.0986E-12, 5.0742E-12, 4.3094E-12,    F06260
     *     3.7190E-12, 3.2520E-12, 2.8756E-12, 2.5680E-12, 2.3139E-12,    F06270
     *     2.1025E-12, 1.9257E-12, 1.7777E-12, 1.6539E-12, 1.5508E-12,    F06280
     *     1.4657E-12, 1.3966E-12, 1.3417E-12, 1.2998E-12, 1.2700E-12,    F06290
     *     1.2514E-12, 1.2437E-12, 1.2463E-12, 1.2592E-12, 1.2823E-12,    F06300
     *     1.3157E-12, 1.3596E-12, 1.4144E-12, 1.4806E-12, 1.5588E-12,    F06310
     *     1.6497E-12, 1.7544E-12, 1.8738E-12, 2.0094E-12, 2.1626E-12,    F06320
     *     2.3354E-12, 2.5297E-12, 2.7483E-12, 2.9941E-12, 3.2708E-12,    F06330
     *     3.5833E-12, 3.9374E-12, 4.3415E-12, 4.8079E-12, 5.3602E-12,    F06340
     *     5.9816E-12, 6.7436E-12, 7.6368E-12, 8.6812E-12, 9.8747E-12/    F06350
      DATA S1651/                                                         F06360
     *     1.1350E-11, 1.3181E-11, 1.5406E-11, 1.7868E-11, 2.0651E-11,    F06370
     *     2.4504E-11, 2.9184E-11, 3.4159E-11, 3.9979E-11, 4.8704E-11,    F06380
     *     5.7856E-11, 6.7576E-11, 7.9103E-11, 9.4370E-11, 1.1224E-10,    F06390
     *     1.3112E-10, 1.5674E-10, 1.8206E-10, 2.0576E-10, 2.3187E-10,    F06400
     *     2.7005E-10, 3.0055E-10, 3.3423E-10, 3.6956E-10, 3.8737E-10,    F06410
     *     4.2630E-10, 4.5154E-10, 4.8383E-10, 5.3582E-10, 5.8109E-10,    F06420
     *     6.3741E-10, 6.3874E-10, 6.3870E-10, 6.5818E-10, 6.5056E-10,    F06430
     *     6.5291E-10, 6.3159E-10, 6.3984E-10, 6.4549E-10, 6.5444E-10,    F06440
     *     6.7035E-10, 6.7665E-10, 6.9124E-10, 6.8451E-10, 6.9255E-10,    F06450
     *     6.9923E-10, 7.0396E-10, 6.7715E-10, 6.0371E-10, 5.3774E-10/    F06460
      DATA S1701/                                                         F06470
     *     4.6043E-10, 3.7635E-10, 2.9484E-10, 2.2968E-10, 1.8185E-10,    F06480
     *     1.4191E-10, 1.1471E-10, 9.4790E-11, 7.9613E-11, 6.7989E-11,    F06490
     *     5.9391E-11, 5.2810E-11, 4.7136E-11, 4.2618E-11, 3.8313E-11,    F06500
     *     3.4686E-11, 3.1669E-11, 2.9110E-11, 2.6871E-11, 2.5074E-11,    F06510
     *     2.4368E-11, 2.3925E-11, 2.4067E-11, 2.4336E-11, 2.4704E-11,    F06520
     *     2.5823E-11, 2.7177E-11, 2.9227E-11, 3.1593E-11, 3.5730E-11,    F06530
     *     4.0221E-11, 4.3994E-11, 4.8448E-11, 5.3191E-11, 5.8552E-11,    F06540
     *     6.3458E-11, 6.6335E-11, 7.2457E-11, 7.9091E-11, 8.2234E-11,    F06550
     *     8.7668E-11, 8.7951E-11, 9.2952E-11, 9.6157E-11, 9.5926E-11,    F06560
     *     1.0120E-10, 1.0115E-10, 9.9577E-11, 9.6633E-11, 9.2891E-11/    F06570
      DATA S1751/                                                         F06580
     *     9.3315E-11, 9.5584E-11, 1.0064E-10, 1.0509E-10, 1.1455E-10,    F06590
     *     1.2443E-10, 1.2963E-10, 1.2632E-10, 1.1308E-10, 1.0186E-10,    F06600
     *     8.5880E-11, 6.7863E-11, 5.1521E-11, 3.7780E-11, 2.8842E-11,    F06610
     *     2.2052E-11, 1.7402E-11, 1.4406E-11, 1.1934E-11, 1.0223E-11,    F06620
     *     8.9544E-12, 7.9088E-12, 7.0675E-12, 6.2222E-12, 5.6051E-12,    F06630
     *     5.0502E-12, 4.5578E-12, 4.2636E-12, 3.9461E-12, 3.7599E-12,    F06640
     *     3.5215E-12, 3.2467E-12, 3.0018E-12, 2.6558E-12, 2.3928E-12,    F06650
     *     2.0707E-12, 1.7575E-12, 1.5114E-12, 1.2941E-12, 1.1004E-12,    F06660
     *     9.5175E-13, 8.2894E-13, 7.3253E-13, 6.5551E-13, 5.9098E-13,    F06670
     *     5.3548E-13, 4.8697E-13, 4.4413E-13, 4.0600E-13, 3.7188E-13/    F06680
      DATA S1801/                                                         F06690
     *     3.4121E-13, 3.1356E-13, 2.8856E-13, 2.6590E-13, 2.4533E-13,    F06700
     *     2.2663E-13, 2.0960E-13, 1.9407E-13, 1.7990E-13, 1.6695E-13,    F06710
     *     1.5512E-13, 1.4429E-13, 1.3437E-13, 1.2527E-13, 1.1693E-13,    F06720
     *     1.0927E-13, 1.0224E-13, 9.5767E-14, 8.9816E-14, 8.4335E-14,    F06730
     *     7.9285E-14, 7.4626E-14, 7.0325E-14, 6.6352E-14, 6.2676E-14,    F06740
     *     5.9274E-14, 5.6121E-14, 5.3195E-14, 5.0479E-14, 4.7953E-14,    F06750
     *     4.5602E-14, 4.3411E-14, 4.1367E-14, 3.9456E-14, 3.7670E-14,    F06760
     *     3.5996E-14, 3.4427E-14, 3.2952E-14, 3.1566E-14, 3.0261E-14,    F06770
     *     2.9030E-14, 2.7868E-14, 2.6770E-14, 2.5730E-14, 2.4745E-14,    F06780
     *     2.3809E-14, 2.2921E-14, 2.2076E-14, 2.1271E-14, 2.0504E-14/    F06790
      DATA S1851/                                                         F06800
     *     1.9772E-14, 1.9073E-14, 1.8404E-14, 1.7764E-14, 1.7151E-14,    F06810
     *     1.6564E-14, 1.6000E-14, 1.5459E-14, 1.4939E-14, 1.4439E-14,    F06820
     *     1.3958E-14, 1.3495E-14, 1.3049E-14, 1.2620E-14, 1.2206E-14,    F06830
     *     1.1807E-14, 1.1422E-14, 1.1050E-14, 1.0691E-14, 1.0345E-14,    F06840
     *     1.0010E-14, 9.6870E-15, 9.3747E-15, 9.0727E-15, 8.7808E-15,    F06850
     *     8.4986E-15, 8.2257E-15, 7.9617E-15, 7.7064E-15, 7.4594E-15,    F06860
     *     7.2204E-15, 6.9891E-15, 6.7653E-15, 6.5488E-15, 6.3392E-15,    F06870
     *     6.1363E-15, 5.9399E-15, 5.7499E-15, 5.5659E-15, 5.3878E-15,    F06880
     *     5.2153E-15, 5.0484E-15, 4.8868E-15, 4.7303E-15, 4.5788E-15,    F06890
     *     4.4322E-15, 4.2902E-15, 4.1527E-15, 4.0196E-15, 3.8907E-15/    F06900
      DATA S1901/                                                         F06910
     *     3.7659E-15, 3.6451E-15, 3.5281E-15, 3.4149E-15, 3.3052E-15,    F06920
     *     3.1991E-15, 3.0963E-15, 2.9967E-15, 2.9004E-15, 2.8071E-15,    F06930
     *     2.7167E-15, 2.6293E-15, 2.5446E-15, 2.4626E-15, 2.3833E-15,    F06940
     *     2.3064E-15, 2.2320E-15, 2.1600E-15, 2.0903E-15, 2.0228E-15,    F06950
     *     1.9574E-15, 1.8942E-15, 1.8329E-15, 1.7736E-15, 1.7163E-15,    F06960
     *     1.6607E-15, 1.6069E-15, 1.5548E-15, 1.5044E-15, 1.4557E-15,    F06970
     *     1.4084E-15, 1.3627E-15, 1.3185E-15, 1.2757E-15, 1.2342E-15,    F06980
     *     1.1941E-15, 1.1552E-15, 1.1177E-15, 1.0813E-15, 1.0461E-15,    F06990
     *     1.0120E-15, 9.7900E-16, 9.4707E-16, 9.1618E-16, 8.8628E-16,    F07000
     *     8.5734E-16, 8.2933E-16, 8.0223E-16, 7.7600E-16, 7.5062E-16/    F07010
      DATA S1951/                                                         F07020
     *     7.2606E-16, 7.0229E-16, 6.7929E-16, 6.5703E-16, 6.3550E-16,    F07030
     *     6.1466E-16, 5.9449E-16, 5.7498E-16, 5.5610E-16, 5.3783E-16,    F07040
     *     5.2015E-16, 5.0305E-16, 4.8650E-16, 4.7049E-16, 4.5500E-16,    F07050
     *     4.4002E-16, 4.2552E-16, 4.1149E-16, 3.9792E-16, 3.8479E-16,    F07060
     *     3.7209E-16, 3.5981E-16, 3.4792E-16, 3.3642E-16, 3.2530E-16,    F07070
     *     3.1454E-16, 3.0413E-16, 2.9406E-16, 2.8432E-16, 2.7490E-16,    F07080
     *     2.6579E-16, 2.5697E-16, 2.4845E-16, 2.4020E-16, 2.3223E-16,    F07090
     *     2.2451E-16, 2.1705E-16, 2.0984E-16, 2.0286E-16, 1.9611E-16,    F07100
     *     1.8958E-16, 1.8327E-16, 1.7716E-16, 1.7126E-16, 1.6555E-16,    F07110
     *     1.6003E-16, 1.5469E-16, 1.4952E-16, 1.4453E-16, 1.3970E-16/    F07120
      DATA S2001/                                                         F07130
     *     1.3503E-16/                                                    F07140
C                                                                         F07150
      END                                                                 F07160
C
C     --------------------------------------------------------------
C
      SUBROUTINE SL260 (V1C,V2C,DVC,NPTC,C)                               F07170
C                                                                         F07180
      IMPLICIT REAL*8           (V)                                     ! F07190
C                                                                         F07200
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F07210
      COMMON /S260/ V1S,V2S,DVS,NPTS,S(2003)                              F07220
      DIMENSION C(*)                                                      F07230
C                                                                         F07240
      DVC = DVS                                                           F07250
      V1C = V1ABS-DVC                                                     F07260
      V2C = V2ABS+DVC                                                     F07270
C                                                                         F07280
      I1 = (V1C-V1S)/DVS                                                  F07290
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F07360
         I = I1+J                                                         F07370
         C(J) = 0.                                                        F07380
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F07390
         C(J) = S(I)                                                      F07400
   10 CONTINUE                                                            F07410
C                                                                         F07420
      RETURN                                                              F07430
C                                                                         F07440
      END                                                                 F07450
C
C     --------------------------------------------------------------
C
      BLOCK DATA BS260                                                    F07460
C                                                                         F07470
      IMPLICIT REAL*8           (V)                                     ! F07480
C                                                                         F07490
C               06/28/82                                                  F07500
C               UNITS OF (CM**3/MOL) * 1.E-20                             F07510
C                                                                         F07520
      COMMON /S260/ V1,V2,DV,NPT,                                         F07530
     *              S0000( 2),S0001(50),S0051(50),S0101(50),S0151(50),    F07540
     *              S0201(50),S0251(50),S0301(50),S0351(50),S0401(50),    F07550
     *              S0451(50),S0501(50),S0551(50),S0601(50),S0651(50),    F07560
     *              S0701(50),S0751(50),S0801(50),S0851(50),S0901(50),    F07570
     *              S0951(50),S1001(50),S1051(50),S1101(50),S1151(50),    F07580
     *              S1201(50),S1251(50),S1301(50),S1351(50),S1401(50),    F07590
     *              S1451(50),S1501(50),S1551(50),S1601(50),S1651(50),    F07600
     *              S1701(50),S1751(50),S1801(50),S1851(50),S1901(50),    F07610
     *              S1951(50),S2001( 1)                                   F07620
C                                                                         F07630
       DATA V1,V2,DV,NPT / -20.0, 20000.0, 10.0, 2003/                    F07640
C                                                                         F07650
C                                                                         F07660
      DATA S0000/                                                         F07670
     *     1.7750E-01, 1.7045E-01/                                        F07680
      DATA S0001/                                                         F07690
     *     1.6457E-01, 1.7045E-01, 1.7750E-01, 2.0036E-01, 2.1347E-01,    F07700
     *     2.2454E-01, 2.3428E-01, 2.3399E-01, 2.3022E-01, 2.0724E-01,    F07710
     *     1.9712E-01, 1.8317E-01, 1.6724E-01, 1.4780E-01, 1.2757E-01,    F07720
     *     1.1626E-01, 1.0098E-01, 8.9033E-02, 7.9770E-02, 6.7416E-02,    F07730
     *     5.9588E-02, 5.1117E-02, 4.6218E-02, 4.2179E-02, 3.4372E-02,    F07740
     *     2.9863E-02, 2.5252E-02, 2.2075E-02, 1.9209E-02, 1.5816E-02,    F07750
     *     1.3932E-02, 1.1943E-02, 1.0079E-02, 8.7667E-03, 7.4094E-03,    F07760
     *     6.4967E-03, 5.5711E-03, 4.8444E-03, 4.2552E-03, 3.6953E-03,    F07770
     *     3.2824E-03, 2.9124E-03, 2.6102E-03, 2.3370E-03, 2.1100E-03,    F07780
     *     1.9008E-03, 1.7145E-03, 1.5573E-03, 1.4206E-03, 1.2931E-03/    F07790
      DATA S0051/                                                         F07800
     *     1.1803E-03, 1.0774E-03, 9.8616E-04, 9.0496E-04, 8.3071E-04,    F07810
     *     7.6319E-04, 7.0149E-04, 6.4637E-04, 5.9566E-04, 5.4987E-04,    F07820
     *     5.0768E-04, 4.6880E-04, 4.3317E-04, 4.0037E-04, 3.7064E-04,    F07830
     *     3.4325E-04, 3.1809E-04, 2.9501E-04, 2.7382E-04, 2.5430E-04,    F07840
     *     2.3630E-04, 2.1977E-04, 2.0452E-04, 1.9042E-04, 1.7740E-04,    F07850
     *     1.6544E-04, 1.5442E-04, 1.4425E-04, 1.3486E-04, 1.2618E-04,    F07860
     *     1.1817E-04, 1.1076E-04, 1.0391E-04, 9.7563E-05, 9.1696E-05,    F07870
     *     8.6272E-05, 8.1253E-05, 7.6607E-05, 7.2302E-05, 6.8311E-05,    F07880
     *     6.4613E-05, 6.1183E-05, 5.8001E-05, 5.5048E-05, 5.2307E-05,    F07890
     *     4.9761E-05, 4.7395E-05, 4.5197E-05, 4.3155E-05, 4.1256E-05/    F07900
      DATA S0101/                                                         F07910
     *     3.9491E-05, 3.7849E-05, 3.6324E-05, 3.4908E-05, 3.3594E-05,    F07920
     *     3.2374E-05, 3.1244E-05, 3.0201E-05, 2.9240E-05, 2.8356E-05,    F07930
     *     2.7547E-05, 2.6814E-05, 2.6147E-05, 2.5551E-05, 2.5029E-05,    F07940
     *     2.4582E-05, 2.4203E-05, 2.3891E-05, 2.3663E-05, 2.3531E-05,    F07950
     *     2.3483E-05, 2.3516E-05, 2.3694E-05, 2.4032E-05, 2.4579E-05,    F07960
     *     2.5234E-05, 2.6032E-05, 2.7119E-05, 2.8631E-05, 3.0848E-05,    F07970
     *     3.3262E-05, 3.6635E-05, 4.0732E-05, 4.5923E-05, 5.3373E-05,    F07980
     *     6.1875E-05, 7.2031E-05, 8.5980E-05, 9.8642E-05, 1.1469E-04,    F07990
     *     1.3327E-04, 1.5390E-04, 1.7513E-04, 2.0665E-04, 2.3609E-04,    F08000
     *     2.6220E-04, 2.8677E-04, 3.2590E-04, 3.8624E-04, 4.1570E-04/    F08010
      DATA S0151/                                                         F08020
     *     4.5207E-04, 4.9336E-04, 5.4500E-04, 5.8258E-04, 5.8086E-04,    F08030
     *     5.6977E-04, 5.3085E-04, 4.8020E-04, 4.3915E-04, 4.0343E-04,    F08040
     *     3.7853E-04, 3.7025E-04, 3.9637E-04, 4.4675E-04, 4.7072E-04,    F08050
     *     4.9022E-04, 5.2076E-04, 5.3676E-04, 5.2755E-04, 4.8244E-04,    F08060
     *     4.5473E-04, 4.3952E-04, 3.9614E-04, 3.4086E-04, 2.9733E-04,    F08070
     *     2.6367E-04, 2.3767E-04, 2.0427E-04, 1.7595E-04, 1.5493E-04,    F08080
     *     1.3851E-04, 1.1874E-04, 1.0735E-04, 9.0490E-05, 8.1149E-05,    F08090
     *     7.4788E-05, 6.5438E-05, 5.8248E-05, 4.8076E-05, 4.3488E-05,    F08100
     *     3.7856E-05, 3.3034E-05, 2.9592E-05, 2.6088E-05, 2.3497E-05,    F08110
     *     2.0279E-05, 1.7526E-05, 1.5714E-05, 1.3553E-05, 1.2145E-05/    F08120
      DATA S0201/                                                         F08130
     *     1.0802E-05, 9.7681E-06, 8.8196E-06, 7.8291E-06, 7.1335E-06,    F08140
     *     6.4234E-06, 5.8391E-06, 5.3532E-06, 4.9079E-06, 4.5378E-06,    F08150
     *     4.1716E-06, 3.8649E-06, 3.5893E-06, 3.3406E-06, 3.1199E-06,    F08160
     *     2.9172E-06, 2.7348E-06, 2.5644E-06, 2.4086E-06, 2.2664E-06,    F08170
     *     2.1359E-06, 2.0159E-06, 1.9051E-06, 1.8031E-06, 1.7074E-06,    F08180
     *     1.6185E-06, 1.5356E-06, 1.4584E-06, 1.3861E-06, 1.3179E-06,    F08190
     *     1.2545E-06, 1.1951E-06, 1.1395E-06, 1.0873E-06, 1.0384E-06,    F08200
     *     9.9250E-07, 9.4935E-07, 9.0873E-07, 8.7050E-07, 8.3446E-07,    F08210
     *     8.0046E-07, 7.6834E-07, 7.3800E-07, 7.0931E-07, 6.8217E-07,    F08220
     *     6.5648E-07, 6.3214E-07, 6.0909E-07, 5.8725E-07, 5.6655E-07/    F08230
      DATA S0251/                                                         F08240
     *     5.4693E-07, 5.2835E-07, 5.1077E-07, 4.9416E-07, 4.7853E-07,    F08250
     *     4.6381E-07, 4.5007E-07, 4.3728E-07, 4.2550E-07, 4.1450E-07,    F08260
     *     4.0459E-07, 3.9532E-07, 3.8662E-07, 3.7855E-07, 3.7041E-07,    F08270
     *     3.6254E-07, 3.5420E-07, 3.4617E-07, 3.3838E-07, 3.3212E-07,    F08280
     *     3.2655E-07, 3.1865E-07, 3.1203E-07, 3.0670E-07, 3.0252E-07,    F08290
     *     2.9749E-07, 2.9184E-07, 2.8795E-07, 2.8501E-07, 2.8202E-07,    F08300
     *     2.7856E-07, 2.7509E-07, 2.7152E-07, 2.6844E-07, 2.6642E-07,    F08310
     *     2.6548E-07, 2.6617E-07, 2.6916E-07, 2.7372E-07, 2.8094E-07,    F08320
     *     2.9236E-07, 3.1035E-07, 3.2854E-07, 3.5481E-07, 3.9377E-07,    F08330
     *     4.4692E-07, 5.0761E-07, 5.7715E-07, 6.7725E-07, 8.0668E-07/    F08340
      DATA S0301/                                                         F08350
     *     9.3716E-07, 1.0797E-06, 1.1689E-06, 1.3217E-06, 1.4814E-06,    F08360
     *     1.5627E-06, 1.6519E-06, 1.7601E-06, 1.9060E-06, 2.0474E-06,    F08370
     *     2.0716E-06, 2.0433E-06, 1.9752E-06, 1.8466E-06, 1.7526E-06,    F08380
     *     1.6657E-06, 1.5870E-06, 1.5633E-06, 1.6520E-06, 1.8471E-06,    F08390
     *     1.9953E-06, 2.0975E-06, 2.2016E-06, 2.2542E-06, 2.3081E-06,    F08400
     *     2.3209E-06, 2.2998E-06, 2.3056E-06, 2.2757E-06, 2.2685E-06,    F08410
     *     2.2779E-06, 2.2348E-06, 2.2445E-06, 2.3174E-06, 2.4284E-06,    F08420
     *     2.5290E-06, 2.7340E-06, 2.9720E-06, 3.2332E-06, 3.5392E-06,    F08430
     *     3.9013E-06, 4.3334E-06, 4.9088E-06, 5.3428E-06, 5.9142E-06,    F08440
     *     6.6106E-06, 7.4709E-06, 8.5019E-06, 9.6835E-06, 1.0984E-05/    F08450
      DATA S0351/                                                         F08460
     *     1.2831E-05, 1.4664E-05, 1.7080E-05, 2.0103E-05, 2.4148E-05,    F08470
     *     2.7948E-05, 3.2855E-05, 3.9046E-05, 4.6429E-05, 5.6633E-05,    F08480
     *     6.6305E-05, 7.6048E-05, 8.7398E-05, 1.0034E-04, 1.1169E-04,    F08490
     *     1.2813E-04, 1.3354E-04, 1.3952E-04, 1.4204E-04, 1.4615E-04,    F08500
     *     1.5144E-04, 1.5475E-04, 1.6561E-04, 1.7135E-04, 1.6831E-04,    F08510
     *     1.6429E-04, 1.6353E-04, 1.6543E-04, 1.5944E-04, 1.5404E-04,    F08520
     *     1.5458E-04, 1.6287E-04, 1.7277E-04, 1.8387E-04, 1.7622E-04,    F08530
     *     1.6360E-04, 1.5273E-04, 1.3667E-04, 1.2364E-04, 9.7576E-05,    F08540
     *     7.9140E-05, 6.4241E-05, 5.1826E-05, 4.1415E-05, 3.1347E-05,    F08550
     *     2.5125E-05, 2.0027E-05, 1.6362E-05, 1.3364E-05, 1.1117E-05/    F08560
      DATA S0401/                                                         F08570
     *     9.4992E-06, 8.1581E-06, 7.1512E-06, 6.2692E-06, 5.5285E-06,    F08580
     *     4.9000E-06, 4.3447E-06, 3.8906E-06, 3.4679E-06, 3.1089E-06,    F08590
     *     2.8115E-06, 2.5496E-06, 2.2982E-06, 2.0861E-06, 1.8763E-06,    F08600
     *     1.7035E-06, 1.5548E-06, 1.4107E-06, 1.2839E-06, 1.1706E-06,    F08610
     *     1.0709E-06, 9.8099E-07, 8.9901E-07, 8.2394E-07, 7.5567E-07,    F08620
     *     6.9434E-07, 6.3867E-07, 5.8845E-07, 5.4263E-07, 5.0033E-07,    F08630
     *     4.6181E-07, 4.2652E-07, 3.9437E-07, 3.6497E-07, 3.3781E-07,    F08640
     *     3.1292E-07, 2.9011E-07, 2.6915E-07, 2.4989E-07, 2.3215E-07,    F08650
     *     2.1582E-07, 2.0081E-07, 1.8700E-07, 1.7432E-07, 1.6264E-07,    F08660
     *     1.5191E-07, 1.4207E-07, 1.3306E-07, 1.2484E-07, 1.1737E-07/    F08670
      DATA S0451/                                                         F08680
     *     1.1056E-07, 1.0451E-07, 9.9060E-08, 9.4135E-08, 8.9608E-08,    F08690
     *     8.5697E-08, 8.1945E-08, 7.8308E-08, 7.4808E-08, 7.1686E-08,    F08700
     *     6.8923E-08, 6.5869E-08, 6.3308E-08, 6.0840E-08, 5.8676E-08,    F08710
     *     5.6744E-08, 5.5016E-08, 5.3813E-08, 5.2792E-08, 5.2097E-08,    F08720
     *     5.1737E-08, 5.1603E-08, 5.1656E-08, 5.1989E-08, 5.2467E-08,    F08730
     *     5.2918E-08, 5.3589E-08, 5.4560E-08, 5.5869E-08, 5.7403E-08,    F08740
     *     5.8968E-08, 6.0973E-08, 6.3432E-08, 6.6245E-08, 6.9353E-08,    F08750
     *     7.2686E-08, 7.6541E-08, 8.0991E-08, 8.5950E-08, 9.1429E-08,    F08760
     *     9.7851E-08, 1.0516E-07, 1.1349E-07, 1.2295E-07, 1.3335E-07,    F08770
     *     1.4488E-07, 1.5864E-07, 1.7412E-07, 1.9140E-07, 2.1078E-07/    F08780
      DATA S0501/                                                         F08790
     *     2.3369E-07, 2.5996E-07, 2.8848E-07, 3.2169E-07, 3.5991E-07,    F08800
     *     4.0566E-07, 4.5969E-07, 5.3094E-07, 6.1458E-07, 7.1155E-07,    F08810
     *     8.3045E-07, 9.9021E-07, 1.2042E-06, 1.4914E-06, 1.8145E-06,    F08820
     *     2.2210E-06, 2.7831E-06, 3.4533E-06, 4.4446E-06, 5.1989E-06,    F08830
     *     6.2289E-06, 7.1167E-06, 8.3949E-06, 9.6417E-06, 1.0313E-05,    F08840
     *     1.0485E-05, 1.0641E-05, 1.0898E-05, 1.0763E-05, 1.0506E-05,    F08850
     *     1.0497E-05, 1.1696E-05, 1.2654E-05, 1.3029E-05, 1.3175E-05,    F08860
     *     1.4264E-05, 1.4985E-05, 1.4999E-05, 1.4317E-05, 1.4616E-05,    F08870
     *     1.4963E-05, 1.5208E-05, 1.4942E-05, 1.3879E-05, 1.3087E-05,    F08880
     *     1.1727E-05, 1.0515E-05, 9.0073E-06, 7.3133E-06, 6.1181E-06/    F08890
      DATA S0551/                                                         F08900
     *     5.0623E-06, 4.1105E-06, 3.3915E-06, 2.6711E-06, 2.1464E-06,    F08910
     *     1.7335E-06, 1.4302E-06, 1.1847E-06, 9.9434E-07, 8.2689E-07,    F08920
     *     7.0589E-07, 6.0750E-07, 5.3176E-07, 4.6936E-07, 4.1541E-07,    F08930
     *     3.6625E-07, 3.2509E-07, 2.9156E-07, 2.6308E-07, 2.3819E-07,    F08940
     *     2.1421E-07, 1.9366E-07, 1.7626E-07, 1.5982E-07, 1.4567E-07,    F08950
     *     1.3354E-07, 1.2097E-07, 1.1029E-07, 1.0063E-07, 9.2003E-08,    F08960
     *     8.4245E-08, 7.7004E-08, 7.0636E-08, 6.4923E-08, 5.9503E-08,    F08970
     *     5.4742E-08, 5.0450E-08, 4.6470E-08, 4.2881E-08, 3.9550E-08,    F08980
     *     3.6541E-08, 3.3803E-08, 3.1279E-08, 2.8955E-08, 2.6858E-08,    F08990
     *     2.4905E-08, 2.3146E-08, 2.1539E-08, 2.0079E-08, 1.8746E-08/    F09000
      DATA S0601/                                                         F09010
     *     1.7517E-08, 1.6396E-08, 1.5369E-08, 1.4426E-08, 1.3543E-08,    F09020
     *     1.2724E-08, 1.1965E-08, 1.1267E-08, 1.0617E-08, 1.0010E-08,    F09030
     *     9.4662E-09, 8.9553E-09, 8.4988E-09, 8.0807E-09, 7.7043E-09,    F09040
     *     7.3721E-09, 7.0707E-09, 6.8047E-09, 6.5702E-09, 6.3634E-09,    F09050
     *     6.1817E-09, 6.0239E-09, 5.8922E-09, 5.7824E-09, 5.7019E-09,    F09060
     *     5.6368E-09, 5.5940E-09, 5.5669E-09, 5.5583E-09, 5.5653E-09,    F09070
     *     5.5837E-09, 5.6243E-09, 5.6883E-09, 5.7800E-09, 5.8964E-09,    F09080
     *     6.0429E-09, 6.2211E-09, 6.4282E-09, 6.6634E-09, 6.9306E-09,    F09090
     *     7.2336E-09, 7.5739E-09, 7.9562E-09, 8.3779E-09, 8.8575E-09,    F09100
     *     9.3992E-09, 1.0004E-08, 1.0684E-08, 1.1450E-08, 1.2320E-08/    F09110
      DATA S0651/                                                         F09120
     *     1.3311E-08, 1.4455E-08, 1.5758E-08, 1.7254E-08, 1.8927E-08,    F09130
     *     2.0930E-08, 2.3348E-08, 2.6074E-08, 2.9221E-08, 3.2770E-08,    F09140
     *     3.7485E-08, 4.2569E-08, 4.8981E-08, 5.5606E-08, 6.2393E-08,    F09150
     *     7.1901E-08, 8.2921E-08, 9.5513E-08, 1.1111E-07, 1.3143E-07,    F09160
     *     1.5971E-07, 1.8927E-07, 2.2643E-07, 2.7860E-07, 3.2591E-07,    F09170
     *     3.7024E-07, 4.2059E-07, 4.9432E-07, 5.5543E-07, 5.7498E-07,    F09180
     *     5.9210E-07, 6.1005E-07, 6.1577E-07, 5.9193E-07, 5.6602E-07,    F09190
     *     5.7403E-07, 6.0050E-07, 6.4723E-07, 6.7073E-07, 7.5415E-07,    F09200
     *     8.0982E-07, 8.7658E-07, 9.1430E-07, 9.4459E-07, 9.8347E-07,    F09210
     *     9.8768E-07, 1.0153E-06, 1.0066E-06, 1.0353E-06, 1.0353E-06/    F09220
      DATA S0701/                                                         F09230
     *     1.0722E-06, 1.1138E-06, 1.1923E-06, 1.2947E-06, 1.4431E-06,    F09240
     *     1.6537E-06, 1.8662E-06, 2.2473E-06, 2.6464E-06, 3.1041E-06,    F09250
     *     3.4858E-06, 4.0167E-06, 4.6675E-06, 5.0983E-06, 5.7997E-06,    F09260
     *     6.0503E-06, 6.4687E-06, 6.5396E-06, 6.7986E-06, 7.0244E-06,    F09270
     *     7.2305E-06, 7.6732E-06, 7.9783E-06, 7.9846E-06, 7.7617E-06,    F09280
     *     7.7657E-06, 7.7411E-06, 7.8816E-06, 7.8136E-06, 8.0051E-06,    F09290
     *     8.5799E-06, 9.1659E-06, 9.8646E-06, 9.4920E-06, 8.7670E-06,    F09300
     *     8.2034E-06, 7.2297E-06, 6.2324E-06, 4.9315E-06, 3.9128E-06,    F09310
     *     3.1517E-06, 2.4469E-06, 1.8815E-06, 1.4627E-06, 1.1698E-06,    F09320
     *     9.4686E-07, 7.8486E-07, 6.6970E-07, 5.8811E-07, 5.2198E-07/    F09330
      DATA S0751/                                                         F09340
     *     4.6809E-07, 4.1671E-07, 3.7006E-07, 3.3066E-07, 2.9387E-07,    F09350
     *     2.6415E-07, 2.3409E-07, 2.0991E-07, 1.9132E-07, 1.7519E-07,    F09360
     *     1.5939E-07, 1.4368E-07, 1.3050E-07, 1.1883E-07, 1.0772E-07,    F09370
     *     9.6884E-08, 8.7888E-08, 7.8956E-08, 7.1024E-08, 6.3824E-08,    F09380
     *     5.7256E-08, 5.1769E-08, 4.7037E-08, 4.2901E-08, 3.8970E-08,    F09390
     *     3.5467E-08, 3.2502E-08, 2.9827E-08, 2.7389E-08, 2.5111E-08,    F09400
     *     2.3056E-08, 2.1267E-08, 1.9610E-08, 1.8133E-08, 1.6775E-08,    F09410
     *     1.5491E-08, 1.4329E-08, 1.3265E-08, 1.2300E-08, 1.1420E-08,    F09420
     *     1.0593E-08, 9.8475E-09, 9.1585E-09, 8.5256E-09, 7.9525E-09,    F09430
     *     7.4226E-09, 6.9379E-09, 6.4950E-09, 6.0911E-09, 5.7242E-09/    F09440
      DATA S0801/                                                         F09450
     *     5.3877E-09, 5.0821E-09, 4.8051E-09, 4.5554E-09, 4.3315E-09,    F09460
     *     4.1336E-09, 3.9632E-09, 3.8185E-09, 3.7080E-09, 3.6296E-09,    F09470
     *     3.5804E-09, 3.5776E-09, 3.6253E-09, 3.7115E-09, 3.8151E-09,    F09480
     *     3.9804E-09, 4.1742E-09, 4.3581E-09, 4.5306E-09, 4.7736E-09,    F09490
     *     5.1297E-09, 5.5291E-09, 5.9125E-09, 6.4956E-09, 7.0362E-09,    F09500
     *     7.5318E-09, 7.9947E-09, 8.6438E-09, 9.7227E-09, 1.0130E-08,    F09510
     *     1.0549E-08, 1.1064E-08, 1.1702E-08, 1.2043E-08, 1.1781E-08,    F09520
     *     1.1838E-08, 1.1917E-08, 1.2131E-08, 1.2476E-08, 1.3611E-08,    F09530
     *     1.4360E-08, 1.5057E-08, 1.6247E-08, 1.7284E-08, 1.8420E-08,    F09540
     *     1.8352E-08, 1.8722E-08, 1.9112E-08, 1.9092E-08, 1.9311E-08/    F09550
      DATA S0851/                                                         F09560
     *     1.9411E-08, 1.9884E-08, 2.0508E-08, 2.1510E-08, 2.3143E-08,    F09570
     *     2.5050E-08, 2.7596E-08, 3.1231E-08, 3.6260E-08, 4.3410E-08,    F09580
     *     5.2240E-08, 6.3236E-08, 7.7522E-08, 9.8688E-08, 1.1859E-07,    F09590
     *     1.4341E-07, 1.6798E-07, 1.9825E-07, 2.2898E-07, 2.6257E-07,    F09600
     *     2.9884E-07, 3.3247E-07, 3.4936E-07, 3.5583E-07, 3.7150E-07,    F09610
     *     3.6580E-07, 3.7124E-07, 3.7030E-07, 4.1536E-07, 4.6656E-07,    F09620
     *     4.6677E-07, 4.7507E-07, 4.9653E-07, 5.3795E-07, 5.4957E-07,    F09630
     *     5.2238E-07, 5.4690E-07, 5.6569E-07, 5.9844E-07, 5.9835E-07,    F09640
     *     5.6522E-07, 5.4123E-07, 4.7904E-07, 4.2851E-07, 3.5603E-07,    F09650
     *     2.8932E-07, 2.3655E-07, 1.8592E-07, 1.4943E-07, 1.1971E-07/    F09660
      DATA S0901/                                                         F09670
     *     9.8482E-08, 8.3675E-08, 7.1270E-08, 6.2496E-08, 5.4999E-08,    F09680
     *     4.9821E-08, 4.5387E-08, 4.1340E-08, 3.7453E-08, 3.3298E-08,    F09690
     *     3.0120E-08, 2.7032E-08, 2.4236E-08, 2.1500E-08, 1.8988E-08,    F09700
     *     1.7414E-08, 1.5706E-08, 1.4192E-08, 1.3204E-08, 1.1759E-08,    F09710
     *     1.0737E-08, 9.6309E-09, 8.8179E-09, 8.2619E-09, 7.2264E-09,    F09720
     *     6.4856E-09, 5.8037E-09, 5.2093E-09, 4.7205E-09, 4.1749E-09,    F09730
     *     3.7852E-09, 3.3915E-09, 3.0089E-09, 2.7335E-09, 2.4398E-09,    F09740
     *     2.2031E-09, 1.9786E-09, 1.7890E-09, 1.6266E-09, 1.4830E-09,    F09750
     *     1.3576E-09, 1.2518E-09, 1.1587E-09, 1.0726E-09, 9.9106E-10,    F09760
     *     9.1673E-10, 8.5084E-10, 7.9147E-10, 7.2882E-10, 6.7342E-10/    F09770
      DATA S0951/                                                         F09780
     *     6.2593E-10, 5.8294E-10, 5.4435E-10, 5.0997E-10, 4.7806E-10,    F09790
     *     4.4931E-10, 4.2357E-10, 4.0023E-10, 3.7909E-10, 3.5999E-10,    F09800
     *     3.4285E-10, 3.2776E-10, 3.1468E-10, 3.0377E-10, 2.9479E-10,    F09810
     *     2.8877E-10, 2.8512E-10, 2.8617E-10, 2.8976E-10, 3.0001E-10,    F09820
     *     3.1718E-10, 3.3898E-10, 3.5857E-10, 3.8358E-10, 4.3131E-10,    F09830
     *     4.5741E-10, 4.6948E-10, 4.7594E-10, 4.9529E-10, 5.1563E-10,    F09840
     *     4.9475E-10, 4.8369E-10, 4.8829E-10, 5.0047E-10, 5.0203E-10,    F09850
     *     5.1954E-10, 5.5352E-10, 5.9928E-10, 6.7148E-10, 7.1121E-10,    F09860
     *     7.4317E-10, 7.6039E-10, 7.8313E-10, 8.0684E-10, 7.8553E-10,    F09870
     *     7.8312E-10, 7.8537E-10, 7.8872E-10, 8.0185E-10, 8.1004E-10/    F09880
      DATA S1001/                                                         F09890
     *     8.2608E-10, 8.2525E-10, 8.3857E-10, 8.7920E-10, 9.2451E-10,    F09900
     *     9.8661E-10, 1.0629E-09, 1.1659E-09, 1.2922E-09, 1.4387E-09,    F09910
     *     1.6254E-09, 1.8425E-09, 2.1428E-09, 2.5477E-09, 3.0379E-09,    F09920
     *     3.7570E-09, 4.4354E-09, 5.1802E-09, 6.2769E-09, 7.4894E-09,    F09930
     *     8.7474E-09, 9.8037E-09, 1.1582E-08, 1.3293E-08, 1.4471E-08,    F09940
     *     1.5025E-08, 1.5580E-08, 1.6228E-08, 1.6413E-08, 1.6020E-08,    F09950
     *     1.6393E-08, 1.7545E-08, 1.9590E-08, 2.1449E-08, 2.3856E-08,    F09960
     *     2.7050E-08, 3.0214E-08, 3.3733E-08, 3.6487E-08, 3.9353E-08,    F09970
     *     4.2660E-08, 4.6385E-08, 4.9955E-08, 5.5313E-08, 6.0923E-08,    F09980
     *     6.8948E-08, 7.3649E-08, 8.2602E-08, 9.2212E-08, 9.9080E-08/    F09990
      DATA S1051/                                                         F10000
     *     1.1319E-07, 1.1790E-07, 1.2941E-07, 1.3199E-07, 1.3914E-07,    F10010
     *     1.4843E-07, 1.5300E-07, 1.6419E-07, 1.7095E-07, 1.6988E-07,    F10020
     *     1.6494E-07, 1.6327E-07, 1.6067E-07, 1.6909E-07, 1.7118E-07,    F10030
     *     1.8106E-07, 1.9857E-07, 2.1696E-07, 2.3385E-07, 2.2776E-07,    F10040
     *     2.1402E-07, 1.9882E-07, 1.7362E-07, 1.4308E-07, 1.1158E-07,    F10050
     *     8.8781E-08, 6.8689E-08, 5.2062E-08, 4.0427E-08, 3.2669E-08,    F10060
     *     2.7354E-08, 2.3200E-08, 2.0580E-08, 1.8676E-08, 1.7329E-08,    F10070
     *     1.6621E-08, 1.6433E-08, 1.6953E-08, 1.7134E-08, 1.7948E-08,    F10080
     *     1.9107E-08, 1.9875E-08, 2.1416E-08, 2.1556E-08, 2.2265E-08,    F10090
     *     2.2171E-08, 2.2534E-08, 2.3029E-08, 2.2828E-08, 2.3143E-08/    F10100
      DATA S1101/                                                         F10110
     *     2.2965E-08, 2.2223E-08, 2.1108E-08, 2.0265E-08, 1.9516E-08,    F10120
     *     1.9941E-08, 2.0312E-08, 2.1080E-08, 2.2611E-08, 2.4210E-08,    F10130
     *     2.6069E-08, 2.5097E-08, 2.3318E-08, 2.1543E-08, 1.8942E-08,    F10140
     *     1.5960E-08, 1.2386E-08, 9.9340E-09, 7.7502E-09, 5.9462E-09,    F10150
     *     4.5113E-09, 3.5523E-09, 2.8844E-09, 2.3394E-09, 1.9584E-09,    F10160
     *     1.6749E-09, 1.4624E-09, 1.2809E-09, 1.1359E-09, 1.0087E-09,    F10170
     *     9.0166E-10, 8.1079E-10, 7.2219E-10, 6.4922E-10, 5.8803E-10,    F10180
     *     5.3290E-10, 4.8590E-10, 4.4111E-10, 4.0184E-10, 3.6644E-10,    F10190
     *     3.3529E-10, 3.0789E-10, 2.8286E-10, 2.6089E-10, 2.4125E-10,    F10200
     *     2.2355E-10, 2.0783E-10, 1.9370E-10, 1.8088E-10, 1.6948E-10/    F10210
      DATA S1151/                                                         F10220
     *     1.5929E-10, 1.5013E-10, 1.4193E-10, 1.3470E-10, 1.2841E-10,    F10230
     *     1.2307E-10, 1.1865E-10, 1.1502E-10, 1.1243E-10, 1.1099E-10,    F10240
     *     1.1066E-10, 1.1216E-10, 1.1529E-10, 1.2171E-10, 1.3128E-10,    F10250
     *     1.4153E-10, 1.5962E-10, 1.8048E-10, 2.0936E-10, 2.3165E-10,    F10260
     *     2.5746E-10, 2.9600E-10, 3.3707E-10, 3.5267E-10, 3.5953E-10,    F10270
     *     3.6822E-10, 3.8363E-10, 3.8286E-10, 3.5883E-10, 3.6154E-10,    F10280
     *     3.6653E-10, 3.8507E-10, 4.0250E-10, 4.4435E-10, 4.9889E-10,    F10290
     *     5.6932E-10, 6.3599E-10, 7.0281E-10, 7.5777E-10, 8.1279E-10,    F10300
     *     8.8910E-10, 9.3400E-10, 1.0076E-09, 1.0945E-09, 1.1898E-09,    F10310
     *     1.3108E-09, 1.4725E-09, 1.7028E-09, 1.9619E-09, 2.3527E-09/    F10320
      DATA S1201/                                                         F10330
     *     2.6488E-09, 3.0327E-09, 3.4396E-09, 3.8797E-09, 4.4115E-09,    F10340
     *     4.6853E-09, 4.9553E-09, 4.9551E-09, 5.1062E-09, 5.0996E-09,    F10350
     *     5.1119E-09, 5.2283E-09, 5.8297E-09, 6.3439E-09, 6.2675E-09,    F10360
     *     6.3296E-09, 6.5173E-09, 7.1685E-09, 7.0528E-09, 6.8856E-09,    F10370
     *     7.3182E-09, 7.6990E-09, 8.3461E-09, 8.1946E-09, 7.7153E-09,    F10380
     *     7.2411E-09, 6.4511E-09, 5.7336E-09, 4.6105E-09, 3.6962E-09,    F10390
     *     2.9944E-09, 2.4317E-09, 1.9399E-09, 1.5331E-09, 1.2633E-09,    F10400
     *     1.0613E-09, 9.0136E-10, 7.9313E-10, 7.1543E-10, 6.6485E-10,    F10410
     *     6.4225E-10, 6.3980E-10, 6.4598E-10, 6.7428E-10, 7.0270E-10,    F10420
     *     7.4694E-10, 7.7946E-10, 7.9395E-10, 7.8716E-10, 7.6933E-10/    F10430
      DATA S1251/                                                         F10440
     *     7.6220E-10, 7.4825E-10, 7.4805E-10, 7.6511E-10, 7.6492E-10,    F10450
     *     7.4103E-10, 7.1979E-10, 7.1686E-10, 7.3403E-10, 7.1142E-10,    F10460
     *     7.0212E-10, 7.1548E-10, 7.5253E-10, 8.0444E-10, 8.2378E-10,    F10470
     *     7.8004E-10, 7.1712E-10, 6.4978E-10, 5.7573E-10, 4.8675E-10,    F10480
     *     3.7945E-10, 3.0118E-10, 2.4241E-10, 1.9100E-10, 1.4816E-10,    F10490
     *     1.1567E-10, 9.4183E-11, 7.7660E-11, 6.5270E-11, 5.6616E-11,    F10500
     *     4.9576E-11, 4.4137E-11, 3.9459E-11, 3.5759E-11, 3.2478E-11,    F10510
     *     2.9419E-11, 2.6703E-11, 2.4365E-11, 2.2412E-11, 2.0606E-11,    F10520
     *     1.9067E-11, 1.7800E-11, 1.6695E-11, 1.5729E-11, 1.4887E-11,    F10530
     *     1.4135E-11, 1.3519E-11, 1.2992E-11, 1.2563E-11, 1.2223E-11/    F10540
      DATA S1301/                                                         F10550
     *     1.1962E-11, 1.1775E-11, 1.1657E-11, 1.1605E-11, 1.1619E-11,    F10560
     *     1.1697E-11, 1.1839E-11, 1.2046E-11, 1.2319E-11, 1.2659E-11,    F10570
     *     1.3070E-11, 1.3553E-11, 1.4113E-11, 1.4754E-11, 1.5480E-11,    F10580
     *     1.6298E-11, 1.7214E-11, 1.8236E-11, 1.9372E-11, 2.0635E-11,    F10590
     *     2.2036E-11, 2.3590E-11, 2.5317E-11, 2.7242E-11, 2.9400E-11,    F10600
     *     3.1849E-11, 3.4654E-11, 3.7923E-11, 4.1695E-11, 4.6055E-11,    F10610
     *     5.0940E-11, 5.5624E-11, 6.0667E-11, 6.6261E-11, 7.2692E-11,    F10620
     *     7.9711E-11, 8.7976E-11, 9.6884E-11, 1.0775E-10, 1.2093E-10,    F10630
     *     1.3531E-10, 1.5404E-10, 1.7315E-10, 1.9862E-10, 2.3341E-10,    F10640
     *     2.7014E-10, 3.1716E-10, 3.6957E-10, 4.3233E-10, 5.2566E-10/    F10650
      DATA S1351/                                                         F10660
     *     6.2251E-10, 7.2149E-10, 8.3958E-10, 9.5931E-10, 1.1388E-09,    F10670
     *     1.2973E-09, 1.4442E-09, 1.5638E-09, 1.6974E-09, 1.8489E-09,    F10680
     *     1.9830E-09, 2.1720E-09, 2.3662E-09, 2.6987E-09, 3.1697E-09,    F10690
     *     3.6907E-09, 4.2625E-09, 4.7946E-09, 5.3848E-09, 6.0897E-09,    F10700
     *     6.4730E-09, 7.1483E-09, 7.7432E-09, 8.0851E-09, 8.5013E-09,    F10710
     *     8.5909E-09, 9.1890E-09, 9.3124E-09, 9.5936E-09, 9.8787E-09,    F10720
     *     9.9036E-09, 9.6712E-09, 9.2036E-09, 9.0466E-09, 8.9380E-09,    F10730
     *     9.1815E-09, 9.5092E-09, 1.0027E-08, 1.0876E-08, 1.1744E-08,    F10740
     *     1.1853E-08, 1.1296E-08, 1.0134E-08, 8.8245E-09, 7.3930E-09,    F10750
     *     5.7150E-09, 4.4884E-09, 3.4027E-09, 2.6054E-09, 2.0790E-09/    F10760
      DATA S1401/                                                         F10770
     *     1.7267E-09, 1.4724E-09, 1.2722E-09, 1.1234E-09, 1.0186E-09,    F10780
     *     9.4680E-10, 8.8854E-10, 8.5127E-10, 8.3157E-10, 8.2226E-10,    F10790
     *     8.3395E-10, 8.3294E-10, 8.4725E-10, 8.8814E-10, 9.3697E-10,    F10800
     *     1.0112E-09, 1.0412E-09, 1.0948E-09, 1.1810E-09, 1.2267E-09,    F10810
     *     1.3690E-09, 1.4512E-09, 1.5568E-09, 1.6552E-09, 1.7321E-09,    F10820
     *     1.8797E-09, 1.9210E-09, 1.9686E-09, 1.9917E-09, 1.9357E-09,    F10830
     *     1.8486E-09, 1.7575E-09, 1.7113E-09, 1.7163E-09, 1.7623E-09,    F10840
     *     1.8536E-09, 1.9765E-09, 2.1334E-09, 2.3237E-09, 2.3259E-09,    F10850
     *     2.1833E-09, 1.9785E-09, 1.7308E-09, 1.4596E-09, 1.1198E-09,    F10860
     *     8.7375E-10, 6.5381E-10, 4.8677E-10, 3.6756E-10, 2.9155E-10/    F10870
      DATA S1451/                                                         F10880
     *     2.3735E-10, 1.9590E-10, 1.6638E-10, 1.4549E-10, 1.2947E-10,    F10890
     *     1.1511E-10, 1.0548E-10, 9.6511E-11, 9.0469E-11, 8.5170E-11,    F10900
     *     7.7804E-11, 7.1971E-11, 6.6213E-11, 6.1063E-11, 5.5881E-11,    F10910
     *     5.0508E-11, 4.5932E-11, 4.1997E-11, 3.7672E-11, 3.3972E-11,    F10920
     *     3.0318E-11, 2.6769E-11, 2.3874E-11, 2.1336E-11, 1.9073E-11,    F10930
     *     1.7313E-11, 1.5904E-11, 1.4684E-11, 1.3698E-11, 1.2873E-11,    F10940
     *     1.2175E-11, 1.1542E-11, 1.1024E-11, 1.0602E-11, 1.0267E-11,    F10950
     *     1.0012E-11, 9.8379E-12, 9.7482E-12, 9.7564E-12, 9.8613E-12,    F10960
     *     1.0092E-11, 1.0418E-11, 1.0868E-11, 1.1585E-11, 1.2351E-11,    F10970
     *     1.3372E-11, 1.4841E-11, 1.6457E-11, 1.8681E-11, 2.0550E-11/    F10980
      DATA S1501/                                                         F10990
     *     2.2912E-11, 2.5958E-11, 2.9137E-11, 3.2368E-11, 3.4848E-11,    F11000
     *     3.8462E-11, 4.2190E-11, 4.5629E-11, 4.9022E-11, 5.4232E-11,    F11010
     *     6.1900E-11, 7.1953E-11, 8.5368E-11, 9.9699E-11, 1.1734E-10,    F11020
     *     1.4185E-10, 1.7017E-10, 1.9813E-10, 2.3859E-10, 2.7304E-10,    F11030
     *     3.0971E-10, 3.5129E-10, 3.9405E-10, 4.5194E-10, 4.8932E-10,    F11040
     *     5.2436E-10, 5.4098E-10, 5.5542E-10, 5.7794E-10, 5.6992E-10,    F11050
     *     5.8790E-10, 6.1526E-10, 6.8034E-10, 6.7956E-10, 6.6864E-10,    F11060
     *     6.9329E-10, 7.2971E-10, 7.6546E-10, 7.5078E-10, 7.8406E-10,    F11070
     *     8.3896E-10, 9.0111E-10, 9.1994E-10, 8.7189E-10, 8.1426E-10,    F11080
     *     7.3097E-10, 6.3357E-10, 5.1371E-10, 4.0936E-10, 3.2918E-10/    F11090
      DATA S1551/                                                         F11100
     *     2.6255E-10, 2.0724E-10, 1.6879E-10, 1.4165E-10, 1.1989E-10,    F11110
     *     1.0125E-10, 8.9629E-11, 7.8458E-11, 6.8826E-11, 6.0935E-11,    F11120
     *     5.5208E-11, 5.2262E-11, 5.0260E-11, 4.8457E-11, 4.7888E-11,    F11130
     *     4.8032E-11, 5.0838E-11, 5.4668E-11, 5.5790E-11, 6.0056E-11,    F11140
     *     6.3811E-11, 6.8848E-11, 7.4590E-11, 7.8249E-11, 8.3371E-11,    F11150
     *     8.3641E-11, 8.6591E-11, 8.9599E-11, 9.3487E-11, 1.0066E-10,    F11160
     *     1.0765E-10, 1.0851E-10, 1.0619E-10, 1.0557E-10, 1.0460E-10,    F11170
     *     1.0796E-10, 1.0523E-10, 1.0674E-10, 1.1261E-10, 1.1431E-10,    F11180
     *     1.1408E-10, 1.0901E-10, 9.9105E-11, 8.8077E-11, 6.9928E-11,    F11190
     *     5.4595E-11, 4.5401E-11, 3.6313E-11, 2.6986E-11, 1.9463E-11/    F11200
      DATA S1601/                                                         F11210
     *     1.4577E-11, 1.1583E-11, 9.5492E-12, 8.0770E-12, 6.9642E-12,    F11220
     *     6.0966E-12, 5.4046E-12, 4.8431E-12, 4.3815E-12, 3.9987E-12,    F11230
     *     3.6790E-12, 3.4113E-12, 3.1868E-12, 2.9992E-12, 2.8434E-12,    F11240
     *     2.7153E-12, 2.6120E-12, 2.5311E-12, 2.4705E-12, 2.4290E-12,    F11250
     *     2.4053E-12, 2.3988E-12, 2.4087E-12, 2.4349E-12, 2.4771E-12,    F11260
     *     2.5355E-12, 2.6103E-12, 2.7019E-12, 2.8110E-12, 2.9383E-12,    F11270
     *     3.0848E-12, 3.2518E-12, 3.4405E-12, 3.6527E-12, 3.8902E-12,    F11280
     *     4.1555E-12, 4.4510E-12, 4.7801E-12, 5.1462E-12, 5.5539E-12,    F11290
     *     6.0086E-12, 6.5171E-12, 7.0884E-12, 7.7357E-12, 8.4831E-12,    F11300
     *     9.3096E-12, 1.0282E-11, 1.1407E-11, 1.2690E-11, 1.4148E-11/    F11310
      DATA S1651/                                                         F11320
     *     1.5888E-11, 1.7992E-11, 2.0523E-11, 2.3342E-11, 2.6578E-11,    F11330
     *     3.0909E-11, 3.6228E-11, 4.2053E-11, 4.9059E-11, 5.9273E-11,    F11340
     *     7.0166E-11, 8.2298E-11, 9.7071E-11, 1.1673E-10, 1.4010E-10,    F11350
     *     1.6621E-10, 2.0127E-10, 2.3586E-10, 2.7050E-10, 3.0950E-10,    F11360
     *     3.6584E-10, 4.1278E-10, 4.6591E-10, 5.2220E-10, 5.5246E-10,    F11370
     *     6.1500E-10, 6.5878E-10, 7.1167E-10, 7.9372E-10, 8.6975E-10,    F11380
     *     9.6459E-10, 9.7368E-10, 9.8142E-10, 1.0202E-09, 1.0200E-09,    F11390
     *     1.0356E-09, 1.0092E-09, 1.0269E-09, 1.0366E-09, 1.0490E-09,    F11400
     *     1.0717E-09, 1.0792E-09, 1.1016E-09, 1.0849E-09, 1.0929E-09,    F11410
     *     1.0971E-09, 1.0969E-09, 1.0460E-09, 9.2026E-10, 8.1113E-10/    F11420
      DATA S1701/                                                         F11430
     *     6.8635E-10, 5.5369E-10, 4.2908E-10, 3.3384E-10, 2.6480E-10,    F11440
     *     2.0810E-10, 1.6915E-10, 1.4051E-10, 1.1867E-10, 1.0158E-10,    F11450
     *     8.8990E-11, 7.9175E-11, 7.0440E-11, 6.3453E-11, 5.7009E-11,    F11460
     *     5.1662E-11, 4.7219E-11, 4.3454E-11, 4.0229E-11, 3.7689E-11,    F11470
     *     3.6567E-11, 3.5865E-11, 3.5955E-11, 3.5928E-11, 3.6298E-11,    F11480
     *     3.7629E-11, 3.9300E-11, 4.1829E-11, 4.4806E-11, 5.0534E-11,    F11490
     *     5.6672E-11, 6.2138E-11, 6.8678E-11, 7.6111E-11, 8.4591E-11,    F11500
     *     9.2634E-11, 9.8085E-11, 1.0830E-10, 1.1949E-10, 1.2511E-10,    F11510
     *     1.3394E-10, 1.3505E-10, 1.4342E-10, 1.4874E-10, 1.4920E-10,    F11520
     *     1.5872E-10, 1.5972E-10, 1.5821E-10, 1.5425E-10, 1.4937E-10/    F11530
      DATA S1751/                                                         F11540
     *     1.5089E-10, 1.5521E-10, 1.6325E-10, 1.6924E-10, 1.8265E-10,    F11550
     *     1.9612E-10, 2.0176E-10, 1.9359E-10, 1.7085E-10, 1.5197E-10,    F11560
     *     1.2646E-10, 9.8552E-11, 7.4530E-11, 5.5052E-11, 4.2315E-11,    F11570
     *     3.2736E-11, 2.6171E-11, 2.1909E-11, 1.8286E-11, 1.5752E-11,    F11580
     *     1.3859E-11, 1.2288E-11, 1.1002E-11, 9.7534E-12, 8.8412E-12,    F11590
     *     8.0169E-12, 7.2855E-12, 6.8734E-12, 6.4121E-12, 6.1471E-12,    F11600
     *     5.7780E-12, 5.3478E-12, 4.9652E-12, 4.4043E-12, 3.9862E-12,    F11610
     *     3.4684E-12, 2.9681E-12, 2.5791E-12, 2.2339E-12, 1.9247E-12,    F11620
     *     1.6849E-12, 1.4863E-12, 1.3291E-12, 1.2021E-12, 1.0947E-12,    F11630
     *     1.0015E-12, 9.1935E-13, 8.4612E-13, 7.8036E-13, 7.2100E-13/    F11640
      DATA S1801/                                                         F11650
     *     6.6718E-13, 6.1821E-13, 5.7353E-13, 5.3269E-13, 4.9526E-13,    F11660
     *     4.6093E-13, 4.2937E-13, 4.0034E-13, 3.7361E-13, 3.4895E-13,    F11670
     *     3.2621E-13, 3.0520E-13, 2.8578E-13, 2.6782E-13, 2.5120E-13,    F11680
     *     2.3581E-13, 2.2154E-13, 2.0832E-13, 1.9605E-13, 1.8466E-13,    F11690
     *     1.7408E-13, 1.6425E-13, 1.5511E-13, 1.4661E-13, 1.3869E-13,    F11700
     *     1.3131E-13, 1.2444E-13, 1.1803E-13, 1.1205E-13, 1.0646E-13,    F11710
     *     1.0124E-13, 9.6358E-14, 9.1789E-14, 8.7509E-14, 8.3498E-14,    F11720
     *     7.9735E-14, 7.6202E-14, 7.2882E-14, 6.9760E-14, 6.6822E-14,    F11730
     *     6.4053E-14, 6.1442E-14, 5.8978E-14, 5.6650E-14, 5.4448E-14,    F11740
     *     5.2364E-14, 5.0389E-14, 4.8516E-14, 4.6738E-14, 4.5048E-14/    F11750
      DATA S1851/                                                         F11760
     *     4.3441E-14, 4.1911E-14, 4.0453E-14, 3.9063E-14, 3.7735E-14,    F11770
     *     3.6467E-14, 3.5254E-14, 3.4093E-14, 3.2980E-14, 3.1914E-14,    F11780
     *     3.0891E-14, 2.9909E-14, 2.8965E-14, 2.8058E-14, 2.7185E-14,    F11790
     *     2.6344E-14, 2.5535E-14, 2.4755E-14, 2.4002E-14, 2.3276E-14,    F11800
     *     2.2576E-14, 2.1899E-14, 2.1245E-14, 2.0613E-14, 2.0002E-14,    F11810
     *     1.9411E-14, 1.8839E-14, 1.8285E-14, 1.7749E-14, 1.7230E-14,    F11820
     *     1.6727E-14, 1.6240E-14, 1.5768E-14, 1.5310E-14, 1.4867E-14,    F11830
     *     1.4436E-14, 1.4019E-14, 1.3614E-14, 1.3221E-14, 1.2840E-14,    F11840
     *     1.2471E-14, 1.2112E-14, 1.1764E-14, 1.1425E-14, 1.1097E-14,    F11850
     *     1.0779E-14, 1.0469E-14, 1.0169E-14, 9.8775E-15, 9.5943E-15/    F11860
      DATA S1901/                                                         F11870
     *     9.3193E-15, 9.0522E-15, 8.7928E-15, 8.5409E-15, 8.2962E-15,    F11880
     *     8.0586E-15, 7.8278E-15, 7.6036E-15, 7.3858E-15, 7.1742E-15,    F11890
     *     6.9687E-15, 6.7691E-15, 6.5752E-15, 6.3868E-15, 6.2038E-15,    F11900
     *     6.0260E-15, 5.8533E-15, 5.6856E-15, 5.5226E-15, 5.3642E-15,    F11910
     *     5.2104E-15, 5.0610E-15, 4.9158E-15, 4.7748E-15, 4.6378E-15,    F11920
     *     4.5047E-15, 4.3753E-15, 4.2497E-15, 4.1277E-15, 4.0091E-15,    F11930
     *     3.8939E-15, 3.7820E-15, 3.6733E-15, 3.5677E-15, 3.4651E-15,    F11940
     *     3.3655E-15, 3.2686E-15, 3.1746E-15, 3.0832E-15, 2.9944E-15,    F11950
     *     2.9082E-15, 2.8244E-15, 2.7431E-15, 2.6640E-15, 2.5872E-15,    F11960
     *     2.5126E-15, 2.4401E-15, 2.3697E-15, 2.3014E-15, 2.2349E-15/    F11970
      DATA S1951/                                                         F11980
     *     2.1704E-15, 2.1077E-15, 2.0468E-15, 1.9877E-15, 1.9302E-15,    F11990
     *     1.8744E-15, 1.8202E-15, 1.7675E-15, 1.7164E-15, 1.6667E-15,    F12000
     *     1.6184E-15, 1.5716E-15, 1.5260E-15, 1.4818E-15, 1.4389E-15,    F12010
     *     1.3971E-15, 1.3566E-15, 1.3172E-15, 1.2790E-15, 1.2419E-15,    F12020
     *     1.2058E-15, 1.1708E-15, 1.1368E-15, 1.1037E-15, 1.0716E-15,    F12030
     *     1.0405E-15, 1.0102E-15, 9.8079E-16, 9.5224E-16, 9.2451E-16,    F12040
     *     8.9758E-16, 8.7142E-16, 8.4602E-16, 8.2136E-16, 7.9740E-16,    F12050
     *     7.7414E-16, 7.5154E-16, 7.2961E-16, 7.0830E-16, 6.8761E-16,    F12060
     *     6.6752E-16, 6.4801E-16, 6.2906E-16, 6.1066E-16, 5.9280E-16,    F12070
     *     5.7545E-16, 5.5860E-16, 5.4224E-16, 5.2636E-16, 5.1094E-16/    F12080
      DATA S2001/                                                         F12090
     *     4.9596E-16/                                                    F12100
C                                                                         F12110
      END                                                                 F12120
C
C     --------------------------------------------------------------
C
      SUBROUTINE FRN296 (V1C,V2C,DVC,NPTC,C)                              F12130
C                                                                         F12140
      IMPLICIT REAL*8           (V)                                     ! F12150
C                                                                         F12160
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F12170
      COMMON /FH2O/ V1S,V2S,DVS,NPTS,S(2003)                              F12180
      DIMENSION C(*)                                                      F12190
C                                                                         F12200
      DVC = DVS                                                           F12210
      V1C = V1ABS-DVC                                                     F12220
      V2C = V2ABS+DVC                                                     F12230
C                                                                         F12240
      I1 = (V1C-V1S)/DVS                                                  F12250
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F12320
         I = I1+J                                                         F12330
         C(J) = 0.                                                        F12340
         IF ((I.GE.1).AND.(I.LE.NPTS)) THEN                               F12350
            C(J) = S(I)                                                   F12360
         ENDIF                                                            F12370
   10 CONTINUE                                                            F12380
C                                                                         F12390
      RETURN                                                              F12400
C                                                                         F12410
      END                                                                 F12420
C
C     --------------------------------------------------------------
C
      BLOCK DATA BFH2O                                                    F12430
C                                                                         F12440

c     Adjustments made to coefficients from 780 - 900 cm-1 to match
c     foreign component of CKD_2.3 (as listed in block data BFH2Oa).

      IMPLICIT REAL*8           (V)                                     ! F12450
C                                                                         F12460
C               06/28/82                                                  F12470
C               UNITS OF (CM**3/MOL)*1.E-20                               F12480
C                                                                         F12490
      COMMON /FH2O/ V1,V2,DV,NPT,                                         F12500
     *              F0000( 2),F0001(50),F0051(50),F0101(50),F0151(50),    F12510
     *              F0201(50),F0251(50),F0301(50),F0351(50),F0401(50),    F12520
     *              F0451(50),F0501(50),F0551(50),F0601(50),F0651(50),    F12530
     *              F0701(50),F0751(50),F0801(50),F0851(50),F0901(50),    F12540
     *              F0951(50),F1001(50),F1051(50),F1101(50),F1151(50),    F12550
     *              F1201(50),F1251(50),F1301(50),F1351(50),F1401(50),    F12560
     *              F1451(50),F1501(50),F1551(50),F1601(50),F1651(50),    F12570
     *              F1701(50),F1751(50),F1801(50),F1851(50),F1901(50),    F12580
     *              F1951(50),F2001(1)                                    F12590
C                                                                         F12600
       DATA V1,V2,DV,NPT / -20.0, 20000.0, 10.0, 2003/                    F12610
C                                                                         F12620
      DATA F0000/                                                         F12630
     *     1.2859E-02, 1.1715E-02/                                        F12640
      DATA F0001/                                                         F12650
     *     1.1038E-02, 1.1715E-02, 1.2859E-02, 1.5326E-02, 1.6999E-02,    F12660
     *     1.8321E-02, 1.9402E-02, 1.9570E-02, 1.9432E-02, 1.7572E-02,    F12670
     *     1.6760E-02, 1.5480E-02, 1.3984E-02, 1.2266E-02, 1.0467E-02,    F12680
     *     9.4526E-03, 8.0485E-03, 6.9484E-03, 6.1416E-03, 5.0941E-03,    F12690
     *     4.4836E-03, 3.8133E-03, 3.4608E-03, 3.1487E-03, 2.4555E-03,    F12700
     *     2.0977E-03, 1.7266E-03, 1.4920E-03, 1.2709E-03, 9.8081E-04,    F12710
     *     8.5063E-04, 6.8822E-04, 5.3809E-04, 4.4679E-04, 3.3774E-04,    F12720
     *     2.7979E-04, 2.1047E-04, 1.6511E-04, 1.2993E-04, 9.3033E-05,    F12730
     *     7.4360E-05, 5.6428E-05, 4.5442E-05, 3.4575E-05, 2.7903E-05,    F12740
     *     2.1374E-05, 1.6075E-05, 1.3022E-05, 1.0962E-05, 8.5959E-06/    F12750
      DATA F0051/                                                         F12760
     *     6.9125E-06, 5.3808E-06, 4.3586E-06, 3.6394E-06, 2.9552E-06,    F12770
     *     2.3547E-06, 1.8463E-06, 1.6036E-06, 1.3483E-06, 1.1968E-06,    F12780
     *     1.0333E-06, 8.4484E-07, 6.7195E-07, 5.0947E-07, 4.2343E-07,    F12790
     *     3.4453E-07, 2.7830E-07, 2.3063E-07, 1.9951E-07, 1.7087E-07,    F12800
     *     1.4393E-07, 1.2575E-07, 1.0750E-07, 8.2325E-08, 5.7524E-08,    F12810
     *     4.4482E-08, 3.8106E-08, 3.4315E-08, 2.9422E-08, 2.5069E-08,    F12820
     *     2.2402E-08, 1.9349E-08, 1.6152E-08, 1.2208E-08, 8.9660E-09,    F12830
     *     7.1322E-09, 6.1028E-09, 5.2938E-09, 4.5350E-09, 3.4977E-09,    F12840
     *     2.9511E-09, 2.4734E-09, 2.0508E-09, 1.8507E-09, 1.6373E-09,    F12850
     *     1.5171E-09, 1.3071E-09, 1.2462E-09, 1.2148E-09, 1.2590E-09/    F12860
      DATA F0101/                                                         F12870
     *     1.3153E-09, 1.3301E-09, 1.4483E-09, 1.6944E-09, 2.0559E-09,    F12880
     *     2.2954E-09, 2.6221E-09, 3.2606E-09, 4.2392E-09, 5.2171E-09,    F12890
     *     6.2553E-09, 8.2548E-09, 9.5842E-09, 1.1280E-08, 1.3628E-08,    F12900
     *     1.7635E-08, 2.1576E-08, 2.4835E-08, 3.0014E-08, 3.8485E-08,    F12910
     *     4.7440E-08, 5.5202E-08, 7.0897E-08, 9.6578E-08, 1.3976E-07,    F12920
     *     1.8391E-07, 2.3207E-07, 2.9960E-07, 4.0408E-07, 5.9260E-07,    F12930
     *     7.8487E-07, 1.0947E-06, 1.4676E-06, 1.9325E-06, 2.6587E-06,    F12940
     *     3.4534E-06, 4.4376E-06, 5.8061E-06, 7.0141E-06, 8.4937E-06,    F12950
     *     1.0186E-05, 1.2034E-05, 1.3837E-05, 1.6595E-05, 1.9259E-05,    F12960
     *     2.1620E-05, 2.3681E-05, 2.7064E-05, 3.2510E-05, 3.5460E-05/    F12970
      DATA F0151/                                                         F12980
     *     3.9109E-05, 4.2891E-05, 4.7757E-05, 5.0981E-05, 5.0527E-05,    F12990
     *     4.8618E-05, 4.4001E-05, 3.7982E-05, 3.2667E-05, 2.7794E-05,    F13000
     *     2.4910E-05, 2.4375E-05, 2.7316E-05, 3.2579E-05, 3.5499E-05,    F13010
     *     3.8010E-05, 4.1353E-05, 4.3323E-05, 4.3004E-05, 3.9790E-05,    F13020
     *     3.7718E-05, 3.6360E-05, 3.2386E-05, 2.7409E-05, 2.3626E-05,    F13030
     *     2.0631E-05, 1.8371E-05, 1.5445E-05, 1.2989E-05, 1.1098E-05,    F13040
     *     9.6552E-06, 8.0649E-06, 7.2365E-06, 5.9137E-06, 5.2759E-06,    F13050
     *     4.8860E-06, 4.1321E-06, 3.5918E-06, 2.7640E-06, 2.4892E-06,    F13060
     *     2.1018E-06, 1.7848E-06, 1.5855E-06, 1.3569E-06, 1.1986E-06,    F13070
     *     9.4693E-07, 7.4097E-07, 6.3443E-07, 4.8131E-07, 4.0942E-07/    F13080
      DATA F0201/                                                         F13090
     *     3.3316E-07, 2.8488E-07, 2.3461E-07, 1.7397E-07, 1.4684E-07,    F13100
     *     1.0953E-07, 8.5396E-08, 6.9261E-08, 5.4001E-08, 4.5430E-08,    F13110
     *     3.2791E-08, 2.5995E-08, 2.0225E-08, 1.5710E-08, 1.3027E-08,    F13120
     *     1.0229E-08, 8.5277E-09, 6.5249E-09, 5.0117E-09, 3.9906E-09,    F13130
     *     3.2332E-09, 2.7847E-09, 2.4570E-09, 2.3359E-09, 2.0599E-09,    F13140
     *     1.8436E-09, 1.6559E-09, 1.4910E-09, 1.2794E-09, 9.8229E-10,    F13150
     *     8.0054E-10, 6.0769E-10, 4.5646E-10, 3.3111E-10, 2.4428E-10,    F13160
     *     1.8007E-10, 1.3291E-10, 9.7974E-11, 7.8271E-11, 6.3833E-11,    F13170
     *     5.4425E-11, 4.6471E-11, 4.0209E-11, 3.5227E-11, 3.1212E-11,    F13180
     *     2.8840E-11, 2.7762E-11, 2.7935E-11, 3.2012E-11, 3.9525E-11/    F13190
      DATA F0251/                                                         F13200
     *     5.0303E-11, 6.8027E-11, 9.3954E-11, 1.2986E-10, 1.8478E-10,    F13210
     *     2.5331E-10, 3.4827E-10, 4.6968E-10, 6.2380E-10, 7.9106E-10,    F13220
     *     1.0026E-09, 1.2102E-09, 1.4146E-09, 1.6154E-09, 1.7510E-09,    F13230
     *     1.8575E-09, 1.8742E-09, 1.8700E-09, 1.8582E-09, 1.9657E-09,    F13240
     *     2.1204E-09, 2.0381E-09, 2.0122E-09, 2.0436E-09, 2.1213E-09,    F13250
     *     2.0742E-09, 1.9870E-09, 2.0465E-09, 2.1556E-09, 2.2222E-09,    F13260
     *     2.1977E-09, 2.1047E-09, 1.9334E-09, 1.7357E-09, 1.5754E-09,    F13270
     *     1.4398E-09, 1.4018E-09, 1.5459E-09, 1.7576E-09, 2.1645E-09,    F13280
     *     2.9480E-09, 4.4439E-09, 5.8341E-09, 8.0757E-09, 1.1658E-08,    F13290
     *     1.6793E-08, 2.2694E-08, 2.9468E-08, 3.9278E-08, 5.2145E-08/    F13300
      DATA F0301/                                                         F13310
     *     6.4378E-08, 7.7947E-08, 8.5321E-08, 9.7848E-08, 1.0999E-07,    F13320
     *     1.1489E-07, 1.2082E-07, 1.2822E-07, 1.4053E-07, 1.5238E-07,    F13330
     *     1.5454E-07, 1.5018E-07, 1.4048E-07, 1.2359E-07, 1.0858E-07,    F13340
     *     9.3486E-08, 8.1638E-08, 7.7690E-08, 8.4625E-08, 1.0114E-07,    F13350
     *     1.1430E-07, 1.2263E-07, 1.3084E-07, 1.3380E-07, 1.3573E-07,    F13360
     *     1.3441E-07, 1.2962E-07, 1.2638E-07, 1.1934E-07, 1.1371E-07,    F13370
     *     1.0871E-07, 9.8843E-08, 9.1877E-08, 9.1050E-08, 9.3213E-08,    F13380
     *     9.2929E-08, 1.0155E-07, 1.1263E-07, 1.2370E-07, 1.3636E-07,    F13390
     *     1.5400E-07, 1.7656E-07, 2.1329E-07, 2.3045E-07, 2.5811E-07,    F13400
     *     2.9261E-07, 3.4259E-07, 4.0770E-07, 4.8771E-07, 5.8081E-07/    F13410
      DATA F0351/                                                         F13420
     *     7.2895E-07, 8.7482E-07, 1.0795E-06, 1.3384E-06, 1.7208E-06,    F13430
     *     2.0677E-06, 2.5294E-06, 3.1123E-06, 3.7900E-06, 4.7752E-06,    F13440
     *     5.6891E-06, 6.6261E-06, 7.6246E-06, 8.7730E-06, 9.6672E-06,    F13450
     *     1.0980E-05, 1.1287E-05, 1.1670E-05, 1.1635E-05, 1.1768E-05,    F13460
     *     1.2039E-05, 1.2253E-05, 1.3294E-05, 1.4005E-05, 1.3854E-05,    F13470
     *     1.3420E-05, 1.3003E-05, 1.2645E-05, 1.1715E-05, 1.1258E-05,    F13480
     *     1.1516E-05, 1.2494E-05, 1.3655E-05, 1.4931E-05, 1.4649E-05,    F13490
     *     1.3857E-05, 1.3120E-05, 1.1791E-05, 1.0637E-05, 8.2760E-06,    F13500
     *     6.5821E-06, 5.1959E-06, 4.0158E-06, 3.0131E-06, 2.0462E-06,    F13510
     *     1.4853E-06, 1.0365E-06, 7.3938E-07, 4.9752E-07, 3.4148E-07/    F13520
      DATA F0401/                                                         F13530
     *     2.4992E-07, 1.8363E-07, 1.4591E-07, 1.1380E-07, 9.0588E-08,    F13540
     *     7.3697E-08, 6.0252E-08, 5.1868E-08, 4.2660E-08, 3.6163E-08,    F13550
     *     3.2512E-08, 2.9258E-08, 2.4238E-08, 2.1209E-08, 1.6362E-08,    F13560
     *     1.3871E-08, 1.2355E-08, 9.6940E-09, 7.7735E-09, 6.2278E-09,    F13570
     *     5.2282E-09, 4.3799E-09, 3.5545E-09, 2.7527E-09, 2.0950E-09,    F13580
     *     1.6344E-09, 1.2689E-09, 1.0403E-09, 8.4880E-10, 6.3461E-10,    F13590
     *     4.7657E-10, 3.5220E-10, 2.7879E-10, 2.3021E-10, 1.6167E-10,    F13600
     *     1.1732E-10, 8.9206E-11, 7.0596E-11, 5.8310E-11, 4.4084E-11,    F13610
     *     3.1534E-11, 2.5068E-11, 2.2088E-11, 2.2579E-11, 2.2637E-11,    F13620
     *     2.5705E-11, 3.2415E-11, 4.6116E-11, 6.5346E-11, 9.4842E-11/    F13630
      DATA F0451/                                                         F13640
     *     1.2809E-10, 1.8211E-10, 2.4052E-10, 3.0270E-10, 3.5531E-10,    F13650
     *     4.2402E-10, 4.6730E-10, 4.7942E-10, 4.6813E-10, 4.5997E-10,    F13660
     *     4.5788E-10, 4.0311E-10, 3.7367E-10, 3.3149E-10, 2.9281E-10,    F13670
     *     2.5231E-10, 2.1152E-10, 1.9799E-10, 1.8636E-10, 1.9085E-10,    F13680
     *     2.0786E-10, 2.2464E-10, 2.3785E-10, 2.5684E-10, 2.7499E-10,    F13690
     *     2.6962E-10, 2.6378E-10, 2.6297E-10, 2.6903E-10, 2.7035E-10,    F13700
     *     2.5394E-10, 2.5655E-10, 2.7184E-10, 2.9013E-10, 3.0585E-10,    F13710
     *     3.0791E-10, 3.1667E-10, 3.4343E-10, 3.7365E-10, 4.0269E-10,    F13720
     *     4.7260E-10, 5.6584E-10, 6.9791E-10, 8.6569E-10, 1.0393E-09,    F13730
     *     1.2067E-09, 1.5047E-09, 1.8583E-09, 2.2357E-09, 2.6498E-09/    F13740
      DATA F0501/                                                         F13750
     *     3.2483E-09, 3.9927E-09, 4.6618E-09, 5.5555E-09, 6.6609E-09,    F13760
     *     8.2139E-09, 1.0285E-08, 1.3919E-08, 1.8786E-08, 2.5150E-08,    F13770
     *     3.3130E-08, 4.5442E-08, 6.3370E-08, 9.0628E-08, 1.2118E-07,    F13780
     *     1.5927E-07, 2.1358E-07, 2.7825E-07, 3.7671E-07, 4.4894E-07,    F13790
     *     5.4442E-07, 6.2240E-07, 7.3004E-07, 8.3384E-07, 8.7933E-07,    F13800
     *     8.8080E-07, 8.6939E-07, 8.6541E-07, 8.2055E-07, 7.7278E-07,    F13810
     *     7.5989E-07, 8.6909E-07, 9.7945E-07, 1.0394E-06, 1.0646E-06,    F13820
     *     1.1509E-06, 1.2017E-06, 1.1915E-06, 1.1259E-06, 1.1549E-06,    F13830
     *     1.1938E-06, 1.2356E-06, 1.2404E-06, 1.1716E-06, 1.1149E-06,    F13840
     *     1.0073E-06, 8.9845E-07, 7.6639E-07, 6.1517E-07, 5.0887E-07/    F13850
      DATA F0551/                                                         F13860
     *     4.1269E-07, 3.2474E-07, 2.5698E-07, 1.8893E-07, 1.4009E-07,    F13870
     *     1.0340E-07, 7.7724E-08, 5.7302E-08, 4.2178E-08, 2.9603E-08,    F13880
     *     2.1945E-08, 1.6301E-08, 1.2806E-08, 1.0048E-08, 7.8970E-09,    F13890
     *     6.1133E-09, 4.9054E-09, 4.1985E-09, 3.6944E-09, 3.2586E-09,    F13900
     *     2.7362E-09, 2.3647E-09, 2.1249E-09, 1.8172E-09, 1.6224E-09,    F13910
     *     1.5158E-09, 1.2361E-09, 1.0682E-09, 9.2312E-10, 7.9220E-10,    F13920
     *     6.8174E-10, 5.6147E-10, 4.8268E-10, 4.1534E-10, 3.3106E-10,    F13930
     *     2.8275E-10, 2.4584E-10, 2.0742E-10, 1.7840E-10, 1.4664E-10,    F13940
     *     1.2390E-10, 1.0497E-10, 8.5038E-11, 6.7008E-11, 5.6355E-11,    F13950
     *     4.3323E-11, 3.6914E-11, 3.2262E-11, 3.0749E-11, 3.0318E-11/    F13960
      DATA F0601/                                                         F13970
     *     2.9447E-11, 2.9918E-11, 3.0668E-11, 3.1315E-11, 3.0329E-11,    F13980
     *     2.8259E-11, 2.6065E-11, 2.3578E-11, 2.0469E-11, 1.6908E-11,    F13990
     *     1.4912E-11, 1.1867E-11, 9.9730E-12, 8.1014E-12, 6.7528E-12,    F14000
     *     6.3133E-12, 5.8599E-12, 6.0145E-12, 6.5105E-12, 7.0537E-12,    F14010
     *     7.4973E-12, 7.8519E-12, 8.5039E-12, 9.1995E-12, 1.0694E-11,    F14020
     *     1.1659E-11, 1.2685E-11, 1.3087E-11, 1.3222E-11, 1.2634E-11,    F14030
     *     1.1077E-11, 9.6259E-12, 8.3202E-12, 7.4857E-12, 6.8069E-12,    F14040
     *     6.7496E-12, 7.3116E-12, 8.0171E-12, 8.6394E-12, 9.2659E-12,    F14050
     *     1.0048E-11, 1.0941E-11, 1.2226E-11, 1.3058E-11, 1.5193E-11,    F14060
     *     1.8923E-11, 2.3334E-11, 2.8787E-11, 3.6693E-11, 4.8295E-11/    F14070
      DATA F0651/                                                         F14080
     *     6.4260E-11, 8.8269E-11, 1.1865E-10, 1.5961E-10, 2.0605E-10,    F14090
     *     2.7349E-10, 3.7193E-10, 4.8216E-10, 6.1966E-10, 7.7150E-10,    F14100
     *     1.0195E-09, 1.2859E-09, 1.6535E-09, 2.0316E-09, 2.3913E-09,    F14110
     *     3.0114E-09, 3.7495E-09, 4.6504E-09, 5.9145E-09, 7.6840E-09,    F14120
     *     1.0304E-08, 1.3010E-08, 1.6441E-08, 2.1475E-08, 2.5892E-08,    F14130
     *     2.9788E-08, 3.3820E-08, 4.0007E-08, 4.4888E-08, 4.5765E-08,    F14140
     *     4.6131E-08, 4.6239E-08, 4.4849E-08, 4.0729E-08, 3.6856E-08,    F14150
     *     3.6164E-08, 3.7606E-08, 4.1457E-08, 4.3750E-08, 5.1150E-08,    F14160
     *     5.6054E-08, 6.1586E-08, 6.4521E-08, 6.6494E-08, 6.9024E-08,    F14170
     *     6.8893E-08, 7.0901E-08, 6.9760E-08, 7.1485E-08, 7.0740E-08/    F14180
      DATA F0701/                                                         F14190
     *     7.3764E-08, 7.6618E-08, 8.4182E-08, 9.3838E-08, 1.0761E-07,    F14200
     *     1.2851E-07, 1.4748E-07, 1.8407E-07, 2.2109E-07, 2.6392E-07,    F14210
     *     2.9887E-07, 3.4493E-07, 4.0336E-07, 4.3551E-07, 4.9231E-07,    F14220
     *     5.0728E-07, 5.3781E-07, 5.3285E-07, 5.4496E-07, 5.5707E-07,    F14230
     *     5.6944E-07, 6.1123E-07, 6.4317E-07, 6.4581E-07, 6.1999E-07,    F14240
     *     6.0191E-07, 5.7762E-07, 5.7241E-07, 5.7013E-07, 6.0160E-07,    F14250
     *     6.6905E-07, 7.4095E-07, 8.2121E-07, 8.0947E-07, 7.6145E-07,    F14260
     *     7.2193E-07, 6.3722E-07, 5.4316E-07, 4.2186E-07, 3.2528E-07,    F14270
     *     2.5207E-07, 1.8213E-07, 1.2658E-07, 8.6746E-08, 6.0216E-08,    F14280
     *     4.1122E-08, 2.8899E-08, 2.1740E-08, 1.7990E-08, 1.5593E-08/    F14290
      DATA F0751/                                                         F14300
     *     1.3970E-08, 1.2238E-08, 1.0539E-08, 9.2386E-09, 7.8481E-09,    F14310
     *     6.8704E-09, 5.7615E-09, 5.0434E-09, 4.6886E-09, 4.3770E-09,    F14320
     *     3.9768E-09, 3.5202E-09, 3.1854E-09, 2.9009E-09, 2.5763E-09,    F14330
     *     2.2135E-09, 1.9455E-09, 1.6248E-09, 1.3368E-09, 1.0842E-09,    F14340
     *     8.4254E-10, 6.7414E-10, 5.4667E-10, 4.5005E-10, 3.4932E-10,    F14350
     *     2.6745E-10, 2.2053E-10, 1.8162E-10, 1.4935E-10, 1.1618E-10,    F14360
     *     9.1888E-11, 8.0672E-11, 6.8746E-11, 6.2668E-11, 5.5715E-11,    F14370
     *     4.5074E-11, 3.7669E-11, 3.2082E-11, 2.8085E-11, 2.4838E-11,    F14380
     *     1.9791E-11, 1.6964E-11, 1.3887E-11, 1.1179E-11, 9.7499E-12,    F14390
     *     7.8255E-12, 6.3698E-12, 5.3265E-12, 4.6588E-12, 4.4498E-12/    F14400
      DATA F0801/                                                         F14410
     *     3.9984E-12, 3.7513E-12, 3.7176E-12, 3.9148E-12, 4.2702E-12,    F14420
     *     5.0090E-12, 6.5801E-12, 8.7787E-12, 1.2718E-11, 1.8375E-11,    F14430
     *     2.5304E-11, 3.5403E-11, 4.8842E-11, 6.4840E-11, 8.0911E-11,    F14440
     *     1.0136E-10, 1.2311E-10, 1.4203E-10, 1.5869E-10, 1.8093E-10,    F14450
     *     2.1370E-10, 2.5228E-10, 2.8816E-10, 3.4556E-10, 3.9860E-10,    F14460
     *     4.4350E-10, 4.7760E-10, 5.2357E-10, 6.0827E-10, 6.3635E-10,    F14470
     *     6.5886E-10, 6.8753E-10, 7.2349E-10, 7.2789E-10, 6.8232E-10,    F14480
     *     6.6081E-10, 6.4232E-10, 6.3485E-10, 6.4311E-10, 7.2235E-10,    F14490
     *     7.7263E-10, 8.1668E-10, 9.0324E-10, 9.7643E-10, 1.0535E-09,    F14500
     *     1.0195E-09, 1.0194E-09, 1.0156E-09, 9.6792E-10, 9.2725E-10/    F14510
      DATA F0851/                                                         F14520
     *     8.7347E-10, 8.4484E-10, 8.2647E-10, 8.4363E-10, 9.1261E-10,    F14530
     *     1.0051E-09, 1.1511E-09, 1.4037E-09, 1.8066E-09, 2.4483E-09,    F14540
     *     3.2739E-09, 4.3194E-09, 5.6902E-09, 7.7924E-09, 9.7376E-09,    F14550
     *     1.2055E-08, 1.4303E-08, 1.6956E-08, 1.9542E-08, 2.2233E-08,    F14560
     *     2.5186E-08, 2.7777E-08, 2.8943E-08, 2.8873E-08, 2.9417E-08,    F14570
     *     2.7954E-08, 2.7524E-08, 2.7040E-08, 3.1254E-08, 3.6843E-08,    F14580
     *     3.7797E-08, 3.8713E-08, 4.0135E-08, 4.2824E-08, 4.3004E-08,    F14590
     *     4.0279E-08, 4.2781E-08, 4.5220E-08, 4.8948E-08, 5.0172E-08,    F14600
     *     4.8499E-08, 4.7182E-08, 4.2204E-08, 3.7701E-08, 3.0972E-08,    F14610
     *     2.4654E-08, 1.9543E-08, 1.4609E-08, 1.1171E-08, 8.3367E-09/    F14620
      DATA F0901/                                                         F14630
     *     6.3791E-09, 5.0790E-09, 4.0655E-09, 3.3658E-09, 2.7882E-09,    F14640
     *     2.4749E-09, 2.2287E-09, 2.0217E-09, 1.8191E-09, 1.5897E-09,    F14650
     *     1.4191E-09, 1.2448E-09, 1.0884E-09, 9.3585E-10, 7.9429E-10,    F14660
     *     7.3214E-10, 6.5008E-10, 5.7549E-10, 5.4300E-10, 4.7251E-10,    F14670
     *     4.3451E-10, 3.8446E-10, 3.5589E-10, 3.4432E-10, 2.8209E-10,    F14680
     *     2.4620E-10, 2.1278E-10, 1.8406E-10, 1.6314E-10, 1.3261E-10,    F14690
     *     1.1696E-10, 9.6865E-11, 7.6814E-11, 6.6411E-11, 5.0903E-11,    F14700
     *     4.0827E-11, 3.0476E-11, 2.3230E-11, 1.7707E-11, 1.3548E-11,    F14710
     *     1.0719E-11, 9.3026E-12, 8.7967E-12, 8.3136E-12, 7.3918E-12,    F14720
     *     6.5293E-12, 5.9243E-12, 5.3595E-12, 3.5266E-12, 2.2571E-12/    F14730
      DATA F0951/                                                         F14740
     *     1.6150E-12, 1.1413E-12, 8.4998E-13, 7.0803E-13, 5.1747E-13,    F14750
     *     4.0694E-13, 3.6528E-13, 3.3670E-13, 3.1341E-13, 2.9390E-13,    F14760
     *     2.8680E-13, 3.1283E-13, 3.7294E-13, 5.0194E-13, 6.7919E-13,    F14770
     *     1.0455E-12, 1.5230E-12, 2.3932E-12, 3.4231E-12, 5.0515E-12,    F14780
     *     7.3193E-12, 9.9406E-12, 1.2193E-11, 1.4742E-11, 1.9269E-11,    F14790
     *     2.1816E-11, 2.2750E-11, 2.2902E-11, 2.3888E-11, 2.4902E-11,    F14800
     *     2.2160E-11, 2.0381E-11, 1.9903E-11, 2.0086E-11, 1.9304E-11,    F14810
     *     2.0023E-11, 2.2244E-11, 2.5450E-11, 3.1228E-11, 3.4560E-11,    F14820
     *     3.6923E-11, 3.7486E-11, 3.8124E-11, 3.8317E-11, 3.4737E-11,    F14830
     *     3.3037E-11, 3.1724E-11, 2.9840E-11, 2.8301E-11, 2.5857E-11/    F14840
      DATA F1001/                                                         F14850
     *     2.3708E-11, 1.9452E-11, 1.6232E-11, 1.5174E-11, 1.4206E-11,    F14860
     *     1.4408E-11, 1.5483E-11, 1.8642E-11, 2.3664E-11, 3.0181E-11,    F14870
     *     4.0160E-11, 5.2287E-11, 7.2754E-11, 1.0511E-10, 1.4531E-10,    F14880
     *     2.0998E-10, 2.6883E-10, 3.3082E-10, 4.2638E-10, 5.3132E-10,    F14890
     *     6.3617E-10, 7.1413E-10, 8.5953E-10, 9.9715E-10, 1.0796E-09,    F14900
     *     1.0978E-09, 1.1052E-09, 1.1095E-09, 1.0641E-09, 9.7881E-10,    F14910
     *     9.6590E-10, 1.0332E-09, 1.1974E-09, 1.3612E-09, 1.5829E-09,    F14920
     *     1.8655E-09, 2.1465E-09, 2.4779E-09, 2.7370E-09, 2.9915E-09,    F14930
     *     3.3037E-09, 3.6347E-09, 3.9587E-09, 4.4701E-09, 5.0122E-09,    F14940
     *     5.8044E-09, 6.1916E-09, 6.9613E-09, 7.7863E-09, 8.2820E-09/    F14950
      DATA F1051/                                                         F14960
     *     9.4359E-09, 9.7387E-09, 1.0656E-08, 1.0746E-08, 1.1210E-08,    F14970
     *     1.1905E-08, 1.2194E-08, 1.3145E-08, 1.3738E-08, 1.3634E-08,    F14980
     *     1.3011E-08, 1.2511E-08, 1.1805E-08, 1.2159E-08, 1.2390E-08,    F14990
     *     1.3625E-08, 1.5678E-08, 1.7886E-08, 1.9933E-08, 1.9865E-08,    F15000
     *     1.9000E-08, 1.7812E-08, 1.5521E-08, 1.2593E-08, 9.5635E-09,    F15010
     *     7.2987E-09, 5.2489E-09, 3.5673E-09, 2.4206E-09, 1.6977E-09,    F15020
     *     1.2456E-09, 9.3744E-10, 7.8379E-10, 6.9960E-10, 6.6451E-10,    F15030
     *     6.8521E-10, 7.4234E-10, 8.6658E-10, 9.4972E-10, 1.0791E-09,    F15040
     *     1.2359E-09, 1.3363E-09, 1.5025E-09, 1.5368E-09, 1.6152E-09,    F15050
     *     1.6184E-09, 1.6557E-09, 1.7035E-09, 1.6916E-09, 1.7237E-09/    F15060
      DATA F1101/                                                         F15070
     *     1.7175E-09, 1.6475E-09, 1.5335E-09, 1.4272E-09, 1.3282E-09,    F15080
     *     1.3459E-09, 1.4028E-09, 1.5192E-09, 1.7068E-09, 1.9085E-09,    F15090
     *     2.1318E-09, 2.1020E-09, 1.9942E-09, 1.8654E-09, 1.6391E-09,    F15100
     *     1.3552E-09, 1.0186E-09, 7.8540E-10, 5.7022E-10, 3.9247E-10,    F15110
     *     2.5441E-10, 1.6699E-10, 1.1132E-10, 6.8989E-11, 4.5255E-11,    F15120
     *     3.1106E-11, 2.3161E-11, 1.7618E-11, 1.4380E-11, 1.1601E-11,    F15130
     *     9.7148E-12, 8.4519E-12, 6.5392E-12, 5.4113E-12, 4.7624E-12,    F15140
     *     4.0617E-12, 3.6173E-12, 2.8608E-12, 2.2724E-12, 1.7436E-12,    F15150
     *     1.3424E-12, 1.0358E-12, 7.3064E-13, 5.4500E-13, 4.0551E-13,    F15160
     *     2.8642E-13, 2.1831E-13, 1.6860E-13, 1.2086E-13, 1.0150E-13/    F15170
      DATA F1151/                                                         F15180
     *     9.3550E-14, 8.4105E-14, 7.3051E-14, 6.9796E-14, 7.9949E-14,    F15190
     *     1.0742E-13, 1.5639E-13, 2.1308E-13, 3.1226E-13, 4.6853E-13,    F15200
     *     6.6917E-13, 1.0088E-12, 1.4824E-12, 2.2763E-12, 3.3917E-12,    F15210
     *     4.4585E-12, 6.3187E-12, 8.4189E-12, 1.1302E-11, 1.3431E-11,    F15220
     *     1.5679E-11, 1.9044E-11, 2.2463E-11, 2.3605E-11, 2.3619E-11,    F15230
     *     2.3505E-11, 2.3805E-11, 2.2549E-11, 1.9304E-11, 1.8382E-11,    F15240
     *     1.7795E-11, 1.8439E-11, 1.9146E-11, 2.1966E-11, 2.6109E-11,    F15250
     *     3.1883E-11, 3.7872E-11, 4.3966E-11, 4.8789E-11, 5.3264E-11,    F15260
     *     5.9705E-11, 6.3744E-11, 7.0163E-11, 7.9114E-11, 8.8287E-11,    F15270
     *     9.9726E-11, 1.1498E-10, 1.3700E-10, 1.6145E-10, 1.9913E-10/    F15280
      DATA F1201/                                                         F15290
     *     2.2778E-10, 2.6216E-10, 2.9770E-10, 3.3405E-10, 3.7821E-10,    F15300
     *     3.9552E-10, 4.1322E-10, 4.0293E-10, 4.0259E-10, 3.8853E-10,    F15310
     *     3.7842E-10, 3.8551E-10, 4.4618E-10, 5.0527E-10, 5.0695E-10,    F15320
     *     5.1216E-10, 5.1930E-10, 5.5794E-10, 5.3320E-10, 5.2008E-10,    F15330
     *     5.6888E-10, 6.1883E-10, 6.9006E-10, 6.9505E-10, 6.6768E-10,    F15340
     *     6.3290E-10, 5.6753E-10, 5.0327E-10, 3.9830E-10, 3.1147E-10,    F15350
     *     2.4416E-10, 1.8860E-10, 1.3908E-10, 9.9156E-11, 7.3779E-11,    F15360
     *     5.6048E-11, 4.2457E-11, 3.4505E-11, 2.9881E-11, 2.7865E-11,    F15370
     *     2.8471E-11, 3.1065E-11, 3.4204E-11, 3.9140E-11, 4.3606E-11,    F15380
     *     4.9075E-11, 5.3069E-11, 5.5236E-11, 5.5309E-11, 5.3832E-11/    F15390
      DATA F1251/                                                         F15400
     *     5.3183E-11, 5.1783E-11, 5.2042E-11, 5.4422E-11, 5.5656E-11,    F15410
     *     5.4409E-11, 5.2659E-11, 5.1696E-11, 5.1726E-11, 4.9003E-11,    F15420
     *     4.9050E-11, 5.1700E-11, 5.6818E-11, 6.3129E-11, 6.6542E-11,    F15430
     *     6.4367E-11, 5.9908E-11, 5.4470E-11, 4.7903E-11, 3.9669E-11,    F15440
     *     2.9651E-11, 2.2286E-11, 1.6742E-11, 1.1827E-11, 7.7739E-12,    F15450
     *     4.8805E-12, 3.1747E-12, 2.0057E-12, 1.2550E-12, 8.7434E-13,    F15460
     *     6.2755E-13, 4.9752E-13, 4.0047E-13, 3.5602E-13, 3.0930E-13,    F15470
     *     2.4903E-13, 1.9316E-13, 1.4995E-13, 1.2059E-13, 8.7242E-14,    F15480
     *     6.4511E-14, 5.3300E-14, 4.3741E-14, 3.4916E-14, 2.6560E-14,    F15490
     *     1.6923E-14, 1.1816E-14, 6.7071E-15, 3.6474E-15, 2.0686E-15/    F15500
      DATA F1301/                                                         F15510
     *     1.1925E-15, 6.8948E-16, 3.9661E-16, 2.2576E-16, 1.2669E-16,    F15520
     *     6.9908E-17, 3.7896E-17, 2.0280E-17, 1.1016E-17, 6.7816E-18,    F15530
     *     6.0958E-18, 8.9913E-18, 1.7201E-17, 3.4964E-17, 7.0722E-17,    F15540
     *     1.4020E-16, 2.7167E-16, 5.1478E-16, 9.5500E-16, 1.7376E-15,    F15550
     *     3.1074E-15, 5.4789E-15, 9.5640E-15, 1.6635E-14, 2.9145E-14,    F15560
     *     5.2179E-14, 8.8554E-14, 1.4764E-13, 2.3331E-13, 3.5996E-13,    F15570
     *     5.2132E-13, 6.3519E-13, 7.3174E-13, 8.3752E-13, 9.8916E-13,    F15580
     *     1.1515E-12, 1.4034E-12, 1.6594E-12, 2.1021E-12, 2.7416E-12,    F15590
     *     3.4135E-12, 4.5517E-12, 5.5832E-12, 7.2303E-12, 9.9484E-12,    F15600
     *     1.2724E-11, 1.6478E-11, 2.0588E-11, 2.5543E-11, 3.3625E-11/    F15610
      DATA F1351/                                                         F15620
     *     4.1788E-11, 5.0081E-11, 6.0144E-11, 6.9599E-11, 8.4408E-11,    F15630
     *     9.7143E-11, 1.0805E-10, 1.1713E-10, 1.2711E-10, 1.3727E-10,    F15640
     *     1.4539E-10, 1.6049E-10, 1.7680E-10, 2.0557E-10, 2.4967E-10,    F15650
     *     3.0096E-10, 3.5816E-10, 4.0851E-10, 4.6111E-10, 5.2197E-10,    F15660
     *     5.5043E-10, 6.0324E-10, 6.4983E-10, 6.7498E-10, 7.0545E-10,    F15670
     *     7.0680E-10, 7.5218E-10, 7.5723E-10, 7.7840E-10, 8.0081E-10,    F15680
     *     8.0223E-10, 7.7271E-10, 7.1676E-10, 6.7819E-10, 6.4753E-10,    F15690
     *     6.5844E-10, 7.0163E-10, 7.7503E-10, 8.8152E-10, 9.9022E-10,    F15700
     *     1.0229E-09, 9.9296E-10, 8.9911E-10, 7.7813E-10, 6.3785E-10,    F15710
     *     4.7491E-10, 3.5280E-10, 2.4349E-10, 1.6502E-10, 1.1622E-10/    F15720
      DATA F1401/                                                         F15730
     *     8.6715E-11, 6.7360E-11, 5.3910E-11, 4.5554E-11, 4.1300E-11,    F15740
     *     3.9728E-11, 3.9000E-11, 3.9803E-11, 4.1514E-11, 4.3374E-11,    F15750
     *     4.6831E-11, 4.8921E-11, 5.1995E-11, 5.7242E-11, 6.2759E-11,    F15760
     *     7.0801E-11, 7.4555E-11, 7.9754E-11, 8.7616E-11, 9.1171E-11,    F15770
     *     1.0349E-10, 1.1047E-10, 1.2024E-10, 1.2990E-10, 1.3725E-10,    F15780
     *     1.5005E-10, 1.5268E-10, 1.5535E-10, 1.5623E-10, 1.5009E-10,    F15790
     *     1.4034E-10, 1.3002E-10, 1.2225E-10, 1.1989E-10, 1.2411E-10,    F15800
     *     1.3612E-10, 1.5225E-10, 1.7202E-10, 1.9471E-10, 1.9931E-10,    F15810
     *     1.9079E-10, 1.7478E-10, 1.5259E-10, 1.2625E-10, 9.3332E-11,    F15820
     *     6.8796E-11, 4.6466E-11, 2.9723E-11, 1.8508E-11, 1.2106E-11/    F15830
      DATA F1451/                                                         F15840
     *     8.0142E-12, 5.4066E-12, 3.9329E-12, 3.1665E-12, 2.7420E-12,    F15850
     *     2.3996E-12, 2.3804E-12, 2.3242E-12, 2.4476E-12, 2.5331E-12,    F15860
     *     2.3595E-12, 2.2575E-12, 2.1298E-12, 2.0088E-12, 1.8263E-12,    F15870
     *     1.6114E-12, 1.4422E-12, 1.2946E-12, 1.0837E-12, 9.1282E-13,    F15880
     *     7.2359E-13, 5.3307E-13, 3.8837E-13, 2.6678E-13, 1.6769E-13,    F15890
     *     1.0826E-13, 7.2364E-14, 4.5201E-14, 3.0808E-14, 2.2377E-14,    F15900
     *     1.7040E-14, 9.2181E-15, 5.2934E-15, 3.5774E-15, 3.1431E-15,    F15910
     *     3.7647E-15, 5.6428E-15, 9.5139E-15, 1.7322E-14, 2.8829E-14,    F15920
     *     4.7708E-14, 6.9789E-14, 9.7267E-14, 1.4662E-13, 1.9429E-13,    F15930
     *     2.5998E-13, 3.6636E-13, 4.7960E-13, 6.5129E-13, 7.7638E-13/    F15940
      DATA F1501/                                                         F15950
     *     9.3774E-13, 1.1467E-12, 1.3547E-12, 1.5686E-12, 1.6893E-12,    F15960
     *     1.9069E-12, 2.1352E-12, 2.3071E-12, 2.4759E-12, 2.8247E-12,    F15970
     *     3.4365E-12, 4.3181E-12, 5.6107E-12, 7.0017E-12, 8.6408E-12,    F15980
     *     1.0974E-11, 1.3742E-11, 1.6337E-11, 2.0157E-11, 2.3441E-11,    F15990
     *     2.6733E-11, 3.0247E-11, 3.3737E-11, 3.8618E-11, 4.1343E-11,    F16000
     *     4.3870E-11, 4.4685E-11, 4.4881E-11, 4.5526E-11, 4.3628E-11,    F16010
     *     4.4268E-11, 4.6865E-11, 5.3426E-11, 5.4020E-11, 5.3218E-11,    F16020
     *     5.4587E-11, 5.6360E-11, 5.7740E-11, 5.6426E-11, 6.0399E-11,    F16030
     *     6.6981E-11, 7.4319E-11, 7.7977E-11, 7.5539E-11, 7.1610E-11,    F16040
     *     6.4606E-11, 5.5498E-11, 4.3944E-11, 3.3769E-11, 2.5771E-11/    F16050
      DATA F1551/                                                         F16060
     *     1.9162E-11, 1.3698E-11, 1.0173E-11, 7.8925E-12, 6.1938E-12,    F16070
     *     4.7962E-12, 4.0811E-12, 3.3912E-12, 2.8625E-12, 2.4504E-12,    F16080
     *     2.2188E-12, 2.2139E-12, 2.2499E-12, 2.2766E-12, 2.3985E-12,    F16090
     *     2.5459E-12, 2.9295E-12, 3.4196E-12, 3.6155E-12, 4.0733E-12,    F16100
     *     4.4610E-12, 4.9372E-12, 5.4372E-12, 5.7304E-12, 6.1640E-12,    F16110
     *     6.1278E-12, 6.2940E-12, 6.4947E-12, 6.8174E-12, 7.5190E-12,    F16120
     *     8.2608E-12, 8.4971E-12, 8.3484E-12, 8.1888E-12, 7.8552E-12,    F16130
     *     7.8468E-12, 7.5943E-12, 7.9096E-12, 8.6869E-12, 9.1303E-12,    F16140
     *     9.2547E-12, 8.9322E-12, 8.2177E-12, 7.3408E-12, 5.7956E-12,    F16150
     *     4.4470E-12, 3.5881E-12, 2.6748E-12, 1.7074E-12, 9.6700E-13/    F16160
      DATA F1601/                                                         F16170
     *     5.2645E-13, 2.9943E-13, 1.7316E-13, 1.0039E-13, 5.7859E-14,    F16180
     *     3.2968E-14, 1.8499E-14, 1.0192E-14, 5.5015E-15, 2.9040E-15,    F16190
     *     1.4968E-15, 7.5244E-16, 3.6852E-16, 1.7568E-16, 8.1464E-17,    F16200
     *     3.6717E-17, 1.6076E-17, 6.8341E-18, 2.8195E-18, 1.1286E-18,    F16210
     *      .0000E+00,  .0000E+00,  .0000E+00,  .0000E+00,  .0000E+00,    F16220
     *      .0000E+00,  .0000E+00,  .0000E+00,  .0000E+00, 1.4070E-18,    F16230
     *     3.0405E-18, 6.4059E-18, 1.3169E-17, 2.6443E-17, 5.1917E-17,    F16240
     *     9.9785E-17, 1.8802E-16, 3.4788E-16, 6.3328E-16, 1.1370E-15,    F16250
     *     2.0198E-15, 3.5665E-15, 6.3053E-15, 1.1309E-14, 2.1206E-14,    F16260
     *     3.2858E-14, 5.5165E-14, 8.6231E-14, 1.2776E-13, 1.7780E-13/    F16270
      DATA F1651/                                                         F16280
     *     2.5266E-13, 3.6254E-13, 5.1398E-13, 6.8289E-13, 8.7481E-13,    F16290
     *     1.1914E-12, 1.6086E-12, 2.0469E-12, 2.5761E-12, 3.4964E-12,    F16300
     *     4.4980E-12, 5.5356E-12, 6.7963E-12, 8.5720E-12, 1.0700E-11,    F16310
     *     1.2983E-11, 1.6270E-11, 1.9609E-11, 2.2668E-11, 2.5963E-11,    F16320
     *     3.0918E-11, 3.4930E-11, 3.9330E-11, 4.4208E-11, 4.6431E-11,    F16330
     *     5.1141E-11, 5.4108E-11, 5.8077E-11, 6.5050E-11, 7.2126E-11,    F16340
     *     8.1064E-11, 8.1973E-11, 8.1694E-11, 8.3081E-11, 8.0240E-11,    F16350
     *     7.9225E-11, 7.6256E-11, 7.8468E-11, 8.0041E-11, 8.1585E-11,    F16360
     *     8.3485E-11, 8.3774E-11, 8.5870E-11, 8.6104E-11, 8.8516E-11,    F16370
     *     9.0814E-11, 9.2522E-11, 8.8913E-11, 7.8381E-11, 6.8568E-11/    F16380
      DATA F1701/                                                         F16390
     *     5.6797E-11, 4.4163E-11, 3.2369E-11, 2.3259E-11, 1.6835E-11,    F16400
     *     1.1733E-11, 8.5273E-12, 6.3805E-12, 4.8983E-12, 3.8831E-12,    F16410
     *     3.2610E-12, 2.8577E-12, 2.5210E-12, 2.2913E-12, 2.0341E-12,    F16420
     *     1.8167E-12, 1.6395E-12, 1.4890E-12, 1.3516E-12, 1.2542E-12,    F16430
     *     1.2910E-12, 1.3471E-12, 1.4689E-12, 1.5889E-12, 1.6989E-12,    F16440
     *     1.8843E-12, 2.0902E-12, 2.3874E-12, 2.7294E-12, 3.3353E-12,    F16450
     *     4.0186E-12, 4.5868E-12, 5.2212E-12, 5.8856E-12, 6.5991E-12,    F16460
     *     7.2505E-12, 7.6637E-12, 8.5113E-12, 9.4832E-12, 9.9678E-12,    F16470
     *     1.0723E-11, 1.0749E-11, 1.1380E-11, 1.1774E-11, 1.1743E-11,    F16480
     *     1.2493E-11, 1.2559E-11, 1.2332E-11, 1.1782E-11, 1.1086E-11/    F16490
      DATA F1751/                                                         F16500
     *     1.0945E-11, 1.1178E-11, 1.2083E-11, 1.3037E-11, 1.4730E-11,    F16510
     *     1.6450E-11, 1.7403E-11, 1.7004E-11, 1.5117E-11, 1.3339E-11,    F16520
     *     1.0844E-11, 8.0915E-12, 5.6615E-12, 3.7196E-12, 2.5194E-12,    F16530
     *     1.6569E-12, 1.1201E-12, 8.2335E-13, 6.0270E-13, 4.8205E-13,    F16540
     *     4.1313E-13, 3.6243E-13, 3.2575E-13, 2.7730E-13, 2.5292E-13,    F16550
     *     2.3062E-13, 2.1126E-13, 2.1556E-13, 2.1213E-13, 2.2103E-13,    F16560
     *     2.1927E-13, 2.0794E-13, 1.9533E-13, 1.6592E-13, 1.4521E-13,    F16570
     *     1.1393E-13, 8.3772E-14, 6.2077E-14, 4.3337E-14, 2.7165E-14,    F16580
     *     1.6821E-14, 9.5407E-15, 5.3093E-15, 3.0320E-15, 1.7429E-15,    F16590
     *     9.9828E-16, 5.6622E-16, 3.1672E-16, 1.7419E-16, 9.3985E-17/    F16600
      DATA F1801/                                                         F16610
     *     4.9656E-17, 2.5652E-17, 1.2942E-17, 6.3695E-18, 3.0554E-18,    F16620
     *     1.4273E-18, -0.       , -0.       , -0.       , -0.       ,    F16630
     *     -0.       , 0.        , 0.        , 0.        , 0.        ,    F16640
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16650
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16660
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16670
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16680
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16690
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16700
     *     0.        , 0.        , 0.        , 0.        , 0.        /    F16710
      DATA F1851/                                                         F16720
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16730
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16740
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16750
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16760
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16770
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16780
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16790
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16800
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16810
     *     0.        , 0.        , 0.        , 0.        , 0.        /    F16820
      DATA F1901/                                                         F16830
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16840
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16850
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16860
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16870
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16880
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16890
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16900
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16910
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16920
     *     0.        , 0.        , 0.        , 0.        , 0.        /    F16930
      DATA F1951/                                                         F16940
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16950
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16960
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16970
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16980
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F16990
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F17000
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F17010
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F17020
     *     0.        , 0.        , 0.        , 0.        , 0.        ,    F17030
     *     0.        , 0.        , 0.        , 0.        , 0.        /    F17040
      DATA F2001/                                                         F17050
     *     0.        /                                                    F17060
C                                                                         F17070
      END                                                                 F17080
C
C     --------------------------------------------------------------
C
      SUBROUTINE FRNCO2 (V1C,V2C,DVC,NPTC,C)                              F17090
C                                                                         F17100
      IMPLICIT REAL*8           (V)                                     ! F17110
C                                                                         F17120
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F17130
      COMMON /FCO2/ V1S,V2S,DVS,NPTS,S(1003)                              F17140
      DIMENSION C(*)                                                      F17150
C                                                                         F17160
      DVC = DVS                                                           F17170
      V1C = V1ABS-DVC                                                     F17180
      V2C = V2ABS+DVC                                                     F17190
C                                                                         F17200
      I1 = (V1C-V1S)/DVS                                                  F17210
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F17280
         I = I1+J                                                         F17290
         C(J) = 0.                                                        F17300
         IF ((I.GE.1).AND.(I.LE.NPTS)) THEN                               F17310
            C(J) = S(I)                                                   F17320
         ENDIF                                                            F17330
   10 CONTINUE                                                            F17340
C                                                                         F17350
      RETURN                                                              F17360
C                                                                         F17370
      END                                                                 F17380
C
C     --------------------------------------------------------------
C
      BLOCK DATA BFCO2                                                    F17390
C                                                                         F17400
      IMPLICIT REAL*8           (V)                                     ! F17410
C                                                                         F17420
C     CO2 CONTINUUM RIDGEWAY 1982                                         F17430
C               UNITS OF (CM**3/MOL)*1.E-20                               F17440
C                                                                         F17450
      COMMON /FCO2/ V1,V2,DV,NPT,                                         F17460
     *              C0000( 2),C0001(50),C0051(50),C0101(50),C0151(50),    F17470
     *              C0201(50),C0251(50),C0301(50),C0351(50),C0401(50),    F17480
     *              C0451(50),C0501(50),C0551(50),C0601(50),C0651(50),    F17490
     *              C0701(50),C0751(50),C0801(50),C0851(50),C0901(50),    F17500
     *              C0951(50),C1001(1)                                    F17510
C                                                                         F17520
       DATA V1,V2,DV,NPT / -20.0, 10000.0, 10.0, 1003/                    F17530
C                                                                         F17540
C                                                                         F17550
      DATA C0000/                                                         F17560
     *     1.1110E-11, 1.0188E-11/                                        F17570
      DATA C0001/                                                         F17580
     *     9.3516E-12, 1.0188E-11, 1.1110E-11, 1.2127E-11, 1.3251E-11,    F17590
     *     1.4495E-11, 1.5872E-11, 1.7400E-11, 1.9097E-11, 2.0985E-11,    F17600
     *     2.3087E-11, 2.5431E-11, 2.8051E-11, 3.0982E-11, 3.4268E-11,    F17610
     *     3.7956E-11, 4.2105E-11, 4.6779E-11, 5.2056E-11, 5.8025E-11,    F17620
     *     6.4791E-11, 7.2477E-11, 8.1226E-11, 9.1209E-11, 1.0263E-10,    F17630
     *     1.1572E-10, 1.3078E-10, 1.4814E-10, 1.6821E-10, 1.9148E-10,    F17640
     *     2.1857E-10, 2.5019E-10, 2.8723E-10, 3.3080E-10, 3.8223E-10,    F17650
     *     4.4321E-10, 5.1583E-10, 6.0274E-10, 7.0725E-10, 8.3363E-10,    F17660
     *     9.8735E-10, 1.1755E-09, 1.4074E-09, 1.6953E-09, 2.0557E-09,    F17670
     *     2.5107E-09, 3.0909E-09, 3.8391E-09, 4.8165E-09, 6.1117E-09/    F17680
      DATA C0051/                                                         F17690
     *     7.8550E-09, 1.0241E-08, 1.3593E-08, 1.8344E-08, 2.5408E-08,    F17700
     *     3.6386E-08, 5.4251E-08, 8.4262E-08, 1.3273E-07, 2.1867E-07,    F17710
     *     3.5007E-07, 6.0011E-07, 1.0797E-06, 1.8254E-06, 3.1621E-06,    F17720
     *     4.0293E-06, 4.3683E-06, 4.4552E-06, 4.2684E-06, 3.9341E-06,    F17730
     *     2.5972E-06, 1.5617E-06, 8.9063E-07, 5.0360E-07, 3.0616E-07,    F17740
     *     1.9066E-07, 1.1904E-07, 7.6078E-08, 4.9304E-08, 3.3335E-08,    F17750
     *     2.3494E-08, 1.7114E-08, 1.2742E-08, 9.6068E-09, 7.3706E-09,    F17760
     *     5.7386E-09, 4.5302E-09, 3.6223E-09, 2.9309E-09, 2.4001E-09,    F17770
     *     1.9927E-09, 1.6877E-09, 1.4602E-09, 1.2764E-09, 1.1317E-09,    F17780
     *     1.0273E-09, 9.1943E-10, 8.0353E-10, 6.8746E-10, 5.9354E-10/    F17790
      DATA C0101/                                                         F17800
     *     5.1722E-10, 4.4975E-10, 4.2350E-10, 4.2282E-10, 4.2610E-10,    F17810
     *     4.5465E-10, 4.6166E-10, 4.3149E-10, 3.7615E-10, 3.1576E-10,    F17820
     *     2.6490E-10, 1.9143E-10, 1.2885E-10, 9.4954E-11, 7.6499E-11,    F17830
     *     6.4581E-11, 5.5923E-11, 4.9200E-11, 4.3813E-11, 3.9533E-11,    F17840
     *     3.6338E-11, 3.4320E-11, 3.3329E-11, 3.2400E-11, 3.1700E-11,    F17850
     *     3.1267E-11, 2.9940E-11, 2.7628E-11, 2.4496E-11, 2.1764E-11,    F17860
     *     1.9306E-11, 1.7352E-11, 1.7292E-11, 1.8733E-11, 2.0224E-11,    F17870
     *     2.2396E-11, 2.4225E-11, 2.4890E-11, 2.3513E-11, 2.0824E-11,    F17880
     *     1.8642E-11, 1.5676E-11, 1.2882E-11, 1.1054E-11, 1.0074E-11,    F17890
     *     9.6324E-12, 9.4910E-12, 9.5134E-12, 9.6427E-12, 9.8552E-12/    F17900
      DATA C0151/                                                         F17910
     *     1.0140E-11, 1.0494E-11, 1.0915E-11, 1.1405E-11, 1.1965E-11,    F17920
     *     1.2601E-11, 1.3316E-11, 1.4116E-11, 1.5006E-11, 1.5997E-11,    F17930
     *     1.7092E-11, 1.8305E-11, 1.9641E-11, 2.1121E-11, 2.2744E-11,    F17940
     *     2.4503E-11, 2.6419E-11, 2.8221E-11, 3.0609E-11, 3.3260E-11,    F17950
     *     3.6247E-11, 3.9581E-11, 4.3279E-11, 4.7376E-11, 5.1932E-11,    F17960
     *     5.7001E-11, 6.2654E-11, 6.8973E-11, 7.6058E-11, 8.4037E-11,    F17970
     *     9.3081E-11, 1.0344E-10, 1.1547E-10, 1.2970E-10, 1.4659E-10,    F17980
     *     1.6724E-10, 1.9481E-10, 2.3520E-10, 2.9424E-10, 3.6319E-10,    F17990
     *     4.2279E-10, 4.8494E-10, 5.2296E-10, 5.6111E-10, 5.8935E-10,    F18000
     *     6.0807E-10, 6.4204E-10, 6.8457E-10, 7.6709E-10, 8.7664E-10/    F18010
      DATA C0201/                                                         F18020
     *     1.0183E-09, 1.2116E-09, 1.4874E-09, 1.8596E-09, 2.2742E-09,    F18030
     *     2.7577E-09, 3.1932E-09, 3.6381E-09, 4.1207E-09, 4.6458E-09,    F18040
     *     5.3065E-09, 6.0741E-09, 7.1942E-09, 8.7103E-09, 1.0713E-08,    F18050
     *     1.3344E-08, 1.6831E-08, 2.1524E-08, 2.7967E-08, 3.7047E-08,    F18060
     *     5.0312E-08, 7.0566E-08, 1.0275E-07, 1.5419E-07, 2.3309E-07,    F18070
     *     3.4843E-07, 5.3194E-07, 8.7207E-07, 1.5075E-06, 2.7077E-06,    F18080
     *     4.7125E-06, 7.1734E-06, 9.2381E-06, 1.1507E-05, 1.3737E-05,    F18090
     *     1.4004E-05, 1.2679E-05, 1.0478E-05, 8.5684E-06, 6.1472E-06,    F18100
     *     3.2424E-06, 1.5291E-06, 8.0390E-07, 4.6767E-07, 2.9170E-07,    F18110
     *     1.9148E-07, 1.3076E-07, 9.2156E-08, 6.6652E-08, 4.9265E-08/    F18120
      DATA C0251/                                                         F18130
     *     3.7094E-08, 2.8380E-08, 2.2019E-08, 1.7297E-08, 1.3738E-08,    F18140
     *     1.1019E-08, 8.9178E-09, 7.2762E-09, 5.9810E-09, 4.9500E-09,    F18150
     *     4.1226E-09, 3.4534E-09, 2.9082E-09, 2.4611E-09, 2.0922E-09,    F18160
     *     1.7864E-09, 1.5313E-09, 1.3176E-09, 1.1379E-09, 9.8612E-10,    F18170
     *     8.5741E-10, 7.4782E-10, 6.5416E-10, 5.7384E-10, 5.0471E-10,    F18180
     *     4.4503E-10, 3.9334E-10, 3.4841E-10, 3.0927E-10, 2.7510E-10,    F18190
     *     2.4519E-10, 2.1893E-10, 1.9587E-10, 1.7555E-10, 1.5762E-10,    F18200
     *     1.4178E-10, 1.2772E-10, 1.1524E-10, 1.0414E-10, 9.4248E-11,    F18210
     *     8.5421E-11, 7.7530E-11, 7.0466E-11, 6.4134E-11, 5.8450E-11,    F18220
     *     5.3342E-11, 4.8746E-11, 4.4607E-11, 4.0874E-11, 3.7507E-11/    F18230
      DATA C0301/                                                         F18240
     *     3.4466E-11, 3.1719E-11, 2.9237E-11, 2.6993E-11, 2.4968E-11,    F18250
     *     2.3139E-11, 2.1494E-11, 2.0022E-11, 1.8709E-11, 1.7541E-11,    F18260
     *     1.6533E-11, 1.5690E-11, 1.5027E-11, 1.4560E-11, 1.4169E-11,    F18270
     *     1.3796E-11, 1.3553E-11, 1.3526E-11, 1.3567E-11, 1.3399E-11,    F18280
     *     1.3149E-11, 1.3049E-11, 1.3078E-11, 1.3093E-11, 1.3168E-11,    F18290
     *     1.3572E-11, 1.4383E-11, 1.5698E-11, 1.7658E-11, 2.0197E-11,    F18300
     *     2.2845E-11, 2.5944E-11, 3.0250E-11, 3.5900E-11, 4.1482E-11,    F18310
     *     4.6602E-11, 5.2453E-11, 5.9754E-11, 6.9308E-11, 8.0696E-11,    F18320
     *     9.5737E-11, 1.1733E-10, 1.4793E-10, 1.9119E-10, 2.5355E-10,    F18330
     *     3.4588E-10, 4.8343E-10, 6.9378E-10, 1.0212E-09, 1.4858E-09/    F18340
      DATA C0351/                                                         F18350
     *     2.0906E-09, 3.0576E-09, 4.6318E-09, 7.1585E-09, 1.1259E-08,    F18360
     *     1.7954E-08, 2.9760E-08, 4.6693E-08, 6.2035E-08, 7.4399E-08,    F18370
     *     9.1705E-08, 9.9448E-08, 9.5181E-08, 8.3050E-08, 7.1756E-08,    F18380
     *     6.6261E-08, 6.0357E-08, 6.6988E-08, 8.3419E-08, 9.8834E-08,    F18390
     *     1.2385E-07, 1.3962E-07, 1.3651E-07, 1.1963E-07, 9.7731E-08,    F18400
     *     8.0083E-08, 5.1660E-08, 2.5778E-08, 1.2600E-08, 6.8779E-09,    F18410
     *     4.1161E-09, 2.6276E-09, 1.7595E-09, 1.2225E-09, 8.7493E-10,    F18420
     *     6.4179E-10, 4.7987E-10, 3.6491E-10, 2.8191E-10, 2.2084E-10,    F18430
     *     1.7507E-10, 1.4025E-10, 1.1344E-10, 9.2580E-11, 7.6170E-11,    F18440
     *     6.3142E-11, 5.2694E-11, 4.4260E-11, 3.7421E-11, 3.1847E-11/    F18450
      DATA C0401/                                                         F18460
     *     2.7263E-11, 2.3352E-11, 2.0081E-11, 1.7332E-11, 1.5000E-11,    F18470
     *     1.2978E-11, 1.1204E-11, 9.7513E-12, 8.5300E-12, 7.4888E-12,    F18480
     *     6.5947E-12, 5.8231E-12, 5.1548E-12, 4.5739E-12, 4.0675E-12,    F18490
     *     3.6250E-12, 3.2371E-12, 2.8963E-12, 2.5964E-12, 2.3316E-12,    F18500
     *     2.0975E-12, 1.8902E-12, 1.7061E-12, 1.5425E-12, 1.3967E-12,    F18510
     *     1.2665E-12, 1.1503E-12, 1.0463E-12, 9.5319E-13, 8.6963E-13,    F18520
     *     7.9461E-13, 7.2718E-13, 6.6654E-13, 6.1201E-13, 5.6296E-13,    F18530
     *     5.1894E-13, 4.7969E-13, 4.4494E-13, 4.1320E-13, 3.8529E-13,    F18540
     *     3.6202E-13, 3.4320E-13, 3.2546E-13, 3.0741E-13, 2.9156E-13,    F18550
     *     2.7819E-13, 2.6576E-13, 2.5327E-13, 2.4319E-13, 2.3770E-13/    F18560
      DATA C0451/                                                         F18570
     *     2.3645E-13, 2.3967E-13, 2.4960E-13, 2.6858E-13, 2.9679E-13,    F18580
     *     3.3247E-13, 3.8487E-13, 4.7576E-13, 6.1833E-13, 8.0740E-13,    F18590
     *     1.0267E-12, 1.2291E-12, 1.4710E-12, 1.7211E-12, 1.8251E-12,    F18600
     *     1.8982E-12, 1.9768E-12, 2.1877E-12, 2.5008E-12, 3.0545E-12,    F18610
     *     4.1513E-12, 5.7469E-12, 7.7913E-12, 1.0873E-11, 1.5538E-11,    F18620
     *     2.2838E-11, 3.4153E-11, 4.9751E-11, 7.0591E-11, 1.0794E-10,    F18630
     *     1.7287E-10, 2.6554E-10, 3.5250E-10, 4.1952E-10, 5.1979E-10,    F18640
     *     5.7649E-10, 5.6168E-10, 5.0014E-10, 4.3670E-10, 4.0057E-10,    F18650
     *     3.5169E-10, 3.7578E-10, 5.5054E-10, 8.8962E-10, 1.2940E-09,    F18660
     *     1.6293E-09, 2.0553E-09, 2.3945E-09, 2.3926E-09, 2.1385E-09/    F18670
      DATA C0501/                                                         F18680
     *     1.7637E-09, 1.4623E-09, 1.0150E-09, 5.5612E-10, 3.5162E-10,    F18690
     *     3.4009E-10, 4.1744E-10, 5.0009E-10, 6.0748E-10, 7.3258E-10,    F18700
     *     7.6553E-10, 7.2066E-10, 6.1317E-10, 5.1585E-10, 3.9136E-10,    F18710
     *     2.2991E-10, 1.2590E-10, 6.9549E-11, 3.8699E-11, 2.2976E-11,    F18720
     *     1.4702E-11, 9.9989E-12, 7.1233E-12, 5.2612E-12, 4.0298E-12,    F18730
     *     3.2395E-12, 2.7932E-12, 2.6331E-12, 2.7835E-12, 3.3167E-12,    F18740
     *     3.3581E-12, 3.3404E-12, 3.1243E-12, 2.8459E-12, 2.4092E-12,    F18750
     *     1.5349E-12, 9.7039E-13, 5.8611E-13, 3.9686E-13, 2.9332E-13,    F18760
     *     2.2795E-13, 1.8432E-13, 1.5287E-13, 1.2898E-13, 1.1019E-13,    F18770
     *     9.5041E-14, 8.2617E-14, 7.2310E-14, 6.3711E-14, 5.6561E-14/    F18780
      DATA C0551/                                                         F18790
     *     5.0763E-14, 4.6525E-14, 4.4418E-14, 4.4681E-14, 4.7199E-14,    F18800
     *     5.0389E-14, 5.3620E-14, 6.0817E-14, 6.0192E-14, 5.5878E-14,    F18810
     *     4.9874E-14, 4.3955E-14, 3.9854E-14, 3.1697E-14, 3.1135E-14,    F18820
     *     3.4683E-14, 3.8789E-14, 4.6932E-14, 5.0213E-14, 4.7156E-14,    F18830
     *     4.2130E-14, 3.5554E-14, 3.0465E-14, 1.9216E-14, 1.1378E-14,    F18840
     *     8.2878E-15, 6.8260E-15, 6.0960E-15, 5.8135E-15, 5.9618E-15,    F18850
     *     6.8295E-15, 9.2943E-15, 1.2572E-14, 1.4837E-14, 1.8595E-14,    F18860
     *     2.1533E-14, 2.2008E-14, 2.1305E-14, 1.9743E-14, 2.0413E-14,    F18870
     *     2.1131E-14, 2.5346E-14, 3.3709E-14, 4.3995E-14, 5.8911E-14,    F18880
     *     7.8451E-14, 1.0537E-13, 1.4559E-13, 2.0405E-13, 2.6734E-13/    F18890
      DATA C0601/                                                         F18900
     *     3.5029E-13, 4.9788E-13, 7.3207E-13, 1.0979E-12, 1.4960E-12,    F18910
     *     1.7906E-12, 2.2171E-12, 2.5369E-12, 2.5873E-12, 2.3871E-12,    F18920
     *     2.0730E-12, 1.9095E-12, 1.6227E-12, 1.3981E-12, 1.5228E-12,    F18930
     *     2.0956E-12, 3.2493E-12, 5.2740E-12, 8.6666E-12, 1.2672E-11,    F18940
     *     1.5725E-11, 1.9496E-11, 2.2858E-11, 2.2939E-11, 2.0597E-11,    F18950
     *     1.7021E-11, 1.4456E-11, 1.0794E-11, 7.1327E-12, 6.5438E-12,    F18960
     *     8.8057E-12, 1.2311E-11, 1.5284E-11, 1.9273E-11, 2.2796E-11,    F18970
     *     2.3156E-11, 2.0914E-11, 1.7298E-11, 1.4424E-11, 1.0127E-11,    F18980
     *     5.2952E-12, 2.5759E-12, 1.4304E-12, 9.4758E-13, 7.9895E-13,    F18990
     *     9.1124E-13, 1.2297E-12, 1.5898E-12, 1.9056E-12, 2.3905E-12/    F19000
      DATA C0651/                                                         F19010
     *     2.6695E-12, 2.6297E-12, 2.3467E-12, 2.0058E-12, 1.6773E-12,    F19020
     *     1.1327E-12, 6.7331E-13, 4.0954E-13, 2.5152E-13, 1.4491E-13,    F19030
     *     9.0916E-14, 6.6510E-14, 5.9022E-14, 6.4403E-14, 8.3126E-14,    F19040
     *     1.2409E-13, 1.5153E-13, 1.6909E-13, 1.7938E-13, 1.9169E-13,    F19050
     *     2.1173E-13, 2.1941E-13, 2.6360E-13, 3.5956E-13, 4.8369E-13,    F19060
     *     5.9657E-13, 7.4062E-13, 8.9452E-13, 8.7899E-13, 8.2012E-13,    F19070
     *     7.4109E-13, 6.9845E-13, 6.3130E-13, 5.6538E-13, 6.9516E-13,    F19080
     *     9.9486E-13, 1.5226E-12, 2.4155E-12, 3.9119E-12, 6.3541E-12,    F19090
     *     1.0075E-11, 1.5903E-11, 2.5091E-11, 3.6282E-11, 4.6076E-11,    F19100
     *     5.6240E-11, 7.1126E-11, 7.0230E-11, 6.3642E-11, 5.3722E-11/    F19110
      DATA C0701/                                                         F19120
     *     4.4651E-11, 3.4409E-11, 1.5287E-11, 7.2479E-12, 3.9218E-12,    F19130
     *     2.3172E-12, 1.4585E-12, 9.6297E-13, 6.6017E-13, 4.6655E-13,    F19140
     *     3.3814E-13, 2.5034E-13, 1.8874E-13, 1.4457E-13, 1.1228E-13,    F19150
     *     8.8284E-14, 7.0188E-14, 5.6365E-14, 4.5685E-14, 3.7357E-14,    F19160
     *     3.0817E-14, 2.5674E-14, 2.1679E-14, 1.8780E-14, 1.7243E-14,    F19170
     *     1.6273E-14, 1.5201E-14, 1.5091E-14, 1.4725E-14, 1.3668E-14,    F19180
     *     1.1940E-14, 1.0097E-14, 8.8905E-15, 7.1475E-15, 5.8080E-15,    F19190
     *     5.5216E-15, 5.9338E-15, 7.1932E-15, 9.9780E-15, 1.6167E-14,    F19200
     *     2.9100E-14, 5.2355E-14, 8.4889E-14, 1.1311E-13, 1.4192E-13,    F19210
     *     1.7648E-13, 1.8657E-13, 1.7498E-13, 1.4877E-13, 1.2578E-13/    F19220
      DATA C0751/                                                         F19230
     *     1.0051E-13, 6.7213E-14, 5.4750E-14, 7.0454E-14, 1.1351E-13,    F19240
     *     1.8015E-13, 2.4825E-13, 3.0875E-13, 3.9200E-13, 4.2550E-13,    F19250
     *     4.0067E-13, 3.4438E-13, 2.8204E-13, 2.2432E-13, 1.3172E-13,    F19260
     *     6.2820E-14, 3.6474E-14, 2.9409E-14, 3.4164E-14, 4.8300E-14,    F19270
     *     6.4140E-14, 7.7284E-14, 9.7973E-14, 1.0969E-13, 1.0580E-13,    F19280
     *     9.2070E-14, 7.5008E-14, 6.1722E-14, 3.8874E-14, 1.9007E-14,    F19290
     *     9.6765E-15, 5.5169E-15, 3.5254E-15, 2.5012E-15, 2.0013E-15,    F19300
     *     1.8810E-15, 2.2143E-15, 3.5332E-15, 5.7552E-15, 7.3359E-15,    F19310
     *     8.3292E-15, 9.9174E-15, 1.0930E-14, 1.1185E-14, 1.0884E-14,    F19320
     *     1.0577E-14, 1.1048E-14, 1.1611E-14, 1.1128E-14, 1.0729E-14/    F19330
      DATA C0801/                                                         F19340
     *     1.0248E-14, 1.0630E-14, 1.1793E-14, 1.3977E-14, 1.9857E-14,    F19350
     *     2.9182E-14, 4.2229E-14, 6.2710E-14, 9.0717E-14, 1.2561E-13,    F19360
     *     1.6951E-13, 2.2520E-13, 3.2470E-13, 4.5178E-13, 6.3104E-13,    F19370
     *     8.7521E-13, 1.1073E-12, 1.3534E-12, 1.6954E-12, 1.7005E-12,    F19380
     *     1.5993E-12, 1.4416E-12, 1.3280E-12, 1.2760E-12, 1.1076E-12,    F19390
     *     1.2850E-12, 1.6208E-12, 1.9527E-12, 2.4941E-12, 2.5077E-12,    F19400
     *     2.3156E-12, 2.0069E-12, 1.6301E-12, 1.2885E-12, 5.9863E-13,    F19410
     *     2.8012E-13, 1.5065E-13, 8.8802E-14, 5.5888E-14, 3.6951E-14,    F19420
     *     2.5393E-14, 1.8001E-14, 1.3093E-14, 9.7308E-15, 7.3665E-15,    F19430
     *     5.6662E-15, 4.4194E-15, 3.4897E-15, 2.7857E-15, 2.2457E-15/    F19440
      DATA C0851/                                                         F19450
     *     1.8264E-15, 1.4973E-15, 1.2365E-15, 1.0280E-15, 8.5996E-16,    F19460
     *     7.2345E-16, 6.1182E-16, 5.1994E-16, 4.4388E-16, 3.8055E-16,    F19470
     *     3.2756E-16, 2.8300E-16, 2.4537E-16, 2.1347E-16, 1.8630E-16,    F19480
     *     1.6307E-16, 1.4314E-16, 1.2599E-16, 1.1117E-16, 9.8344E-17,    F19490
     *     8.7197E-17, 7.7487E-17, 6.9004E-17, 6.1577E-17, 5.5060E-17,    F19500
     *     4.9325E-17, 4.4271E-17, 3.9810E-17, 3.5861E-17, 3.2361E-17,    F19510
     *     2.9252E-17, 2.6487E-17, 2.4023E-17, 2.1826E-17, 1.9862E-17,    F19520
     *     1.8107E-17, 1.6536E-17, 1.5129E-17, 1.3869E-17, 1.2739E-17,    F19530
     *     1.1726E-17, 1.0820E-17, 1.0009E-17, 9.2846E-18, 8.6398E-18,    F19540
     *     8.0682E-18, 7.5641E-18, 7.1229E-18, 6.7411E-18, 6.4161E-18/    F19550
      DATA C0901/                                                         F19560
     *     6.1455E-18, 5.9290E-18, 5.7662E-18, 5.6574E-18, 5.6049E-18,    F19570
     *     5.6112E-18, 5.6811E-18, 5.8200E-18, 6.0364E-18, 6.3405E-18,    F19580
     *     6.7450E-18, 7.2674E-18, 7.9298E-18, 8.7613E-18, 9.8010E-18,    F19590
     *     1.1086E-17, 1.2686E-17, 1.4679E-17, 1.7177E-17, 2.0335E-17,    F19600
     *     2.4384E-17, 2.9538E-17, 3.6416E-17, 4.5520E-17, 5.7788E-17,    F19610
     *     7.4676E-17, 9.8513E-17, 1.3323E-16, 1.8570E-16, 2.6897E-16,    F19620
     *     4.0958E-16, 6.6785E-16, 1.2064E-15, 2.4023E-15, 4.3240E-15,    F19630
     *     6.6353E-15, 8.6393E-15, 1.1433E-14, 1.3946E-14, 1.3611E-14,    F19640
     *     1.2557E-14, 1.0934E-14, 1.0039E-14, 8.5099E-15, 7.9557E-15,    F19650
     *     1.1346E-14, 1.8512E-14, 2.9285E-14, 4.1585E-14, 5.2809E-14/    F19660
      DATA C0951/                                                         F19670
     *     7.0377E-14, 7.8094E-14, 7.3735E-14, 6.5845E-14, 5.5023E-14,    F19680
     *     4.6866E-14, 2.7430E-14, 1.5975E-14, 1.4522E-14, 1.7075E-14,    F19690
     *     2.0408E-14, 2.5119E-14, 3.1194E-14, 3.0280E-14, 2.7676E-14,    F19700
     *     2.3344E-14, 1.9466E-14, 1.4140E-14, 6.2087E-15, 3.0307E-15,    F19710
     *     1.6815E-15, 1.0169E-15, 6.5448E-16, 4.4162E-16, 3.0928E-16,    F19720
     *     2.2320E-16, 1.6511E-16, 1.2471E-16, 9.5881E-17, 7.4850E-17,    F19730
     *     5.9216E-17, 4.7400E-17, 3.8338E-17, 3.1298E-17, 2.5765E-17,    F19740
     *     2.1371E-17, 1.7848E-17, 1.5000E-17, 1.2679E-17, 1.0774E-17,    F19750
     *     9.2002E-18, 7.8922E-18, 6.7987E-18, 5.8800E-18, 5.1042E-18,    F19760
     *     4.4461E-18, 3.8855E-18, 3.4060E-18, 2.9944E-18, 2.6397E-18/    F19770
      DATA C1001/                                                         F19780
     *     2.3331E-18/                                                    F19790
C                                                                         F19800
      END                                                                 F19810
C
C     --------------------------------------------------------------
C
      SUBROUTINE N2R296 (V1C,V2C,DVC,NPTC,C)
C
C     Model used:
C      Borysow, A, and L. Frommhold, "Collision-induced
C         rototranslational absorption spectra of N2-N2
C         pairs for temperatures from 50 to 300 K", The
C         Astrophysical Journal, 311, 1043-1057, 1986.
C
      IMPLICIT REAL*8           (V)
C
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)
      COMMON /N2RT0/ V1S,V2S,DVS,NPTS,S(73)
      DIMENSION C(*)
C
      DVC = DVS
      V1C = V1ABS-DVC
      V2C = V2ABS+DVC
C
      I1 = (V1C-V1S)/DVS
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
c
      DO 10 J = 1, NPTC
         I = I1+J
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
   10 CONTINUE
C
      RETURN
C
      END
C
      BLOCK DATA BN2T0
C
      IMPLICIT REAL*8           (V)
C
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
C
C           THESE DATA ARE FOR 296K
C
      COMMON /N2RT0/ V1N2CR,V2N2CR,DVN2CR,NPTN2C,CT296(73)
C
      DATA V1N2CR,V2N2CR,DVN2CR,NPTN2C / -10., 350., 5.0, 73 /
C
      DATA CT296/
     *     0.4303E-06, 0.4850E-06, 0.4979E-06, 0.4850E-06, 0.4303E-06,
     *     0.3715E-06, 0.3292E-06, 0.3086E-06, 0.2920E-06, 0.2813E-06,
     *     0.2804E-06, 0.2738E-06, 0.2726E-06, 0.2724E-06, 0.2635E-06,
     *     0.2621E-06, 0.2547E-06, 0.2428E-06, 0.2371E-06, 0.2228E-06,
     *     0.2100E-06, 0.1991E-06, 0.1822E-06, 0.1697E-06, 0.1555E-06,
     *     0.1398E-06, 0.1281E-06, 0.1138E-06, 0.1012E-06, 0.9078E-07,
     *     0.7879E-07, 0.6944E-07, 0.6084E-07, 0.5207E-07, 0.4540E-07,
     *     0.3897E-07, 0.3313E-07, 0.2852E-07, 0.2413E-07, 0.2045E-07,
     *     0.1737E-07, 0.1458E-07, 0.1231E-07, 0.1031E-07, 0.8586E-08,
     *     0.7162E-08, 0.5963E-08, 0.4999E-08, 0.4226E-08, 0.3607E-08,
     *     0.3090E-08, 0.2669E-08, 0.2325E-08, 0.2024E-08, 0.1783E-08,
     *     0.1574E-08, 0.1387E-08, 0.1236E-08, 0.1098E-08, 0.9777E-09,
     *     0.8765E-09, 0.7833E-09, 0.7022E-09, 0.6317E-09, 0.5650E-09,
     *     0.5100E-09, 0.4572E-09, 0.4115E-09, 0.3721E-09, 0.3339E-09,
     *     0.3005E-09, 0.2715E-09, 0.2428E-09/
C
      END
C
      SUBROUTINE N2R220 (V1C,V2C,DVC,NPTC,C)
C
C     Model used:
C      Borysow, A, and L. Frommhold, "Collision-induced
C         rototranslational absorption spectra of N2-N2
C         pairs for temperatures from 50 to 300 K", The
C         Astrophysical Journal, 311, 1043-1057, 1986.
C
      IMPLICIT REAL*8           (V)
C
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)
      COMMON /N2RT1/ V1S,V2S,DVS,NPTS,S(73)
      DIMENSION C(*)
C
      DVC = DVS
      V1C = V1ABS-DVC
      V2C = V2ABS+DVC
C
      I1 = (V1C-V1S)/DVS
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
c
      DO 10 J = 1, NPTC
         I = I1+J
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
   10 CONTINUE
C
      RETURN
C
      END
C
      BLOCK DATA BN2T1
C
      IMPLICIT REAL*8           (V)
C
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
C
C         THESE DATA ARE FOR 220K
C
      COMMON /N2RT1/ V1N2CR,V2N2CR,DVN2CR,NPTN2C,CT220(73)
C
      DATA V1N2CR,V2N2CR,DVN2CR,NPTN2C / -10., 350., 5.0, 73 /
C
      DATA CT220/
     *     0.4946E-06, 0.5756E-06, 0.5964E-06, 0.5756E-06, 0.4946E-06,
     *     0.4145E-06, 0.3641E-06, 0.3482E-06, 0.3340E-06, 0.3252E-06,
     *     0.3299E-06, 0.3206E-06, 0.3184E-06, 0.3167E-06, 0.2994E-06,
     *     0.2943E-06, 0.2794E-06, 0.2582E-06, 0.2468E-06, 0.2237E-06,
     *     0.2038E-06, 0.1873E-06, 0.1641E-06, 0.1474E-06, 0.1297E-06,
     *     0.1114E-06, 0.9813E-07, 0.8309E-07, 0.7059E-07, 0.6068E-07,
     *     0.5008E-07, 0.4221E-07, 0.3537E-07, 0.2885E-07, 0.2407E-07,
     *     0.1977E-07, 0.1605E-07, 0.1313E-07, 0.1057E-07, 0.8482E-08,
     *     0.6844E-08, 0.5595E-08, 0.4616E-08, 0.3854E-08, 0.3257E-08,
     *     0.2757E-08, 0.2372E-08, 0.2039E-08, 0.1767E-08, 0.1548E-08,
     *     0.1346E-08, 0.1181E-08, 0.1043E-08, 0.9110E-09, 0.8103E-09,
     *     0.7189E-09, 0.6314E-09, 0.5635E-09, 0.4976E-09, 0.4401E-09,
     *     0.3926E-09, 0.3477E-09, 0.3085E-09, 0.2745E-09, 0.2416E-09,
     *     0.2155E-09, 0.1895E-09, 0.1678E-09, 0.1493E-09, 0.1310E-09,
     *     0.1154E-09, 0.1019E-09, 0.8855E-10/
 
C
      END
C
C     --------------------------------------------------------------
C
      subroutine n2_ver_1 (v1c,v2c,dvc,nptc,c,T)
c
      IMPLICIT REAL*8 (v)                                                    
c
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)
c
      COMMON /n2_f/ V1S,V2S,DVS,NPTS,xn2(118),xn2t(118)
c
      dimension c(*)
c
c     Nitrogen Collision Induced Fundamental

c     Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M. Hartmann,
c        Infrared collision-induced absorption by N2 near 4.3 microns for
c        atmospheric applications: measurements and emprirical modeling, 
c         Appl. Optics, 35, 5911-5917, (1996).
c
      DATA  To/ 296./, xlosmt/ 2.68675e+19/, vmr_n2/ 0.78 /
c
      xktfac = (1./To)-(1./T)
c     
      a1 = 0.8387
      a2 = 0.0754
c
c     correct formulation for consistency with LBLRTM:
c
      factor = (1.e+20 /xlosmt) * (1./vmr_n2) * (a1-a2*(T/To))
c
c     Lafferty et al. reference  assumes that the
c     column amount is that for air 

C                           
      DVC = DVS             
      V1C = V1ABS-DVC       
      V2C = V2ABS+DVC       
C                           
      I1 = (V1C-V1S)/DVS    
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      do 10 j=1,nptc
         i = i1+j
         C(J) = 0.                                                        F41620
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F41630
         VJ = V1C+DVC* REAL(J-1)                                          F41640
c     the radiation field is removed with 1/vj
c
         c(j) = factor * xn2(i)* exp(xn2t(i)*xktfac) / vj
c
 10   end do
 920  format (f10.2,1p,e12.2,0p,f10.2,1p2e12.2)
      return

      end

      BLOCK DATA bn2f                                                   
                                                                        
      IMPLICIT REAL*8 (v)                                                    
                                                                        
      COMMON /n2_f/ V1n2f,V2n2f,DVn2f,NPTn2f,                                      
     *          xn0001(50),xn0051(50),xn0101(18),
     *          xnt0001(50),xnt0051(50),xnt0101(18)
                                                                        
      DATA V1n2f,V2n2f,DVn2f,NPTn2f /2085.000,2670.000, 5.000, 118/              
      DATA xn0001/                                                       
     *      0.000E+00,  2.000E-10,  5.200E-09,  1.020E-08,  1.520E-08,  
     *      2.020E-08,  2.520E-08,  3.020E-08,  4.450E-08,  5.220E-08,  
     *      6.460E-08,  7.750E-08,  9.030E-08,  1.060E-07,  1.210E-07,  
     *      1.370E-07,  1.570E-07,  1.750E-07,  2.010E-07,  2.300E-07,  
     *      2.590E-07,  2.950E-07,  3.260E-07,  3.660E-07,  4.050E-07,  
     *      4.470E-07,  4.920E-07,  5.340E-07,  5.840E-07,  6.240E-07,  
     *      6.670E-07,  7.140E-07,  7.260E-07,  7.540E-07,  7.840E-07,  
     *      8.090E-07,  8.420E-07,  8.620E-07,  8.870E-07,  9.110E-07,  
     *      9.360E-07,  9.760E-07,  1.030E-06,  1.110E-06,  1.230E-06,  
     *      1.390E-06,  1.610E-06,  1.760E-06,  1.940E-06,  1.970E-06/  
      DATA xn0051/                                                       
     *      1.870E-06,  1.750E-06,  1.560E-06,  1.420E-06,  1.350E-06,  
     *      1.320E-06,  1.290E-06,  1.290E-06,  1.290E-06,  1.300E-06,  
     *      1.320E-06,  1.330E-06,  1.340E-06,  1.350E-06,  1.330E-06,  
     *      1.310E-06,  1.290E-06,  1.240E-06,  1.200E-06,  1.160E-06,  
     *      1.100E-06,  1.040E-06,  9.960E-07,  9.380E-07,  8.630E-07,  
     *      7.980E-07,  7.260E-07,  6.550E-07,  5.940E-07,  5.350E-07,  
     *      4.740E-07,  4.240E-07,  3.770E-07,  3.330E-07,  2.960E-07,  
     *      2.630E-07,  2.340E-07,  2.080E-07,  1.850E-07,  1.670E-07,  
     *      1.470E-07,  1.320E-07,  1.200E-07,  1.090E-07,  9.850E-08,  
     *      9.080E-08,  8.180E-08,  7.560E-08,  6.850E-08,  6.140E-08/  
      DATA xn0101/                                                       
     *      5.830E-08,  5.770E-08,  5.000E-08,  4.320E-08,  3.140E-08,  
     *      2.890E-08,  2.640E-08,  2.390E-08,  2.140E-08,  1.890E-08,  
     *      1.640E-08,  1.390E-08,  1.140E-08,  8.900E-09,  6.400E-09,  
     *      3.900E-09,  1.400E-09,  0.000E+00/                          
                                                                        
c     temperature coefficients:

      DATA xnt0001/                                                       
     *      1.040E+03,  1.010E+03,  9.800E+02,  9.500E+02,  9.200E+02,  
     *      8.900E+02,  8.600E+02,  8.300E+02,  8.020E+02,  7.610E+02,  
     *      7.220E+02,  6.790E+02,  6.460E+02,  6.090E+02,  5.620E+02,  
     *      5.110E+02,  4.720E+02,  4.360E+02,  4.060E+02,  3.770E+02,  
     *      3.550E+02,  3.380E+02,  3.190E+02,  2.990E+02,  2.780E+02,  
     *      2.550E+02,  2.330E+02,  2.080E+02,  1.840E+02,  1.490E+02,  
     *      1.070E+02,  6.600E+01,  2.500E+01, -1.300E+01, -4.900E+01,  
     *     -8.200E+01, -1.040E+02, -1.190E+02, -1.300E+02, -1.390E+02,  
     *     -1.440E+02, -1.460E+02, -1.460E+02, -1.470E+02, -1.480E+02,  
     *     -1.500E+02, -1.530E+02, -1.600E+02, -1.690E+02, -1.810E+02/  
      DATA xnt0051/                                                       
     *     -1.890E+02, -1.950E+02, -2.000E+02, -2.050E+02, -2.090E+02,  
     *     -2.110E+02, -2.100E+02, -2.100E+02, -2.090E+02, -2.050E+02,  
     *     -1.990E+02, -1.900E+02, -1.800E+02, -1.680E+02, -1.570E+02,  
     *     -1.430E+02, -1.260E+02, -1.080E+02, -8.900E+01, -6.300E+01,  
     *     -3.200E+01,  1.000E+00,  3.500E+01,  6.500E+01,  9.500E+01,  
     *      1.210E+02,  1.410E+02,  1.520E+02,  1.610E+02,  1.640E+02,  
     *      1.640E+02,  1.610E+02,  1.550E+02,  1.480E+02,  1.430E+02,  
     *      1.370E+02,  1.330E+02,  1.310E+02,  1.330E+02,  1.390E+02,  
     *      1.500E+02,  1.650E+02,  1.870E+02,  2.130E+02,  2.480E+02,  
     *      2.840E+02,  3.210E+02,  3.720E+02,  4.490E+02,  5.140E+02/  
      DATA xnt0101/                                                       
     *      5.690E+02,  6.090E+02,  6.420E+02,  6.730E+02,  7.000E+02,  
     *      7.300E+02,  7.600E+02,  7.900E+02,  8.200E+02,  8.500E+02,  
     *      8.800E+02,  9.100E+02,  9.400E+02,  9.700E+02,  1.000E+03,  
     *      1.030E+03,  1.060E+03,  1.090E+03/                          

      END
C
C     --------------------------------------------------------------
C
      SUBROUTINE XO3CHP (V1C,V2C,DVC,NPTC,C0,C1,C2)                       F20520
C                                                                         F20530
      IMPLICIT REAL*8           (V)                                     ! F20540
C                                                                         F20550
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F20560
      COMMON /O3CHAP/ V1S,V2S,DVS,NPTS,X(3150),Y(3150),Z(3150)            F20570
      DIMENSION C0(*),C1(*),C2(*)                                         F20580
C                                                                         F20590
      DVC = DVS                                                           F20600
      V1C = V1ABS-DVC                                                     F20610
      V2C = V2ABS+DVC                                                     F20620
C                                                                         F20630
      I1 = (V1C-V1S)/DVS                                                  F20640
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F20710
         I = I1+J                                                         F20720
         IF ((I.LT.1).OR.(I.GT.NPTS)) THEN
             C0(J) = 0.
             C1(J)=0.
             C2(J)=0.
         ELSE
C
C            Remove radiation field from diffuse ozone
C
             VJ = V1C+DVC* REAL(J-1)
             C0(J)=X(I)/VJ
             C1(J)=Y(I)/VJ
             C2(J)=Z(I)/VJ
         ENDIF
   10 CONTINUE                                                            F20800
C                                                                         F20810
      RETURN                                                              F20820
C                                                                         F20830
      END                                                                 F20840
C
C     --------------------------------------------------------------
      BLOCK DATA O3CH
C
C     CHAPPUIS AND WULF BAND
C
C     BEGINNING AND ENDING FREQUENCIES FROM DATA (CM-1):
C
C                        9170.0 24565.0
C
C     Added points at beginning and end (X,Y,Z(1:50) and
C     X,Y,Z(3130:3150)).  Zeroed values of Y,Z(1:789) to eliminate
C     ringing from interpolations done in MODTRAN.  Changed coefficients
C     X(32:50,3130:3150) Y(821:841,3130:3150), & Z(821:841,3130:3150) to
C     smooth coefficients to zero.
C     Smoothing coefficient frequencies (cm-1):
C
C             9075.0 -  9165.0  and 24570.0 - 24665.0 for X
C            13020.0 - 13120.0  and 24570.0 - 24665.0 for Y
C            13020.0 - 13120.0  and 24570.0 - 24665.0 for Z
C             
C     
C
C     CROSS-SECTIONS IN CM^2 TIMES 1.0E20
C     FORMULA FOR CROSS SECTION:  X+Y*DT+Z*DT*DT, DT=T-273.15
C     THE OUTPUT OF THIS ROUTINE IS C0=X, CT1=Y AND CT2=Z.
C
      IMPLICIT REAL*8           (V)
      COMMON /O3CHAP/VBEG,VEND,DVINCR,NMAX,X(3150),Y(3150),Z(3150)
C
      DATA VBEG, VEND, DVINCR, NMAX /8920.0, 24665.0, 5.0, 3150/
      DATA (X(I),I=    1,  50)/
     1      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,
     2      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,
     3      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,
     4      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,
     5      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,
     6      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,
     7      0.00000  ,  0.075E-05,  0.150E-05,  0.225E-05,  0.300E-05,
     8      0.400E-05,  0.500E-05,  0.600E-05,  0.700E-05,  0.850E-05,
     9      1.000E-05,  1.200E-05,  1.430E-05,  1.680E-05,  1.980E-05,
     *      2.280E-05,  2.630E-05,  2.980E-05,  3.376E-05,  3.826E-05/
      DATA (X(I),I= 51, 100)/
     1      4.276E-05,  4.775E-05,  5.825E-05,  6.908E-05,  7.299E-05,
     2      7.116E-05,  7.388E-05,  7.965E-05,  7.689E-05,  6.900E-05,
     3      7.008E-05,  6.945E-05,  7.083E-05,  7.053E-05,  6.908E-05,
     4      6.923E-05,  6.770E-05,  7.146E-05,  7.749E-05,  8.464E-05,
     5      8.441E-05,  8.754E-05,  8.795E-05,  9.971E-05,  9.632E-05,
     6      9.539E-05,  1.037E-04,  1.085E-04,  1.058E-04,  1.077E-04,
     7      1.121E-04,  1.193E-04,  1.292E-04,  1.364E-04,  1.526E-04,
     8      1.658E-04,  1.808E-04,  1.861E-04,  1.786E-04,  1.804E-04,
     9      1.885E-04,  1.972E-04,  2.218E-04,  2.408E-04,  2.317E-04,
     $      2.098E-04,  1.938E-04,  1.851E-04,  1.896E-04,  1.875E-04/
      DATA (X(I),I=  101,  150)/
     1      1.708E-04,  1.710E-04,  1.796E-04,  1.865E-04,  1.943E-04,
     2      1.881E-04,  1.885E-04,  2.136E-04,  2.255E-04,  2.267E-04,
     3      2.234E-04,  2.418E-04,  2.695E-04,  2.710E-04,  2.738E-04,
     4      3.066E-04,  3.269E-04,  3.465E-04,  3.986E-04,  4.410E-04,
     5      4.719E-04,  5.051E-04,  5.211E-04,  5.132E-04,  5.125E-04,
     6      5.159E-04,  5.549E-04,  6.562E-04,  7.168E-04,  6.502E-04,
     7      5.140E-04,  4.161E-04,  3.620E-04,  3.264E-04,  3.004E-04,
     8      2.815E-04,  2.650E-04,  2.527E-04,  2.424E-04,  2.292E-04,
     9      2.155E-04,  2.072E-04,  1.992E-04,  1.943E-04,  1.914E-04,
     $      1.855E-04,  1.813E-04,  1.724E-04,  1.687E-04,  1.676E-04/
      DATA (X(I),I= 151, 200)/
     1      1.601E-04,  1.503E-04,  1.518E-04,  1.436E-04,  1.455E-04,
     2      1.448E-04,  1.410E-04,  1.406E-04,  1.425E-04,  1.407E-04,
     3      1.405E-04,  1.436E-04,  1.369E-04,  1.355E-04,  1.331E-04,
     4      1.328E-04,  1.350E-04,  1.394E-04,  1.372E-04,  1.444E-04,
     5      1.490E-04,  1.455E-04,  1.460E-04,  1.523E-04,  1.559E-04,
     6      1.654E-04,  1.766E-04,  1.843E-04,  1.911E-04,  1.881E-04,
     7      1.894E-04,  1.927E-04,  2.043E-04,  2.106E-04,  2.215E-04,
     8      2.268E-04,  2.249E-04,  2.230E-04,  2.302E-04,  2.408E-04,
     9      2.518E-04,  2.625E-04,  2.753E-04,  2.788E-04,  2.701E-04,
     $      2.746E-04,  2.935E-04,  3.173E-04,  3.457E-04,  3.452E-04/
      DATA (X(I),I=  201,  250)/
     1      3.329E-04,  3.443E-04,  3.706E-04,  4.079E-04,  4.403E-04,
     2      4.343E-04,  4.172E-04,  4.448E-04,  5.132E-04,  5.635E-04,
     3      5.590E-04,  5.419E-04,  6.007E-04,  6.912E-04,  7.258E-04,
     4      7.146E-04,  7.529E-04,  8.706E-04,  9.465E-04,  9.923E-04,
     5      1.134E-03,  1.286E-03,  1.351E-03,  1.485E-03,  1.709E-03,
     6      1.897E-03,  2.086E-03,  2.186E-03,  2.195E-03,  2.185E-03,
     7      2.199E-03,  2.336E-03,  2.666E-03,  3.076E-03,  3.075E-03,
     8      2.543E-03,  1.920E-03,  1.498E-03,  1.283E-03,  1.165E-03,
     9      1.070E-03,  9.833E-04,  9.018E-04,  8.207E-04,  7.451E-04,
     $      6.811E-04,  6.178E-04,  5.661E-04,  5.199E-04,  4.868E-04/
      DATA (X(I),I= 251, 300)/
     1      4.541E-04,  4.291E-04,  4.135E-04,  3.990E-04,  3.878E-04,
     2      3.815E-04,  3.722E-04,  3.691E-04,  3.726E-04,  3.711E-04,
     3      3.744E-04,  3.778E-04,  3.808E-04,  3.826E-04,  3.852E-04,
     4      3.919E-04,  3.975E-04,  3.990E-04,  4.053E-04,  4.176E-04,
     5      4.232E-04,  4.291E-04,  4.414E-04,  4.541E-04,  4.723E-04,
     6      4.887E-04,  5.058E-04,  5.226E-04,  5.501E-04,  5.836E-04,
     7      6.059E-04,  6.238E-04,  6.469E-04,  6.711E-04,  7.046E-04,
     8      7.448E-04,  7.794E-04,  8.054E-04,  8.222E-04,  8.371E-04,
     9      8.538E-04,  8.612E-04,  8.698E-04,  8.914E-04,  9.122E-04,
     $      9.305E-04,  9.562E-04,  9.844E-04,  1.018E-03,  1.053E-03/
      DATA (X(I),I=  301,  350)/
     1      1.091E-03,  1.136E-03,  1.187E-03,  1.233E-03,  1.289E-03,
     2      1.336E-03,  1.372E-03,  1.405E-03,  1.435E-03,  1.470E-03,
     3      1.504E-03,  1.517E-03,  1.511E-03,  1.541E-03,  1.619E-03,
     4      1.728E-03,  1.848E-03,  1.955E-03,  2.044E-03,  2.128E-03,
     5      2.254E-03,  2.396E-03,  2.527E-03,  2.660E-03,  2.832E-03,
     6      3.010E-03,  3.182E-03,  3.340E-03,  3.504E-03,  3.673E-03,
     7      3.822E-03,  3.923E-03,  3.997E-03,  4.042E-03,  4.061E-03,
     8      4.035E-03,  3.979E-03,  3.901E-03,  3.785E-03,  3.642E-03,
     9      3.494E-03,  3.339E-03,  3.173E-03,  3.004E-03,  2.849E-03,
     $      2.703E-03,  2.556E-03,  2.432E-03,  2.310E-03,  2.191E-03/
      DATA (X(I),I= 351, 400)/
     1      2.076E-03,  1.969E-03,  1.883E-03,  1.818E-03,  1.753E-03,
     2      1.705E-03,  1.672E-03,  1.643E-03,  1.617E-03,  1.616E-03,
     3      1.629E-03,  1.648E-03,  1.662E-03,  1.667E-03,  1.669E-03,
     4      1.664E-03,  1.655E-03,  1.645E-03,  1.643E-03,  1.642E-03,
     5      1.632E-03,  1.629E-03,  1.632E-03,  1.638E-03,  1.644E-03,
     6      1.647E-03,  1.646E-03,  1.642E-03,  1.638E-03,  1.632E-03,
     7      1.628E-03,  1.626E-03,  1.628E-03,  1.635E-03,  1.642E-03,
     8      1.649E-03,  1.653E-03,  1.656E-03,  1.660E-03,  1.669E-03,
     9      1.685E-03,  1.705E-03,  1.730E-03,  1.755E-03,  1.779E-03,
     $      1.804E-03,  1.830E-03,  1.861E-03,  1.896E-03,  1.931E-03/
      DATA (X(I),I=  401,  450)/
     1      1.962E-03,  1.991E-03,  2.024E-03,  2.068E-03,  2.131E-03,
     2      2.207E-03,  2.285E-03,  2.357E-03,  2.423E-03,  2.490E-03,
     3      2.564E-03,  2.649E-03,  2.743E-03,  2.842E-03,  2.943E-03,
     4      3.044E-03,  3.146E-03,  3.248E-03,  3.350E-03,  3.452E-03,
     5      3.555E-03,  3.664E-03,  3.785E-03,  3.927E-03,  4.083E-03,
     6      4.250E-03,  4.418E-03,  4.570E-03,  4.708E-03,  4.835E-03,
     7      4.961E-03,  5.088E-03,  5.218E-03,  5.348E-03,  5.471E-03,
     8      5.594E-03,  5.713E-03,  5.828E-03,  5.933E-03,  6.026E-03,
     9      6.100E-03,  6.152E-03,  6.186E-03,  6.193E-03,  6.182E-03,
     $      6.149E-03,  6.093E-03,  6.011E-03,  5.914E-03,  5.799E-03/
      DATA (X(I),I= 451, 500)/
     1      5.676E-03,  5.553E-03,  5.438E-03,  5.330E-03,  5.233E-03,
     2      5.151E-03,  5.080E-03,  5.025E-03,  4.987E-03,  4.972E-03,
     3      4.976E-03,  4.991E-03,  5.013E-03,  5.032E-03,  5.043E-03,
     4      5.043E-03,  5.032E-03,  5.010E-03,  4.980E-03,  4.950E-03,
     5      4.913E-03,  4.879E-03,  4.838E-03,  4.786E-03,  4.723E-03,
     6      4.652E-03,  4.578E-03,  4.503E-03,  4.433E-03,  4.366E-03,
     7      4.306E-03,  4.247E-03,  4.191E-03,  4.135E-03,  4.083E-03,
     8      4.035E-03,  3.997E-03,  3.968E-03,  3.945E-03,  3.923E-03,
     9      3.904E-03,  3.886E-03,  3.867E-03,  3.856E-03,  3.848E-03,
     $      3.845E-03,  3.848E-03,  3.860E-03,  3.878E-03,  3.897E-03/
      DATA (X(I),I=  501,  550)/
     1      3.915E-03,  3.941E-03,  3.971E-03,  4.008E-03,  4.057E-03,
     2      4.113E-03,  4.176E-03,  4.243E-03,  4.325E-03,  4.418E-03,
     3      4.518E-03,  4.626E-03,  4.723E-03,  4.812E-03,  4.891E-03,
     4      4.972E-03,  5.058E-03,  5.155E-03,  5.263E-03,  5.382E-03,
     5      5.516E-03,  5.657E-03,  5.810E-03,  5.974E-03,  6.145E-03,
     6      6.331E-03,  6.532E-03,  6.752E-03,  6.993E-03,  7.247E-03,
     7      7.507E-03,  7.768E-03,  8.036E-03,  8.304E-03,  8.579E-03,
     8      8.862E-03,  9.148E-03,  9.442E-03,  9.744E-03,  1.006E-02,
     9      1.038E-02,  1.071E-02,  1.104E-02,  1.137E-02,  1.168E-02,
     $      1.195E-02,  1.220E-02,  1.242E-02,  1.264E-02,  1.283E-02/
      DATA (X(I),I= 551, 600)/
     1      1.303E-02,  1.322E-02,  1.339E-02,  1.356E-02,  1.371E-02,
     2      1.385E-02,  1.398E-02,  1.408E-02,  1.415E-02,  1.417E-02,
     3      1.415E-02,  1.408E-02,  1.395E-02,  1.376E-02,  1.353E-02,
     4      1.326E-02,  1.295E-02,  1.262E-02,  1.228E-02,  1.194E-02,
     5      1.161E-02,  1.128E-02,  1.097E-02,  1.067E-02,  1.038E-02,
     6      1.011E-02,  9.859E-03,  9.625E-03,  9.409E-03,  9.208E-03,
     7      9.022E-03,  8.843E-03,  8.668E-03,  8.505E-03,  8.348E-03,
     8      8.207E-03,  8.088E-03,  7.987E-03,  7.909E-03,  7.842E-03,
     9      7.782E-03,  7.727E-03,  7.675E-03,  7.619E-03,  7.570E-03,
     $      7.526E-03,  7.488E-03,  7.459E-03,  7.440E-03,  7.429E-03/
      DATA (X(I),I=  601,  650)/
     1      7.429E-03,  7.429E-03,  7.440E-03,  7.455E-03,  7.474E-03,
     2      7.500E-03,  7.529E-03,  7.563E-03,  7.593E-03,  7.622E-03,
     3      7.649E-03,  7.675E-03,  7.715E-03,  7.771E-03,  7.846E-03,
     4      7.939E-03,  8.039E-03,  8.147E-03,  8.255E-03,  8.367E-03,
     5      8.482E-03,  8.605E-03,  8.746E-03,  8.903E-03,  9.078E-03,
     6      9.271E-03,  9.472E-03,  9.677E-03,  9.889E-03,  1.011E-02,
     7      1.034E-02,  1.059E-02,  1.085E-02,  1.113E-02,  1.143E-02,
     8      1.174E-02,  1.207E-02,  1.242E-02,  1.277E-02,  1.313E-02,
     9      1.350E-02,  1.388E-02,  1.425E-02,  1.464E-02,  1.503E-02,
     $      1.544E-02,  1.586E-02,  1.628E-02,  1.670E-02,  1.713E-02/
      DATA (X(I),I= 651, 700)/
     1      1.755E-02,  1.796E-02,  1.837E-02,  1.875E-02,  1.911E-02,
     2      1.945E-02,  1.975E-02,  2.002E-02,  2.028E-02,  2.050E-02,
     3      2.070E-02,  2.089E-02,  2.104E-02,  2.117E-02,  2.126E-02,
     4      2.132E-02,  2.135E-02,  2.135E-02,  2.130E-02,  2.123E-02,
     5      2.114E-02,  2.101E-02,  2.087E-02,  2.072E-02,  2.053E-02,
     6      2.032E-02,  2.010E-02,  1.986E-02,  1.963E-02,  1.939E-02,
     7      1.915E-02,  1.891E-02,  1.868E-02,  1.845E-02,  1.821E-02,
     8      1.798E-02,  1.773E-02,  1.746E-02,  1.719E-02,  1.692E-02,
     9      1.666E-02,  1.643E-02,  1.621E-02,  1.598E-02,  1.576E-02,
     $      1.558E-02,  1.542E-02,  1.529E-02,  1.519E-02,  1.509E-02/
      DATA (X(I),I=  701,  750)/
     1      1.501E-02,  1.493E-02,  1.484E-02,  1.477E-02,  1.473E-02,
     2      1.471E-02,  1.469E-02,  1.468E-02,  1.468E-02,  1.470E-02,
     3      1.473E-02,  1.475E-02,  1.476E-02,  1.477E-02,  1.480E-02,
     4      1.484E-02,  1.489E-02,  1.497E-02,  1.507E-02,  1.520E-02,
     5      1.533E-02,  1.543E-02,  1.551E-02,  1.557E-02,  1.563E-02,
     6      1.569E-02,  1.575E-02,  1.581E-02,  1.591E-02,  1.602E-02,
     7      1.614E-02,  1.625E-02,  1.637E-02,  1.654E-02,  1.676E-02,
     8      1.698E-02,  1.719E-02,  1.740E-02,  1.762E-02,  1.784E-02,
     9      1.805E-02,  1.827E-02,  1.852E-02,  1.878E-02,  1.906E-02,
     $      1.935E-02,  1.962E-02,  1.985E-02,  2.005E-02,  2.024E-02/
      DATA (X(I),I= 751, 800)/
     1      2.047E-02,  2.073E-02,  2.106E-02,  2.137E-02,  2.167E-02,
     2      2.194E-02,  2.223E-02,  2.256E-02,  2.287E-02,  2.316E-02,
     3      2.344E-02,  2.373E-02,  2.405E-02,  2.445E-02,  2.490E-02,
     4      2.542E-02,  2.592E-02,  2.640E-02,  2.684E-02,  2.729E-02,
     5      2.778E-02,  2.828E-02,  2.876E-02,  2.914E-02,  2.948E-02,
     6      2.979E-02,  3.010E-02,  3.039E-02,  3.066E-02,  3.091E-02,
     7      3.116E-02,  3.138E-02,  3.153E-02,  3.155E-02,  3.152E-02,
     8      3.146E-02,  3.144E-02,  3.138E-02,  3.126E-02,  3.110E-02,
     9      3.092E-02,  3.073E-02,  3.054E-02,  3.033E-02,  3.008E-02,
     $      2.980E-02,  2.947E-02,  2.910E-02,  2.870E-02,  2.832E-02/
      DATA (X(I),I=  801,  850)/
     1      2.795E-02,  2.765E-02,  2.735E-02,  2.706E-02,  2.680E-02,
     2      2.656E-02,  2.637E-02,  2.620E-02,  2.604E-02,  2.587E-02,
     3      2.570E-02,  2.551E-02,  2.533E-02,  2.520E-02,  2.514E-02,
     4      2.513E-02,  2.513E-02,  2.512E-02,  2.509E-02,  2.506E-02,
     5      2.504E-02,  2.501E-02,  2.498E-02,  2.494E-02,  2.493E-02,
     6      2.497E-02,  2.508E-02,  2.522E-02,  2.535E-02,  2.545E-02,
     7      2.550E-02,  2.558E-02,  2.568E-02,  2.578E-02,  2.587E-02,
     8      2.592E-02,  2.598E-02,  2.605E-02,  2.619E-02,  2.631E-02,
     9      2.621E-02,  2.617E-02,  2.629E-02,  2.642E-02,  2.654E-02,
     $      2.669E-02,  2.685E-02,  2.700E-02,  2.716E-02,  2.734E-02/
      DATA (X(I),I= 851, 900)/
     1      2.752E-02,  2.772E-02,  2.792E-02,  2.813E-02,  2.834E-02,
     2      2.858E-02,  2.885E-02,  2.913E-02,  2.941E-02,  2.973E-02,
     3      3.005E-02,  3.038E-02,  3.075E-02,  3.117E-02,  3.159E-02,
     4      3.202E-02,  3.246E-02,  3.290E-02,  3.335E-02,  3.384E-02,
     5      3.438E-02,  3.493E-02,  3.547E-02,  3.603E-02,  3.660E-02,
     6      3.718E-02,  3.772E-02,  3.826E-02,  3.879E-02,  3.931E-02,
     7      3.987E-02,  4.042E-02,  4.098E-02,  4.151E-02,  4.199E-02,
     8      4.243E-02,  4.287E-02,  4.316E-02,  4.344E-02,  4.369E-02,
     9      4.392E-02,  4.405E-02,  4.417E-02,  4.429E-02,  4.436E-02,
     $      4.436E-02,  4.438E-02,  4.437E-02,  4.427E-02,  4.416E-02/
      DATA (X(I),I=  901,  950)/
     1      4.405E-02,  4.394E-02,  4.383E-02,  4.372E-02,  4.359E-02,
     2      4.344E-02,  4.329E-02,  4.312E-02,  4.299E-02,  4.289E-02,
     3      4.278E-02,  4.269E-02,  4.258E-02,  4.242E-02,  4.227E-02,
     4      4.213E-02,  4.202E-02,  4.194E-02,  4.183E-02,  4.178E-02,
     5      4.179E-02,  4.177E-02,  4.175E-02,  4.174E-02,  4.174E-02,
     6      4.175E-02,  4.177E-02,  4.183E-02,  4.191E-02,  4.199E-02,
     7      4.207E-02,  4.214E-02,  4.219E-02,  4.226E-02,  4.232E-02,
     8      4.239E-02,  4.247E-02,  4.254E-02,  4.261E-02,  4.270E-02,
     9      4.279E-02,  4.287E-02,  4.298E-02,  4.311E-02,  4.323E-02,
     $      4.338E-02,  4.357E-02,  4.375E-02,  4.393E-02,  4.413E-02/
      DATA (X(I),I= 951, 1000)/
     1      4.433E-02,  4.454E-02,  4.475E-02,  4.499E-02,  4.525E-02,
     2      4.551E-02,  4.579E-02,  4.613E-02,  4.645E-02,  4.678E-02,
     3      4.712E-02,  4.749E-02,  4.785E-02,  4.820E-02,  4.856E-02,
     4      4.894E-02,  4.928E-02,  4.963E-02,  4.991E-02,  5.016E-02,
     5      5.042E-02,  5.072E-02,  5.109E-02,  5.144E-02,  5.183E-02,
     6      5.218E-02,  5.251E-02,  5.282E-02,  5.315E-02,  5.351E-02,
     7      5.391E-02,  5.430E-02,  5.471E-02,  5.510E-02,  5.548E-02,
     8      5.588E-02,  5.628E-02,  5.673E-02,  5.722E-02,  5.771E-02,
     9      5.821E-02,  5.874E-02,  5.927E-02,  5.980E-02,  6.034E-02,
     $      6.088E-02,  6.144E-02,  6.197E-02,  6.250E-02,  6.303E-02/
      DATA (X(I),I= 1001, 1050)/
     1      6.352E-02,  6.404E-02,  6.452E-02,  6.493E-02,  6.537E-02,
     2      6.578E-02,  6.617E-02,  6.653E-02,  6.688E-02,  6.722E-02,
     3      6.747E-02,  6.768E-02,  6.788E-02,  6.808E-02,  6.827E-02,
     4      6.842E-02,  6.859E-02,  6.875E-02,  6.884E-02,  6.889E-02,
     5      6.896E-02,  6.900E-02,  6.915E-02,  6.927E-02,  6.942E-02,
     6      6.956E-02,  6.972E-02,  6.990E-02,  7.009E-02,  7.024E-02,
     7      7.043E-02,  7.062E-02,  7.081E-02,  7.102E-02,  7.127E-02,
     8      7.151E-02,  7.175E-02,  7.199E-02,  7.225E-02,  7.248E-02,
     9      7.274E-02,  7.300E-02,  7.325E-02,  7.351E-02,  7.375E-02,
     $      7.402E-02,  7.429E-02,  7.458E-02,  7.488E-02,  7.514E-02/
      DATA (X(I),I= 1051, 1100)/
     1      7.546E-02,  7.575E-02,  7.605E-02,  7.634E-02,  7.667E-02,
     2      7.698E-02,  7.732E-02,  7.767E-02,  7.803E-02,  7.841E-02,
     3      7.879E-02,  7.916E-02,  7.959E-02,  8.001E-02,  8.042E-02,
     4      8.082E-02,  8.118E-02,  8.152E-02,  8.188E-02,  8.224E-02,
     5      8.270E-02,  8.318E-02,  8.367E-02,  8.415E-02,  8.467E-02,
     6      8.518E-02,  8.569E-02,  8.621E-02,  8.673E-02,  8.725E-02,
     7      8.779E-02,  8.831E-02,  8.887E-02,  8.945E-02,  9.003E-02,
     8      9.060E-02,  9.123E-02,  9.187E-02,  9.254E-02,  9.317E-02,
     9      9.382E-02,  9.444E-02,  9.506E-02,  9.570E-02,  9.634E-02,
     $      9.702E-02,  9.769E-02,  9.838E-02,  9.904E-02,  9.968E-02/
      DATA (X(I),I= 1101, 1150)/
     1      1.003E-01,  1.010E-01,  1.016E-01,  1.022E-01,  1.028E-01,
     2      1.033E-01,  1.039E-01,  1.046E-01,  1.053E-01,  1.060E-01,
     3      1.067E-01,  1.075E-01,  1.082E-01,  1.089E-01,  1.096E-01,
     4      1.103E-01,  1.110E-01,  1.117E-01,  1.125E-01,  1.132E-01,
     5      1.139E-01,  1.147E-01,  1.154E-01,  1.162E-01,  1.169E-01,
     6      1.177E-01,  1.184E-01,  1.191E-01,  1.197E-01,  1.203E-01,
     7      1.209E-01,  1.215E-01,  1.221E-01,  1.227E-01,  1.232E-01,
     8      1.238E-01,  1.244E-01,  1.249E-01,  1.254E-01,  1.259E-01,
     9      1.263E-01,  1.268E-01,  1.273E-01,  1.278E-01,  1.283E-01,
     $      1.288E-01,  1.292E-01,  1.296E-01,  1.300E-01,  1.305E-01/
      DATA (X(I),I= 1151, 1200)/
     1      1.309E-01,  1.314E-01,  1.319E-01,  1.324E-01,  1.329E-01,
     2      1.335E-01,  1.340E-01,  1.345E-01,  1.351E-01,  1.356E-01,
     3      1.362E-01,  1.368E-01,  1.374E-01,  1.380E-01,  1.386E-01,
     4      1.392E-01,  1.399E-01,  1.406E-01,  1.413E-01,  1.420E-01,
     5      1.427E-01,  1.434E-01,  1.442E-01,  1.449E-01,  1.457E-01,
     6      1.465E-01,  1.472E-01,  1.479E-01,  1.487E-01,  1.495E-01,
     7      1.502E-01,  1.509E-01,  1.516E-01,  1.523E-01,  1.530E-01,
     8      1.539E-01,  1.547E-01,  1.555E-01,  1.563E-01,  1.571E-01,
     9      1.580E-01,  1.588E-01,  1.596E-01,  1.605E-01,  1.614E-01,
     $      1.623E-01,  1.632E-01,  1.641E-01,  1.649E-01,  1.658E-01/
      DATA (X(I),I= 1201, 1250)/
     1      1.666E-01,  1.675E-01,  1.684E-01,  1.692E-01,  1.701E-01,
     2      1.710E-01,  1.719E-01,  1.728E-01,  1.737E-01,  1.746E-01,
     3      1.756E-01,  1.764E-01,  1.774E-01,  1.783E-01,  1.792E-01,
     4      1.801E-01,  1.810E-01,  1.820E-01,  1.829E-01,  1.838E-01,
     5      1.848E-01,  1.857E-01,  1.866E-01,  1.876E-01,  1.885E-01,
     6      1.893E-01,  1.902E-01,  1.911E-01,  1.920E-01,  1.928E-01,
     7      1.936E-01,  1.945E-01,  1.953E-01,  1.961E-01,  1.969E-01,
     8      1.978E-01,  1.986E-01,  1.994E-01,  2.002E-01,  2.010E-01,
     9      2.018E-01,  2.026E-01,  2.034E-01,  2.041E-01,  2.049E-01,
     $      2.057E-01,  2.065E-01,  2.073E-01,  2.081E-01,  2.089E-01/
      DATA (X(I),I= 1251, 1300)/
     1      2.097E-01,  2.105E-01,  2.113E-01,  2.121E-01,  2.129E-01,
     2      2.137E-01,  2.146E-01,  2.154E-01,  2.163E-01,  2.172E-01,
     3      2.180E-01,  2.190E-01,  2.198E-01,  2.207E-01,  2.216E-01,
     4      2.225E-01,  2.234E-01,  2.243E-01,  2.251E-01,  2.260E-01,
     5      2.269E-01,  2.277E-01,  2.285E-01,  2.294E-01,  2.302E-01,
     6      2.311E-01,  2.320E-01,  2.328E-01,  2.337E-01,  2.346E-01,
     7      2.355E-01,  2.364E-01,  2.372E-01,  2.381E-01,  2.390E-01,
     8      2.398E-01,  2.407E-01,  2.416E-01,  2.424E-01,  2.432E-01,
     9      2.440E-01,  2.448E-01,  2.456E-01,  2.464E-01,  2.473E-01,
     $      2.482E-01,  2.491E-01,  2.500E-01,  2.509E-01,  2.517E-01/
      DATA (X(I),I= 1301, 1350)/
     1      2.525E-01,  2.533E-01,  2.541E-01,  2.550E-01,  2.559E-01,
     2      2.568E-01,  2.577E-01,  2.587E-01,  2.597E-01,  2.607E-01,
     3      2.617E-01,  2.626E-01,  2.636E-01,  2.645E-01,  2.654E-01,
     4      2.663E-01,  2.672E-01,  2.682E-01,  2.692E-01,  2.703E-01,
     5      2.713E-01,  2.724E-01,  2.734E-01,  2.744E-01,  2.754E-01,
     6      2.764E-01,  2.774E-01,  2.784E-01,  2.795E-01,  2.806E-01,
     7      2.816E-01,  2.827E-01,  2.838E-01,  2.850E-01,  2.861E-01,
     8      2.872E-01,  2.884E-01,  2.895E-01,  2.907E-01,  2.918E-01,
     9      2.930E-01,  2.942E-01,  2.954E-01,  2.967E-01,  2.980E-01,
     $      2.993E-01,  3.005E-01,  3.017E-01,  3.029E-01,  3.041E-01/
      DATA (X(I),I= 1351, 1400)/
     1      3.052E-01,  3.064E-01,  3.076E-01,  3.088E-01,  3.100E-01,
     2      3.112E-01,  3.124E-01,  3.136E-01,  3.149E-01,  3.161E-01,
     3      3.173E-01,  3.185E-01,  3.196E-01,  3.208E-01,  3.219E-01,
     4      3.230E-01,  3.242E-01,  3.253E-01,  3.265E-01,  3.277E-01,
     5      3.289E-01,  3.300E-01,  3.312E-01,  3.323E-01,  3.334E-01,
     6      3.345E-01,  3.356E-01,  3.367E-01,  3.378E-01,  3.389E-01,
     7      3.400E-01,  3.410E-01,  3.421E-01,  3.431E-01,  3.441E-01,
     8      3.452E-01,  3.462E-01,  3.472E-01,  3.482E-01,  3.493E-01,
     9      3.503E-01,  3.513E-01,  3.524E-01,  3.534E-01,  3.545E-01,
     $      3.555E-01,  3.566E-01,  3.577E-01,  3.588E-01,  3.599E-01/
      DATA (X(I),I= 1401, 1450)/
     1      3.611E-01,  3.622E-01,  3.632E-01,  3.643E-01,  3.653E-01,
     2      3.664E-01,  3.674E-01,  3.683E-01,  3.692E-01,  3.702E-01,
     3      3.711E-01,  3.721E-01,  3.732E-01,  3.740E-01,  3.751E-01,
     4      3.762E-01,  3.774E-01,  3.781E-01,  3.792E-01,  3.800E-01,
     5      3.811E-01,  3.818E-01,  3.826E-01,  3.837E-01,  3.844E-01,
     6      3.853E-01,  3.863E-01,  3.870E-01,  3.881E-01,  3.889E-01,
     7      3.898E-01,  3.908E-01,  3.916E-01,  3.928E-01,  3.937E-01,
     8      3.946E-01,  3.957E-01,  3.967E-01,  3.976E-01,  3.983E-01,
     9      3.995E-01,  4.002E-01,  4.013E-01,  4.023E-01,  4.032E-01,
     $      4.042E-01,  4.050E-01,  4.062E-01,  4.073E-01,  4.084E-01/
      DATA (X(I),I= 1451, 1500)/
     1      4.095E-01,  4.106E-01,  4.117E-01,  4.130E-01,  4.144E-01,
     2      4.155E-01,  4.165E-01,  4.178E-01,  4.189E-01,  4.200E-01,
     3      4.215E-01,  4.226E-01,  4.241E-01,  4.252E-01,  4.266E-01,
     4      4.280E-01,  4.293E-01,  4.306E-01,  4.321E-01,  4.335E-01,
     5      4.347E-01,  4.362E-01,  4.376E-01,  4.388E-01,  4.403E-01,
     6      4.418E-01,  4.433E-01,  4.450E-01,  4.465E-01,  4.480E-01,
     7      4.495E-01,  4.510E-01,  4.524E-01,  4.540E-01,  4.555E-01,
     8      4.571E-01,  4.585E-01,  4.603E-01,  4.619E-01,  4.634E-01,
     9      4.652E-01,  4.668E-01,  4.683E-01,  4.700E-01,  4.716E-01,
     $      4.734E-01,  4.750E-01,  4.764E-01,  4.782E-01,  4.797E-01/
      DATA (X(I),I= 1501, 1550)/
     1      4.812E-01,  4.830E-01,  4.845E-01,  4.861E-01,  4.875E-01,
     2      4.893E-01,  4.908E-01,  4.920E-01,  4.935E-01,  4.950E-01,
     3      4.964E-01,  4.976E-01,  4.990E-01,  5.003E-01,  5.016E-01,
     4      5.028E-01,  5.043E-01,  5.054E-01,  5.064E-01,  5.073E-01,
     5      5.082E-01,  5.094E-01,  5.104E-01,  5.111E-01,  5.119E-01,
     6      5.126E-01,  5.133E-01,  5.141E-01,  5.146E-01,  5.149E-01,
     7      5.154E-01,  5.159E-01,  5.162E-01,  5.166E-01,  5.167E-01,
     8      5.168E-01,  5.170E-01,  5.171E-01,  5.171E-01,  5.169E-01,
     9      5.166E-01,  5.162E-01,  5.159E-01,  5.155E-01,  5.151E-01,
     $      5.146E-01,  5.139E-01,  5.133E-01,  5.128E-01,  5.120E-01/
      DATA (X(I),I= 1551, 1600)/
     1      5.113E-01,  5.103E-01,  5.092E-01,  5.080E-01,  5.070E-01,
     2      5.059E-01,  5.048E-01,  5.036E-01,  5.022E-01,  5.011E-01,
     3      4.997E-01,  4.984E-01,  4.969E-01,  4.954E-01,  4.939E-01,
     4      4.924E-01,  4.909E-01,  4.893E-01,  4.877E-01,  4.863E-01,
     5      4.845E-01,  4.831E-01,  4.814E-01,  4.798E-01,  4.781E-01,
     6      4.766E-01,  4.748E-01,  4.732E-01,  4.718E-01,  4.703E-01,
     7      4.688E-01,  4.673E-01,  4.658E-01,  4.643E-01,  4.630E-01,
     8      4.615E-01,  4.600E-01,  4.586E-01,  4.572E-01,  4.559E-01,
     9      4.548E-01,  4.536E-01,  4.524E-01,  4.512E-01,  4.501E-01,
     $      4.491E-01,  4.483E-01,  4.475E-01,  4.468E-01,  4.459E-01/
      DATA (X(I),I= 1601, 1650)/
     1      4.450E-01,  4.444E-01,  4.438E-01,  4.431E-01,  4.424E-01,
     2      4.416E-01,  4.412E-01,  4.409E-01,  4.405E-01,  4.401E-01,
     3      4.397E-01,  4.394E-01,  4.392E-01,  4.390E-01,  4.389E-01,
     4      4.386E-01,  4.386E-01,  4.384E-01,  4.384E-01,  4.385E-01,
     5      4.385E-01,  4.385E-01,  4.385E-01,  4.387E-01,  4.387E-01,
     6      4.387E-01,  4.387E-01,  4.387E-01,  4.387E-01,  4.390E-01,
     7      4.391E-01,  4.394E-01,  4.398E-01,  4.398E-01,  4.402E-01,
     8      4.406E-01,  4.410E-01,  4.413E-01,  4.417E-01,  4.421E-01,
     9      4.425E-01,  4.428E-01,  4.432E-01,  4.440E-01,  4.443E-01,
     $      4.448E-01,  4.452E-01,  4.459E-01,  4.467E-01,  4.471E-01/
      DATA (X(I),I= 1651, 1700)/
     1      4.479E-01,  4.486E-01,  4.491E-01,  4.498E-01,  4.505E-01,
     2      4.512E-01,  4.519E-01,  4.525E-01,  4.532E-01,  4.539E-01,
     3      4.547E-01,  4.554E-01,  4.562E-01,  4.569E-01,  4.577E-01,
     4      4.584E-01,  4.592E-01,  4.599E-01,  4.606E-01,  4.614E-01,
     5      4.621E-01,  4.629E-01,  4.636E-01,  4.640E-01,  4.648E-01,
     6      4.655E-01,  4.662E-01,  4.670E-01,  4.675E-01,  4.682E-01,
     7      4.689E-01,  4.697E-01,  4.701E-01,  4.708E-01,  4.712E-01,
     8      4.718E-01,  4.724E-01,  4.729E-01,  4.735E-01,  4.739E-01,
     9      4.742E-01,  4.745E-01,  4.748E-01,  4.751E-01,  4.753E-01,
     $      4.755E-01,  4.757E-01,  4.757E-01,  4.757E-01,  4.756E-01/
      DATA (X(I),I= 1701, 1750)/
     1      4.756E-01,  4.756E-01,  4.753E-01,  4.752E-01,  4.749E-01,
     2      4.747E-01,  4.744E-01,  4.741E-01,  4.737E-01,  4.734E-01,
     3      4.730E-01,  4.725E-01,  4.721E-01,  4.715E-01,  4.708E-01,
     4      4.701E-01,  4.693E-01,  4.686E-01,  4.681E-01,  4.673E-01,
     5      4.663E-01,  4.657E-01,  4.649E-01,  4.641E-01,  4.632E-01,
     6      4.623E-01,  4.615E-01,  4.606E-01,  4.596E-01,  4.588E-01,
     7      4.579E-01,  4.569E-01,  4.561E-01,  4.551E-01,  4.542E-01,
     8      4.532E-01,  4.524E-01,  4.513E-01,  4.506E-01,  4.498E-01,
     9      4.487E-01,  4.479E-01,  4.472E-01,  4.461E-01,  4.454E-01,
     $      4.443E-01,  4.435E-01,  4.428E-01,  4.418E-01,  4.411E-01/
      DATA (X(I),I= 1751, 1800)/
     1      4.400E-01,  4.388E-01,  4.380E-01,  4.368E-01,  4.357E-01,
     2      4.347E-01,  4.338E-01,  4.328E-01,  4.316E-01,  4.305E-01,
     3      4.294E-01,  4.283E-01,  4.272E-01,  4.261E-01,  4.249E-01,
     4      4.235E-01,  4.222E-01,  4.212E-01,  4.201E-01,  4.186E-01,
     5      4.171E-01,  4.159E-01,  4.145E-01,  4.130E-01,  4.115E-01,
     6      4.100E-01,  4.085E-01,  4.070E-01,  4.057E-01,  4.042E-01,
     7      4.028E-01,  4.014E-01,  3.998E-01,  3.982E-01,  3.967E-01,
     8      3.950E-01,  3.935E-01,  3.919E-01,  3.904E-01,  3.892E-01,
     9      3.878E-01,  3.863E-01,  3.848E-01,  3.833E-01,  3.818E-01,
     $      3.803E-01,  3.789E-01,  3.775E-01,  3.761E-01,  3.746E-01/
      DATA (X(I),I= 1801, 1850)/
     1      3.731E-01,  3.718E-01,  3.706E-01,  3.694E-01,  3.681E-01,
     2      3.669E-01,  3.657E-01,  3.646E-01,  3.635E-01,  3.624E-01,
     3      3.613E-01,  3.603E-01,  3.592E-01,  3.581E-01,  3.571E-01,
     4      3.561E-01,  3.550E-01,  3.540E-01,  3.530E-01,  3.520E-01,
     5      3.510E-01,  3.503E-01,  3.495E-01,  3.487E-01,  3.479E-01,
     6      3.472E-01,  3.464E-01,  3.457E-01,  3.450E-01,  3.443E-01,
     7      3.436E-01,  3.429E-01,  3.422E-01,  3.415E-01,  3.409E-01,
     8      3.403E-01,  3.398E-01,  3.392E-01,  3.386E-01,  3.380E-01,
     9      3.375E-01,  3.369E-01,  3.364E-01,  3.358E-01,  3.353E-01,
     $      3.347E-01,  3.342E-01,  3.337E-01,  3.332E-01,  3.327E-01/
      DATA (X(I),I= 1851, 1900)/
     1      3.323E-01,  3.318E-01,  3.313E-01,  3.308E-01,  3.304E-01,
     2      3.300E-01,  3.295E-01,  3.291E-01,  3.286E-01,  3.282E-01,
     3      3.278E-01,  3.275E-01,  3.271E-01,  3.267E-01,  3.264E-01,
     4      3.260E-01,  3.256E-01,  3.251E-01,  3.245E-01,  3.240E-01,
     5      3.235E-01,  3.230E-01,  3.224E-01,  3.219E-01,  3.213E-01,
     6      3.207E-01,  3.202E-01,  3.196E-01,  3.190E-01,  3.185E-01,
     7      3.179E-01,  3.174E-01,  3.169E-01,  3.163E-01,  3.158E-01,
     8      3.152E-01,  3.146E-01,  3.139E-01,  3.132E-01,  3.125E-01,
     9      3.117E-01,  3.110E-01,  3.102E-01,  3.095E-01,  3.087E-01,
     $      3.079E-01,  3.071E-01,  3.063E-01,  3.055E-01,  3.048E-01/
      DATA (X(I),I= 1901, 1950)/
     1      3.039E-01,  3.031E-01,  3.022E-01,  3.014E-01,  3.005E-01,
     2      2.996E-01,  2.988E-01,  2.979E-01,  2.970E-01,  2.961E-01,
     3      2.952E-01,  2.944E-01,  2.935E-01,  2.927E-01,  2.920E-01,
     4      2.913E-01,  2.906E-01,  2.900E-01,  2.893E-01,  2.886E-01,
     5      2.880E-01,  2.874E-01,  2.869E-01,  2.863E-01,  2.858E-01,
     6      2.852E-01,  2.847E-01,  2.842E-01,  2.838E-01,  2.834E-01,
     7      2.830E-01,  2.826E-01,  2.822E-01,  2.818E-01,  2.815E-01,
     8      2.813E-01,  2.811E-01,  2.809E-01,  2.807E-01,  2.805E-01,
     9      2.803E-01,  2.802E-01,  2.803E-01,  2.803E-01,  2.803E-01,
     $      2.803E-01,  2.803E-01,  2.804E-01,  2.804E-01,  2.805E-01/
      DATA (X(I),I= 1951, 2000)/
     1      2.806E-01,  2.807E-01,  2.808E-01,  2.809E-01,  2.810E-01,
     2      2.810E-01,  2.809E-01,  2.808E-01,  2.808E-01,  2.807E-01,
     3      2.806E-01,  2.805E-01,  2.804E-01,  2.801E-01,  2.799E-01,
     4      2.796E-01,  2.794E-01,  2.791E-01,  2.789E-01,  2.785E-01,
     5      2.780E-01,  2.775E-01,  2.770E-01,  2.765E-01,  2.760E-01,
     6      2.755E-01,  2.749E-01,  2.741E-01,  2.734E-01,  2.726E-01,
     7      2.718E-01,  2.710E-01,  2.703E-01,  2.694E-01,  2.684E-01,
     8      2.674E-01,  2.664E-01,  2.654E-01,  2.644E-01,  2.634E-01,
     9      2.624E-01,  2.612E-01,  2.601E-01,  2.589E-01,  2.578E-01,
     $      2.566E-01,  2.555E-01,  2.543E-01,  2.530E-01,  2.517E-01/
      DATA (X(I),I= 2001, 2050)/
     1      2.503E-01,  2.490E-01,  2.477E-01,  2.463E-01,  2.450E-01,
     2      2.436E-01,  2.423E-01,  2.410E-01,  2.396E-01,  2.383E-01,
     3      2.370E-01,  2.356E-01,  2.342E-01,  2.328E-01,  2.314E-01,
     4      2.299E-01,  2.285E-01,  2.271E-01,  2.257E-01,  2.243E-01,
     5      2.230E-01,  2.218E-01,  2.205E-01,  2.192E-01,  2.179E-01,
     6      2.166E-01,  2.153E-01,  2.140E-01,  2.126E-01,  2.113E-01,
     7      2.100E-01,  2.086E-01,  2.073E-01,  2.060E-01,  2.048E-01,
     8      2.036E-01,  2.025E-01,  2.013E-01,  2.001E-01,  1.989E-01,
     9      1.977E-01,  1.967E-01,  1.957E-01,  1.947E-01,  1.937E-01,
     $      1.927E-01,  1.917E-01,  1.907E-01,  1.898E-01,  1.890E-01/
      DATA (X(I),I= 2051, 2100)/
     1      1.883E-01,  1.875E-01,  1.867E-01,  1.859E-01,  1.851E-01,
     2      1.844E-01,  1.836E-01,  1.829E-01,  1.821E-01,  1.813E-01,
     3      1.806E-01,  1.798E-01,  1.791E-01,  1.786E-01,  1.781E-01,
     4      1.776E-01,  1.771E-01,  1.766E-01,  1.761E-01,  1.756E-01,
     5      1.751E-01,  1.747E-01,  1.743E-01,  1.739E-01,  1.735E-01,
     6      1.731E-01,  1.726E-01,  1.722E-01,  1.718E-01,  1.714E-01,
     7      1.710E-01,  1.706E-01,  1.703E-01,  1.698E-01,  1.695E-01,
     8      1.691E-01,  1.686E-01,  1.682E-01,  1.677E-01,  1.673E-01,
     9      1.668E-01,  1.664E-01,  1.660E-01,  1.655E-01,  1.650E-01,
     $      1.646E-01,  1.642E-01,  1.637E-01,  1.633E-01,  1.628E-01/
      DATA (X(I),I= 2101, 2150)/
     1      1.624E-01,  1.619E-01,  1.615E-01,  1.611E-01,  1.607E-01,
     2      1.602E-01,  1.598E-01,  1.593E-01,  1.588E-01,  1.583E-01,
     3      1.578E-01,  1.574E-01,  1.569E-01,  1.564E-01,  1.560E-01,
     4      1.555E-01,  1.551E-01,  1.548E-01,  1.544E-01,  1.540E-01,
     5      1.536E-01,  1.532E-01,  1.528E-01,  1.525E-01,  1.523E-01,
     6      1.520E-01,  1.518E-01,  1.516E-01,  1.514E-01,  1.511E-01,
     7      1.510E-01,  1.511E-01,  1.513E-01,  1.514E-01,  1.516E-01,
     8      1.518E-01,  1.519E-01,  1.521E-01,  1.523E-01,  1.526E-01,
     9      1.528E-01,  1.531E-01,  1.534E-01,  1.537E-01,  1.540E-01,
     $      1.543E-01,  1.547E-01,  1.551E-01,  1.555E-01,  1.560E-01/
      DATA (X(I),I= 2151, 2200)/
     1      1.564E-01,  1.568E-01,  1.572E-01,  1.575E-01,  1.579E-01,
     2      1.581E-01,  1.584E-01,  1.586E-01,  1.589E-01,  1.592E-01,
     3      1.594E-01,  1.596E-01,  1.598E-01,  1.599E-01,  1.600E-01,
     4      1.601E-01,  1.602E-01,  1.603E-01,  1.604E-01,  1.604E-01,
     5      1.602E-01,  1.600E-01,  1.598E-01,  1.596E-01,  1.594E-01,
     6      1.592E-01,  1.590E-01,  1.585E-01,  1.580E-01,  1.574E-01,
     7      1.568E-01,  1.562E-01,  1.555E-01,  1.549E-01,  1.543E-01,
     8      1.535E-01,  1.527E-01,  1.518E-01,  1.510E-01,  1.501E-01,
     9      1.493E-01,  1.484E-01,  1.475E-01,  1.464E-01,  1.453E-01,
     $      1.442E-01,  1.431E-01,  1.420E-01,  1.409E-01,  1.398E-01/
      DATA (X(I),I= 2201, 2250)/
     1      1.386E-01,  1.374E-01,  1.362E-01,  1.350E-01,  1.338E-01,
     2      1.326E-01,  1.314E-01,  1.302E-01,  1.290E-01,  1.278E-01,
     3      1.265E-01,  1.253E-01,  1.241E-01,  1.229E-01,  1.217E-01,
     4      1.204E-01,  1.192E-01,  1.180E-01,  1.169E-01,  1.157E-01,
     5      1.145E-01,  1.133E-01,  1.122E-01,  1.110E-01,  1.099E-01,
     6      1.089E-01,  1.079E-01,  1.069E-01,  1.060E-01,  1.050E-01,
     7      1.040E-01,  1.031E-01,  1.021E-01,  1.013E-01,  1.005E-01,
     8      9.973E-02,  9.897E-02,  9.820E-02,  9.743E-02,  9.664E-02,
     9      9.588E-02,  9.524E-02,  9.462E-02,  9.400E-02,  9.339E-02,
     $      9.279E-02,  9.217E-02,  9.158E-02,  9.098E-02,  9.049E-02/
      DATA (X(I),I= 2251, 2300)/
     1      9.002E-02,  8.958E-02,  8.913E-02,  8.869E-02,  8.827E-02,
     2      8.783E-02,  8.742E-02,  8.712E-02,  8.690E-02,  8.670E-02,
     3      8.648E-02,  8.629E-02,  8.607E-02,  8.588E-02,  8.568E-02,
     4      8.547E-02,  8.525E-02,  8.503E-02,  8.482E-02,  8.462E-02,
     5      8.440E-02,  8.418E-02,  8.397E-02,  8.379E-02,  8.369E-02,
     6      8.359E-02,  8.349E-02,  8.341E-02,  8.332E-02,  8.322E-02,
     7      8.316E-02,  8.305E-02,  8.288E-02,  8.269E-02,  8.251E-02,
     8      8.232E-02,  8.214E-02,  8.195E-02,  8.178E-02,  8.158E-02,
     9      8.133E-02,  8.108E-02,  8.083E-02,  8.057E-02,  8.031E-02,
     $      8.003E-02,  7.976E-02,  7.949E-02,  7.917E-02,  7.874E-02/
      DATA (X(I),I= 2301, 2350)/
     1      7.830E-02,  7.789E-02,  7.744E-02,  7.704E-02,  7.662E-02,
     2      7.620E-02,  7.579E-02,  7.549E-02,  7.519E-02,  7.490E-02,
     3      7.460E-02,  7.432E-02,  7.404E-02,  7.377E-02,  7.347E-02,
     4      7.333E-02,  7.329E-02,  7.323E-02,  7.318E-02,  7.315E-02,
     5      7.310E-02,  7.307E-02,  7.303E-02,  7.304E-02,  7.320E-02,
     6      7.336E-02,  7.351E-02,  7.369E-02,  7.386E-02,  7.402E-02,
     7      7.417E-02,  7.435E-02,  7.458E-02,  7.481E-02,  7.505E-02,
     8      7.529E-02,  7.556E-02,  7.582E-02,  7.607E-02,  7.634E-02,
     9      7.661E-02,  7.695E-02,  7.728E-02,  7.763E-02,  7.797E-02,
     $      7.830E-02,  7.864E-02,  7.895E-02,  7.929E-02,  7.952E-02/
      DATA (X(I),I= 2351, 2400)/
     1      7.975E-02,  7.996E-02,  8.018E-02,  8.039E-02,  8.061E-02,
     2      8.081E-02,  8.102E-02,  8.115E-02,  8.105E-02,  8.096E-02,
     3      8.086E-02,  8.074E-02,  8.062E-02,  8.048E-02,  8.033E-02,
     4      8.019E-02,  7.988E-02,  7.946E-02,  7.907E-02,  7.866E-02,
     5      7.824E-02,  7.784E-02,  7.743E-02,  7.702E-02,  7.658E-02,
     6      7.607E-02,  7.552E-02,  7.497E-02,  7.442E-02,  7.386E-02,
     7      7.328E-02,  7.271E-02,  7.214E-02,  7.149E-02,  7.072E-02,
     8      6.995E-02,  6.917E-02,  6.840E-02,  6.762E-02,  6.684E-02,
     9      6.605E-02,  6.527E-02,  6.453E-02,  6.382E-02,  6.311E-02,
     $      6.240E-02,  6.169E-02,  6.099E-02,  6.027E-02,  5.956E-02/
      DATA (X(I),I= 2401, 2450)/
     1      5.886E-02,  5.817E-02,  5.747E-02,  5.676E-02,  5.607E-02,
     2      5.537E-02,  5.469E-02,  5.400E-02,  5.331E-02,  5.265E-02,
     3      5.205E-02,  5.146E-02,  5.088E-02,  5.028E-02,  4.970E-02,
     4      4.912E-02,  4.853E-02,  4.794E-02,  4.738E-02,  4.685E-02,
     5      4.634E-02,  4.583E-02,  4.532E-02,  4.480E-02,  4.431E-02,
     6      4.380E-02,  4.331E-02,  4.291E-02,  4.260E-02,  4.230E-02,
     7      4.199E-02,  4.168E-02,  4.137E-02,  4.107E-02,  4.078E-02,
     8      4.047E-02,  4.029E-02,  4.016E-02,  4.004E-02,  3.991E-02,
     9      3.980E-02,  3.970E-02,  3.959E-02,  3.949E-02,  3.937E-02,
     $      3.934E-02,  3.933E-02,  3.933E-02,  3.934E-02,  3.935E-02/
      DATA (X(I),I= 2451, 2500)/
     1      3.935E-02,  3.936E-02,  3.936E-02,  3.936E-02,  3.931E-02,
     2      3.925E-02,  3.918E-02,  3.912E-02,  3.905E-02,  3.897E-02,
     3      3.889E-02,  3.881E-02,  3.874E-02,  3.866E-02,  3.855E-02,
     4      3.846E-02,  3.837E-02,  3.826E-02,  3.818E-02,  3.807E-02,
     5      3.795E-02,  3.786E-02,  3.769E-02,  3.748E-02,  3.727E-02,
     6      3.706E-02,  3.686E-02,  3.664E-02,  3.643E-02,  3.622E-02,
     7      3.601E-02,  3.581E-02,  3.561E-02,  3.542E-02,  3.522E-02,
     8      3.503E-02,  3.484E-02,  3.465E-02,  3.446E-02,  3.427E-02,
     9      3.407E-02,  3.386E-02,  3.364E-02,  3.343E-02,  3.322E-02,
     $      3.301E-02,  3.280E-02,  3.259E-02,  3.238E-02,  3.221E-02/
      DATA (X(I),I= 2501, 2550)/
     1      3.209E-02,  3.198E-02,  3.186E-02,  3.175E-02,  3.164E-02,
     2      3.153E-02,  3.143E-02,  3.132E-02,  3.126E-02,  3.136E-02,
     3      3.148E-02,  3.159E-02,  3.170E-02,  3.182E-02,  3.194E-02,
     4      3.206E-02,  3.219E-02,  3.232E-02,  3.253E-02,  3.275E-02,
     5      3.298E-02,  3.320E-02,  3.343E-02,  3.366E-02,  3.389E-02,
     6      3.412E-02,  3.435E-02,  3.455E-02,  3.475E-02,  3.495E-02,
     7      3.515E-02,  3.535E-02,  3.555E-02,  3.573E-02,  3.593E-02,
     8      3.612E-02,  3.625E-02,  3.628E-02,  3.631E-02,  3.634E-02,
     9      3.637E-02,  3.639E-02,  3.641E-02,  3.643E-02,  3.646E-02,
     $      3.646E-02,  3.636E-02,  3.623E-02,  3.611E-02,  3.599E-02/
      DATA (X(I),I= 2551, 2600)/
     1      3.586E-02,  3.572E-02,  3.559E-02,  3.545E-02,  3.531E-02,
     2      3.509E-02,  3.481E-02,  3.453E-02,  3.425E-02,  3.398E-02,
     3      3.369E-02,  3.341E-02,  3.312E-02,  3.284E-02,  3.254E-02,
     4      3.217E-02,  3.180E-02,  3.143E-02,  3.106E-02,  3.069E-02,
     5      3.031E-02,  2.994E-02,  2.956E-02,  2.919E-02,  2.882E-02,
     6      2.845E-02,  2.808E-02,  2.771E-02,  2.734E-02,  2.696E-02,
     7      2.660E-02,  2.622E-02,  2.585E-02,  2.549E-02,  2.512E-02,
     8      2.476E-02,  2.440E-02,  2.404E-02,  2.368E-02,  2.332E-02,
     9      2.297E-02,  2.261E-02,  2.225E-02,  2.193E-02,  2.164E-02,
     $      2.135E-02,  2.106E-02,  2.076E-02,  2.048E-02,  2.019E-02/
      DATA (X(I),I= 2601, 2650)/
     1      1.990E-02,  1.962E-02,  1.935E-02,  1.916E-02,  1.898E-02,
     2      1.881E-02,  1.863E-02,  1.846E-02,  1.828E-02,  1.811E-02,
     3      1.794E-02,  1.777E-02,  1.764E-02,  1.757E-02,  1.749E-02,
     4      1.742E-02,  1.735E-02,  1.728E-02,  1.721E-02,  1.714E-02,
     5      1.708E-02,  1.701E-02,  1.699E-02,  1.699E-02,  1.699E-02,
     6      1.699E-02,  1.699E-02,  1.699E-02,  1.699E-02,  1.699E-02,
     7      1.699E-02,  1.699E-02,  1.700E-02,  1.702E-02,  1.703E-02,
     8      1.704E-02,  1.706E-02,  1.707E-02,  1.708E-02,  1.709E-02,
     9      1.710E-02,  1.710E-02,  1.706E-02,  1.701E-02,  1.696E-02,
     $      1.692E-02,  1.687E-02,  1.683E-02,  1.678E-02,  1.673E-02/
      DATA (X(I),I= 2651, 2700)/
     1      1.668E-02,  1.661E-02,  1.651E-02,  1.642E-02,  1.632E-02,
     2      1.622E-02,  1.612E-02,  1.602E-02,  1.592E-02,  1.582E-02,
     3      1.572E-02,  1.560E-02,  1.545E-02,  1.531E-02,  1.517E-02,
     4      1.503E-02,  1.489E-02,  1.474E-02,  1.460E-02,  1.446E-02,
     5      1.432E-02,  1.420E-02,  1.408E-02,  1.397E-02,  1.386E-02,
     6      1.375E-02,  1.363E-02,  1.352E-02,  1.341E-02,  1.329E-02,
     7      1.318E-02,  1.313E-02,  1.310E-02,  1.308E-02,  1.305E-02,
     8      1.303E-02,  1.300E-02,  1.298E-02,  1.295E-02,  1.292E-02,
     9      1.290E-02,  1.293E-02,  1.297E-02,  1.302E-02,  1.307E-02,
     $      1.311E-02,  1.316E-02,  1.320E-02,  1.325E-02,  1.330E-02/
      DATA (X(I),I= 2701, 2750)/
     1      1.334E-02,  1.341E-02,  1.349E-02,  1.357E-02,  1.366E-02,
     2      1.374E-02,  1.382E-02,  1.390E-02,  1.398E-02,  1.406E-02,
     3      1.414E-02,  1.421E-02,  1.427E-02,  1.433E-02,  1.438E-02,
     4      1.444E-02,  1.450E-02,  1.456E-02,  1.462E-02,  1.467E-02,
     5      1.473E-02,  1.475E-02,  1.474E-02,  1.473E-02,  1.472E-02,
     6      1.471E-02,  1.470E-02,  1.468E-02,  1.467E-02,  1.466E-02,
     7      1.465E-02,  1.461E-02,  1.452E-02,  1.442E-02,  1.433E-02,
     8      1.423E-02,  1.414E-02,  1.405E-02,  1.395E-02,  1.386E-02,
     9      1.377E-02,  1.367E-02,  1.356E-02,  1.344E-02,  1.333E-02,
     $      1.321E-02,  1.310E-02,  1.298E-02,  1.287E-02,  1.275E-02/
      DATA (X(I),I= 2751, 2800)/
     1      1.264E-02,  1.252E-02,  1.236E-02,  1.220E-02,  1.204E-02,
     2      1.187E-02,  1.171E-02,  1.155E-02,  1.138E-02,  1.122E-02,
     3      1.106E-02,  1.089E-02,  1.072E-02,  1.055E-02,  1.038E-02,
     4      1.021E-02,  1.004E-02,  9.872E-03,  9.701E-03,  9.531E-03,
     5      9.361E-03,  9.190E-03,  9.029E-03,  8.896E-03,  8.763E-03,
     6      8.634E-03,  8.503E-03,  8.370E-03,  8.240E-03,  8.108E-03,
     7      7.977E-03,  7.847E-03,  7.717E-03,  7.622E-03,  7.540E-03,
     8      7.457E-03,  7.373E-03,  7.291E-03,  7.209E-03,  7.126E-03,
     9      7.041E-03,  6.961E-03,  6.878E-03,  6.813E-03,  6.782E-03,
     $      6.751E-03,  6.721E-03,  6.690E-03,  6.658E-03,  6.628E-03/
      DATA (X(I),I= 2801, 2850)/
     1      6.599E-03,  6.567E-03,  6.536E-03,  6.508E-03,  6.499E-03,
     2      6.495E-03,  6.493E-03,  6.490E-03,  6.488E-03,  6.484E-03,
     3      6.480E-03,  6.478E-03,  6.474E-03,  6.470E-03,  6.471E-03,
     4      6.473E-03,  6.476E-03,  6.480E-03,  6.482E-03,  6.486E-03,
     5      6.489E-03,  6.493E-03,  6.496E-03,  6.499E-03,  6.501E-03,
     6      6.488E-03,  6.465E-03,  6.444E-03,  6.422E-03,  6.401E-03,
     7      6.381E-03,  6.359E-03,  6.337E-03,  6.316E-03,  6.294E-03,
     8      6.266E-03,  6.203E-03,  6.135E-03,  6.067E-03,  6.002E-03,
     9      5.932E-03,  5.867E-03,  5.797E-03,  5.732E-03,  5.664E-03,
     $      5.596E-03,  5.528E-03,  5.457E-03,  5.388E-03,  5.318E-03/
      DATA (X(I),I= 2851, 2900)/
     1      5.248E-03,  5.177E-03,  5.107E-03,  5.036E-03,  4.966E-03,
     2      4.895E-03,  4.825E-03,  4.781E-03,  4.755E-03,  4.729E-03,
     3      4.703E-03,  4.675E-03,  4.648E-03,  4.622E-03,  4.595E-03,
     4      4.569E-03,  4.542E-03,  4.516E-03,  4.514E-03,  4.517E-03,
     5      4.520E-03,  4.523E-03,  4.527E-03,  4.530E-03,  4.533E-03,
     6      4.539E-03,  4.541E-03,  4.545E-03,  4.550E-03,  4.568E-03,
     7      4.588E-03,  4.609E-03,  4.628E-03,  4.649E-03,  4.669E-03,
     8      4.689E-03,  4.710E-03,  4.731E-03,  4.750E-03,  4.771E-03,
     9      4.796E-03,  4.822E-03,  4.847E-03,  4.869E-03,  4.896E-03,
     $      4.921E-03,  4.945E-03,  4.970E-03,  4.995E-03,  5.020E-03/
      DATA (X(I),I= 2901, 2950)/
     1      5.038E-03,  5.040E-03,  5.038E-03,  5.037E-03,  5.036E-03,
     2      5.036E-03,  5.034E-03,  5.034E-03,  5.031E-03,  5.031E-03,
     3      5.031E-03,  5.019E-03,  4.979E-03,  4.934E-03,  4.892E-03,
     4      4.848E-03,  4.805E-03,  4.763E-03,  4.718E-03,  4.676E-03,
     5      4.632E-03,  4.590E-03,  4.541E-03,  4.475E-03,  4.405E-03,
     6      4.336E-03,  4.268E-03,  4.198E-03,  4.130E-03,  4.060E-03,
     7      3.990E-03,  3.922E-03,  3.852E-03,  3.782E-03,  3.715E-03,
     8      3.646E-03,  3.577E-03,  3.508E-03,  3.439E-03,  3.370E-03,
     9      3.301E-03,  3.232E-03,  3.163E-03,  3.094E-03,  3.026E-03,
     $      2.971E-03,  2.919E-03,  2.868E-03,  2.816E-03,  2.764E-03/
      DATA (X(I),I= 2951, 3000)/
     1      2.712E-03,  2.661E-03,  2.609E-03,  2.557E-03,  2.505E-03,
     2      2.454E-03,  2.416E-03,  2.386E-03,  2.356E-03,  2.326E-03,
     3      2.297E-03,  2.267E-03,  2.237E-03,  2.207E-03,  2.177E-03,
     4      2.148E-03,  2.118E-03,  2.096E-03,  2.087E-03,  2.078E-03,
     5      2.070E-03,  2.061E-03,  2.052E-03,  2.043E-03,  2.034E-03,
     6      2.025E-03,  2.016E-03,  2.007E-03,  2.000E-03,  2.000E-03,
     7      2.001E-03,  2.002E-03,  2.003E-03,  2.004E-03,  2.005E-03,
     8      2.006E-03,  2.007E-03,  2.007E-03,  2.008E-03,  2.009E-03,
     9      2.008E-03,  2.006E-03,  2.003E-03,  2.001E-03,  1.999E-03,
     $      1.997E-03,  1.994E-03,  1.992E-03,  1.990E-03,  1.988E-03/
      DATA (X(I),I= 3001, 3050)/
     1      1.985E-03,  1.980E-03,  1.968E-03,  1.956E-03,  1.944E-03,
     2      1.932E-03,  1.919E-03,  1.907E-03,  1.895E-03,  1.883E-03,
     3      1.871E-03,  1.859E-03,  1.846E-03,  1.827E-03,  1.805E-03,
     4      1.783E-03,  1.761E-03,  1.740E-03,  1.718E-03,  1.696E-03,
     5      1.674E-03,  1.652E-03,  1.631E-03,  1.609E-03,  1.585E-03,
     6      1.557E-03,  1.529E-03,  1.500E-03,  1.472E-03,  1.444E-03,
     7      1.416E-03,  1.388E-03,  1.359E-03,  1.331E-03,  1.303E-03,
     8      1.275E-03,  1.251E-03,  1.228E-03,  1.206E-03,  1.183E-03,
     9      1.161E-03,  1.139E-03,  1.116E-03,  1.094E-03,  1.072E-03,
     $      1.049E-03,  1.027E-03,  1.007E-03,  1.003E-03,  9.991E-04/
      DATA (X(I),I= 3051, 3100)/
     1      9.959E-04,  9.926E-04,  9.891E-04,  9.857E-04,  9.826E-04,
     2      9.790E-04,  9.758E-04,  9.725E-04,  9.692E-04,  9.720E-04,
     3      9.842E-04,  9.967E-04,  1.009E-03,  1.022E-03,  1.034E-03,
     4      1.047E-03,  1.059E-03,  1.071E-03,  1.084E-03,  1.096E-03,
     5      1.109E-03,  1.118E-03,  1.126E-03,  1.133E-03,  1.141E-03,
     6      1.148E-03,  1.156E-03,  1.163E-03,  1.171E-03,  1.178E-03,
     7      1.186E-03,  1.193E-03,  1.200E-03,  1.194E-03,  1.185E-03,
     8      1.176E-03,  1.166E-03,  1.157E-03,  1.148E-03,  1.138E-03,
     9      1.129E-03,  1.120E-03,  1.110E-03,  1.101E-03,  1.091E-03,
     $      1.076E-03,  1.060E-03,  1.044E-03,  1.028E-03,  1.012E-03/
      DATA (X(I),I= 3101, 3150)/
     1      9.955E-04,  9.794E-04,  9.632E-04,  9.473E-04,  9.310E-04,
     2      9.151E-04,  8.982E-04,  8.770E-04,  8.556E-04,  8.339E-04,
     3      8.121E-04,  7.907E-04,  7.689E-04,  7.471E-04,  7.255E-04,
     4      7.038E-04,  6.822E-04,  6.606E-04,  6.376E-04,  6.083E-04,
     5      5.783E-04,  5.478E-04,  5.177E-04,  4.875E-04,  4.574E-04,
     6      4.272E-04,  3.969E-04,  3.668E-04,  3.365E-04,  3.063E-04,
     7      2.750E-04,  2.500E-04,  2.250E-04,  2.000E-04,  1.850E-04,
     8      1.700E-04,  1.550E-04,  1.400E-04,  1.250E-04,  1.100E-04,
     9      0.950E-04,  0.825E-04,  0.700E-04,  0.575E-04,  0.400E-04,
     *      0.275E-04,  0.175E-04,  0.100E-04,  0.040E-04,  0.00000  /
C
      DATA (Y(I),I=    1,   50)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 51, 100)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  101,  150)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 151, 200)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  201,  250)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 251, 300)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  301,  350)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 351, 400)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  401,  450)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 451, 500)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  501,  550)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 551, 600)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  601,  650)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 651, 700)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  701,  750)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I= 751, 800)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Y(I),I=  801,  850)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.025e-05,  0.050e-05,  0.075e-05,  0.100e-05,  0.125e-05,
     6      0.150e-05,  0.175e-05,  0.225e-05,  0.250e-05,  0.300e-05,
     7      0.350e-05,  0.400e-05,  0.500e-05,  0.650e-05,  0.800e-05,
     8      1.100e-05,  1.400e-05,  1.800e-05,  2.250e-05,  2.650e-05,
     9      3.000e-05,  3.104E-05,  3.136E-05,  3.152E-05,  3.186E-05,
     $      3.213E-05,  3.229E-05,  3.206E-05,  3.156E-05,  3.063E-05/
      DATA (Y(I),I= 851, 900)/
     1      3.098E-05,  3.197E-05,  3.271E-05,  3.315E-05,  3.262E-05,
     2      3.201E-05,  3.129E-05,  3.148E-05,  3.206E-05,  3.175E-05,
     3      3.148E-05,  3.167E-05,  3.159E-05,  3.120E-05,  3.117E-05,
     4      3.109E-05,  3.041E-05,  2.995E-05,  3.011E-05,  3.004E-05,
     5      2.972E-05,  2.933E-05,  2.887E-05,  2.801E-05,  2.731E-05,
     6      2.735E-05,  2.656E-05,  2.712E-05,  2.576E-05,  2.565E-05,
     7      2.449E-05,  2.450E-05,  2.454E-05,  2.414E-05,  2.488E-05,
     8      2.413E-05,  2.297E-05,  2.298E-05,  2.335E-05,  2.259E-05,
     9      2.142E-05,  2.257E-05,  2.259E-05,  2.220E-05,  2.259E-05,
     $      2.336E-05,  2.299E-05,  2.262E-05,  2.298E-05,  2.298E-05/
      DATA (Y(I),I=  901,  950)/
     1      2.298E-05,  2.337E-05,  2.414E-05,  2.491E-05,  2.570E-05,
     2      2.492E-05,  2.531E-05,  2.567E-05,  2.647E-05,  2.839E-05,
     3      2.916E-05,  2.917E-05,  2.916E-05,  2.956E-05,  3.033E-05,
     4      3.110E-05,  3.108E-05,  3.071E-05,  3.226E-05,  3.303E-05,
     5      3.302E-05,  3.417E-05,  3.419E-05,  3.458E-05,  3.458E-05,
     6      3.457E-05,  3.613E-05,  3.769E-05,  3.768E-05,  3.845E-05,
     7      3.843E-05,  3.882E-05,  3.920E-05,  3.959E-05,  3.998E-05,
     8      4.037E-05,  4.036E-05,  4.074E-05,  4.074E-05,  4.036E-05,
     9      4.115E-05,  4.192E-05,  4.192E-05,  4.269E-05,  4.270E-05,
     $      4.346E-05,  4.346E-05,  4.346E-05,  4.463E-05,  4.423E-05/
      DATA (Y(I),I= 951, 1000)/
     1      4.502E-05,  4.539E-05,  4.501E-05,  4.465E-05,  4.388E-05,
     2      4.503E-05,  4.578E-05,  4.617E-05,  4.657E-05,  4.657E-05,
     3      4.617E-05,  4.657E-05,  4.655E-05,  4.695E-05,  4.808E-05,
     4      4.807E-05,  4.768E-05,  4.771E-05,  4.807E-05,  4.808E-05,
     5      4.808E-05,  4.731E-05,  4.731E-05,  4.772E-05,  4.808E-05,
     6      4.846E-05,  4.769E-05,  4.807E-05,  4.846E-05,  4.847E-05,
     7      4.845E-05,  4.885E-05,  4.886E-05,  4.848E-05,  4.810E-05,
     8      4.848E-05,  4.769E-05,  4.773E-05,  4.656E-05,  4.690E-05,
     9      4.770E-05,  4.730E-05,  4.690E-05,  4.728E-05,  4.687E-05,
     $      4.728E-05,  4.649E-05,  4.730E-05,  4.690E-05,  4.766E-05/
      DATA (Y(I),I= 1001, 1050)/
     1      4.804E-05,  4.727E-05,  4.649E-05,  4.727E-05,  4.647E-05,
     2      4.723E-05,  4.760E-05,  4.800E-05,  4.836E-05,  4.762E-05,
     3      4.760E-05,  4.834E-05,  4.837E-05,  4.796E-05,  4.796E-05,
     4      4.796E-05,  4.876E-05,  4.955E-05,  5.030E-05,  5.067E-05,
     5      5.106E-05,  5.183E-05,  5.104E-05,  5.185E-05,  5.224E-05,
     6      5.340E-05,  5.301E-05,  5.340E-05,  5.456E-05,  5.532E-05,
     7      5.532E-05,  5.531E-05,  5.608E-05,  5.762E-05,  5.799E-05,
     8      5.877E-05,  5.837E-05,  5.839E-05,  5.878E-05,  5.916E-05,
     9      6.032E-05,  6.109E-05,  6.226E-05,  6.342E-05,  6.305E-05,
     $      6.305E-05,  6.265E-05,  6.226E-05,  6.298E-05,  6.452E-05/
      DATA (Y(I),I= 1051, 1100)/
     1      6.531E-05,  6.572E-05,  6.572E-05,  6.614E-05,  6.611E-05,
     2      6.574E-05,  6.690E-05,  6.648E-05,  6.688E-05,  6.688E-05,
     3      6.764E-05,  6.841E-05,  6.843E-05,  6.918E-05,  6.919E-05,
     4      6.959E-05,  6.955E-05,  6.993E-05,  6.990E-05,  7.031E-05,
     5      6.914E-05,  6.835E-05,  6.912E-05,  7.105E-05,  7.222E-05,
     6      7.261E-05,  7.258E-05,  7.258E-05,  7.258E-05,  7.260E-05,
     7      7.179E-05,  7.257E-05,  7.180E-05,  7.142E-05,  7.181E-05,
     8      7.217E-05,  7.256E-05,  7.217E-05,  7.217E-05,  7.219E-05,
     9      7.256E-05,  7.180E-05,  7.178E-05,  7.255E-05,  7.216E-05,
     $      7.175E-05,  7.175E-05,  7.137E-05,  7.177E-05,  7.214E-05/
      DATA (Y(I),I= 1101, 1150)/
     1      7.175E-05,  7.134E-05,  7.094E-05,  7.014E-05,  7.014E-05,
     2      7.016E-05,  7.013E-05,  6.937E-05,  6.938E-05,  6.900E-05,
     3      6.863E-05,  6.860E-05,  6.859E-05,  6.741E-05,  6.700E-05,
     4      6.700E-05,  6.740E-05,  6.742E-05,  6.704E-05,  6.665E-05,
     5      6.549E-05,  6.432E-05,  6.427E-05,  6.463E-05,  6.581E-05,
     6      6.662E-05,  6.625E-05,  6.552E-05,  6.473E-05,  6.473E-05,
     7      6.471E-05,  6.469E-05,  6.391E-05,  6.433E-05,  6.475E-05,
     8      6.513E-05,  6.554E-05,  6.552E-05,  6.592E-05,  6.552E-05,
     9      6.553E-05,  6.476E-05,  6.440E-05,  6.483E-05,  6.443E-05,
     $      6.443E-05,  6.522E-05,  6.522E-05,  6.599E-05,  6.599E-05/
      DATA (Y(I),I= 1151, 1200)/
     1      6.560E-05,  6.520E-05,  6.521E-05,  6.558E-05,  6.556E-05,
     2      6.557E-05,  6.596E-05,  6.596E-05,  6.596E-05,  6.634E-05,
     3      6.595E-05,  6.555E-05,  6.555E-05,  6.594E-05,  6.515E-05,
     4      6.553E-05,  6.590E-05,  6.551E-05,  6.474E-05,  6.549E-05,
     5      6.589E-05,  6.628E-05,  6.709E-05,  6.672E-05,  6.634E-05,
     6      6.555E-05,  6.478E-05,  6.475E-05,  6.395E-05,  6.472E-05,
     7      6.472E-05,  6.510E-05,  6.509E-05,  6.508E-05,  6.511E-05,
     8      6.434E-05,  6.357E-05,  6.356E-05,  6.356E-05,  6.278E-05,
     9      6.239E-05,  6.161E-05,  6.199E-05,  6.163E-05,  6.161E-05,
     $      6.161E-05,  6.085E-05,  5.931E-05,  5.888E-05,  5.888E-05/
      DATA (Y(I),I= 1201, 1250)/
     1      5.846E-05,  5.805E-05,  5.809E-05,  5.890E-05,  5.892E-05,
     2      5.896E-05,  5.823E-05,  5.746E-05,  5.709E-05,  5.708E-05,
     3      5.707E-05,  5.668E-05,  5.709E-05,  5.666E-05,  5.627E-05,
     4      5.628E-05,  5.552E-05,  5.550E-05,  5.552E-05,  5.588E-05,
     5      5.627E-05,  5.588E-05,  5.627E-05,  5.591E-05,  5.477E-05,
     6      5.515E-05,  5.517E-05,  5.558E-05,  5.599E-05,  5.591E-05,
     7      5.668E-05,  5.507E-05,  5.429E-05,  5.311E-05,  5.279E-05,
     8      5.281E-05,  5.402E-05,  5.406E-05,  5.485E-05,  5.562E-05,
     9      5.562E-05,  5.481E-05,  5.484E-05,  5.483E-05,  5.443E-05,
     $      5.404E-05,  5.365E-05,  5.442E-05,  5.444E-05,  5.521E-05/
      DATA (Y(I),I= 1251, 1300)/
     1      5.481E-05,  5.443E-05,  5.363E-05,  5.286E-05,  5.283E-05,
     2      5.358E-05,  5.356E-05,  5.317E-05,  5.240E-05,  5.279E-05,
     3      5.240E-05,  5.241E-05,  5.321E-05,  5.322E-05,  5.327E-05,
     4      5.329E-05,  5.331E-05,  5.329E-05,  5.290E-05,  5.251E-05,
     5      5.213E-05,  5.135E-05,  5.057E-05,  5.056E-05,  5.056E-05,
     6      5.094E-05,  5.133E-05,  5.133E-05,  5.132E-05,  5.133E-05,
     7      5.054E-05,  4.934E-05,  4.815E-05,  4.733E-05,  4.735E-05,
     8      4.777E-05,  4.819E-05,  4.784E-05,  4.825E-05,  4.790E-05,
     9      4.830E-05,  4.793E-05,  4.873E-05,  4.874E-05,  4.790E-05,
     $      4.709E-05,  4.627E-05,  4.467E-05,  4.427E-05,  4.469E-05/
      DATA (Y(I),I= 1301, 1350)/
     1      4.434E-05,  4.431E-05,  4.473E-05,  4.477E-05,  4.365E-05,
     2      4.369E-05,  4.527E-05,  4.528E-05,  4.603E-05,  4.641E-05,
     3      4.520E-05,  4.480E-05,  4.363E-05,  4.404E-05,  4.444E-05,
     4      4.408E-05,  4.448E-05,  4.452E-05,  4.336E-05,  4.184E-05,
     5      4.223E-05,  4.263E-05,  4.301E-05,  4.262E-05,  4.299E-05,
     6      4.143E-05,  4.105E-05,  3.989E-05,  3.995E-05,  3.919E-05,
     7      3.999E-05,  4.118E-05,  4.039E-05,  3.997E-05,  3.880E-05,
     8      3.762E-05,  3.645E-05,  3.493E-05,  3.337E-05,  3.222E-05,
     9      3.296E-05,  3.412E-05,  3.255E-05,  3.257E-05,  3.103E-05,
     $      3.029E-05,  2.876E-05,  2.878E-05,  2.800E-05,  2.883E-05/
      DATA (Y(I),I= 1351, 1400)/
     1      2.881E-05,  2.806E-05,  2.729E-05,  2.650E-05,  2.574E-05,
     2      2.380E-05,  2.262E-05,  2.108E-05,  2.031E-05,  1.842E-05,
     3      1.765E-05,  1.648E-05,  1.646E-05,  1.685E-05,  1.529E-05,
     4      1.451E-05,  1.258E-05,  1.104E-05,  9.506E-06,  9.546E-06,
     5      8.010E-06,  6.431E-06,  4.851E-06,  4.067E-06,  2.472E-06,
     6      8.919E-07, -2.698E-07, -2.356E-07, -6.024E-07, -1.335E-06,
     7     -2.450E-06, -3.996E-06, -5.582E-06, -6.779E-06, -7.956E-06,
     8     -9.542E-06, -1.072E-05, -1.150E-05, -1.227E-05, -1.305E-05,
     9     -1.382E-05, -1.500E-05, -1.580E-05, -1.738E-05, -1.779E-05,
     $     -1.823E-05, -1.900E-05, -2.050E-05, -2.086E-05, -2.159E-05/
      DATA (Y(I),I= 1401, 1450)/
     1     -2.191E-05, -2.268E-05, -2.308E-05, -2.276E-05, -2.393E-05,
     2     -2.551E-05, -2.669E-05, -2.709E-05, -2.632E-05, -2.709E-05,
     3     -2.709E-05, -2.749E-05, -3.016E-05, -2.709E-05, -2.709E-05,
     4     -2.709E-05, -2.708E-05, -2.709E-05, -2.709E-05, -2.709E-05,
     5     -2.709E-05, -2.709E-05, -2.709E-05, -2.709E-05, -2.709E-05,
     6     -3.113E-05, -2.709E-05, -2.709E-05, -3.477E-05, -2.709E-05,
     7     -2.685E-05, -3.453E-05, -2.684E-05, -2.685E-05, -3.065E-05,
     8     -2.684E-05, -2.684E-05, -3.065E-05, -2.685E-05, -2.685E-05,
     9     -2.685E-05, -2.684E-05, -2.684E-05, -3.065E-05, -2.685E-05,
     $     -3.453E-05, -2.685E-05, -2.685E-05, -2.684E-05, -2.685E-05/
      DATA (Y(I),I= 1451, 1500)/
     1     -2.684E-05, -2.684E-05, -2.685E-05, -3.086E-05, -3.466E-05,
     2     -3.467E-05, -3.453E-05, -3.854E-05, -3.854E-05, -3.072E-05,
     3     -3.854E-05, -3.840E-05, -3.854E-05, -3.840E-05, -3.452E-05,
     4     -4.221E-05, -3.840E-05, -4.221E-05, -4.221E-05, -3.833E-05,
     5     -4.220E-05, -4.221E-05, -3.833E-05, -4.221E-05, -4.220E-05,
     6     -4.221E-05, -4.221E-05, -3.833E-05, -3.833E-05, -3.833E-05,
     7     -3.834E-05, -4.601E-05, -4.601E-05, -3.833E-05, -3.833E-05,
     8     -4.221E-05, -4.221E-05, -4.601E-05, -4.221E-05, -4.221E-05,
     9     -4.221E-05, -4.626E-05, -4.625E-05, -4.650E-05, -5.054E-05,
     $     -4.650E-05, -5.054E-05, -5.040E-05, -4.650E-05, -5.041E-05/
      DATA (Y(I),I= 1501, 1550)/
     1     -5.027E-05, -4.637E-05, -4.636E-05, -5.041E-05, -5.026E-05,
     2     -5.405E-05, -5.404E-05, -5.027E-05, -5.041E-05, -5.041E-05,
     3     -4.637E-05, -5.027E-05, -5.027E-05, -5.432E-05, -5.028E-05,
     4     -5.419E-05, -6.200E-05, -5.418E-05, -6.224E-05, -5.418E-05,
     5     -5.443E-05, -5.442E-05, -6.248E-05, -6.248E-05, -5.480E-05,
     6     -5.480E-05, -6.249E-05, -6.248E-05, -5.443E-05, -5.442E-05,
     7     -5.834E-05, -6.235E-05, -6.222E-05, -6.221E-05, -5.429E-05,
     8     -5.819E-05, -5.452E-05, -5.429E-05, -5.429E-05, -6.220E-05,
     9     -5.453E-05, -5.453E-05, -5.452E-05, -5.453E-05, -5.453E-05,
     $     -6.246E-05, -5.477E-05, -5.845E-05, -5.477E-05, -5.477E-05/
      DATA (Y(I),I= 1551, 1600)/
     1     -5.478E-05, -5.453E-05, -4.671E-05, -4.672E-05, -4.633E-05,
     2     -4.633E-05, -4.647E-05, -4.671E-05, -4.633E-05, -4.647E-05,
     3     -4.269E-05, -4.672E-05, -4.671E-05, -4.671E-05, -4.671E-05,
     4     -4.671E-05, -4.671E-05, -4.695E-05, -4.291E-05, -4.304E-05,
     5     -4.266E-05, -4.280E-05, -3.899E-05, -3.512E-05, -3.875E-05,
     6     -3.488E-05, -3.107E-05, -2.720E-05, -2.720E-05, -2.719E-05,
     7     -2.720E-05, -2.719E-05, -2.720E-05, -2.719E-05, -3.121E-05,
     8     -3.120E-05, -3.121E-05, -2.353E-05, -2.755E-05, -1.961E-05,
     9     -2.730E-05, -2.730E-05, -2.328E-05, -1.156E-05, -1.169E-05,
     $     -1.561E-05, -1.170E-05, -1.938E-05, -1.951E-05, -1.975E-05/
      DATA (Y(I),I= 1601, 1650)/
     1     -1.574E-05, -1.975E-05, -1.951E-05, -1.951E-05, -1.183E-05,
     2     -1.183E-05, -7.924E-06, -7.916E-06, -7.921E-06, -7.921E-06,
     3     -7.924E-06, -7.925E-06, -1.194E-05, -7.919E-06, -1.169E-05,
     4     -7.920E-06, -4.014E-06, -1.193E-05, -8.024E-06, -7.778E-06,
     5     -7.778E-06, -7.778E-06, -7.778E-06, -7.540E-06, -1.145E-05,
     6     -1.145E-05, -1.145E-05, -3.770E-06, -3.770E-06, -7.538E-06,
     7     -3.765E-06, -1.145E-05, -7.543E-06, -3.772E-06, -3.772E-06,
     8     -3.768E-06, -3.769E-06, -3.769E-06, -3.769E-06, -3.767E-06,
     9     -3.766E-06, -3.769E-06, -3.772E-06, -1.159E-05, -1.145E-05,
     $     -7.817E-06, -7.815E-06, -7.816E-06, -1.563E-05, -7.577E-06/
      DATA (Y(I),I= 1651, 1700)/
     1     -7.570E-06, -1.525E-05, -7.328E-06, -7.572E-06, -7.574E-06,
     2     -1.525E-05, -1.120E-05, -1.145E-05, -1.145E-05, -1.145E-05,
     3     -1.145E-05, -1.913E-05, -1.145E-05, -1.145E-05, -1.145E-05,
     4     -1.145E-05, -1.145E-05, -1.145E-05, -1.145E-05, -2.304E-05,
     5     -2.304E-05, -2.304E-05, -2.304E-05, -1.536E-05, -1.536E-05,
     6     -1.536E-05, -2.304E-05, -2.304E-05, -1.924E-05, -1.923E-05,
     7     -2.692E-05, -2.692E-05, -1.924E-05, -2.691E-05, -1.923E-05,
     8     -2.716E-05, -2.328E-05, -2.716E-05, -3.096E-05, -3.096E-05,
     9     -3.096E-05, -3.120E-05, -3.120E-05, -3.107E-05, -2.301E-05,
     $     -3.107E-05, -2.315E-05, -2.315E-05, -2.315E-05, -2.301E-05/
      DATA (Y(I),I= 1701, 1750)/
     1     -2.301E-05, -3.069E-05, -2.301E-05, -3.069E-05, -2.301E-05,
     2     -3.093E-05, -3.055E-05, -3.056E-05, -2.287E-05, -2.301E-05,
     3     -3.069E-05, -3.093E-05, -3.093E-05, -3.069E-05, -2.301E-05,
     4     -2.301E-05, -2.301E-05, -2.301E-05, -3.107E-05, -3.106E-05,
     5     -2.301E-05, -2.314E-05, -3.083E-05, -3.083E-05, -3.107E-05,
     6     -2.706E-05, -3.083E-05, -2.681E-05, -2.301E-05, -2.682E-05,
     7     -2.706E-05, -1.913E-05, -2.682E-05, -1.914E-05, -2.277E-05,
     8     -1.913E-05, -1.522E-05, -1.509E-05, -1.522E-05, -2.290E-05,
     9     -1.509E-05, -2.290E-05, -2.290E-05, -1.522E-05, -1.522E-05,
     $     -1.910E-05, -1.522E-05, -1.522E-05, -1.155E-05, -1.155E-05/
      DATA (Y(I),I= 1751, 1800)/
     1     -1.156E-05, -1.142E-05, -7.541E-06, -7.536E-06, -7.544E-06,
     2     -1.142E-05, -1.522E-05, -7.542E-06, -7.538E-06, -7.539E-06,
     3     -7.541E-06, -7.541E-06, -7.545E-06, -7.538E-06, -3.636E-06,
     4      1.426E-07, -3.494E-06, -3.633E-06, -3.630E-06, -3.631E-06,
     5     -3.634E-06,  4.149E-07, -3.629E-06, -3.629E-06, -3.633E-06,
     6     -3.632E-06, -3.634E-06, -3.628E-06, -3.385E-06,  4.292E-06,
     7      2.437E-07,  2.431E-07,  4.155E-06, -3.766E-06, -3.771E-06,
     8      2.430E-07,  2.438E-07,  3.915E-06,  3.909E-06,  7.824E-06,
     9      3.913E-06,  3.911E-06,  3.913E-06,  3.909E-06,  3.910E-06,
     $      3.910E-06,  3.914E-06,  4.150E-06,  1.184E-05,  1.184E-05/
      DATA (Y(I),I= 1801, 1850)/
     1      1.183E-05,  1.108E-05,  1.031E-05,  9.561E-06,  9.950E-06,
     2      1.227E-05,  1.420E-05,  1.534E-05,  1.572E-05,  1.571E-05,
     3      1.492E-05,  1.490E-05,  1.445E-05,  1.442E-05,  1.525E-05,
     4      1.487E-05,  1.452E-05,  1.604E-05,  1.840E-05,  2.070E-05,
     5      2.146E-05,  2.144E-05,  2.105E-05,  2.104E-05,  2.101E-05,
     6      2.100E-05,  2.175E-05,  2.095E-05,  2.097E-05,  1.981E-05,
     7      1.982E-05,  2.058E-05,  2.095E-05,  2.098E-05,  1.979E-05,
     8      1.860E-05,  1.782E-05,  1.934E-05,  1.972E-05,  2.083E-05,
     9      2.201E-05,  2.281E-05,  2.476E-05,  2.480E-05,  2.480E-05,
     $      2.407E-05,  2.371E-05,  2.295E-05,  2.218E-05,  2.141E-05/
      DATA (Y(I),I= 1851, 1900)/
     1      2.218E-05,  2.294E-05,  2.371E-05,  2.369E-05,  2.328E-05,
     2      2.287E-05,  2.169E-05,  2.128E-05,  2.126E-05,  2.004E-05,
     3      2.004E-05,  2.085E-05,  2.045E-05,  2.045E-05,  2.006E-05,
     4      2.085E-05,  2.082E-05,  2.121E-05,  2.198E-05,  2.275E-05,
     5      2.238E-05,  2.083E-05,  2.046E-05,  1.928E-05,  1.893E-05,
     6      1.815E-05,  1.815E-05,  1.812E-05,  1.929E-05,  2.046E-05,
     7      2.087E-05,  2.244E-05,  2.283E-05,  2.322E-05,  2.327E-05,
     8      2.291E-05,  2.329E-05,  2.482E-05,  2.558E-05,  2.593E-05,
     9      2.514E-05,  2.436E-05,  2.279E-05,  2.282E-05,  2.363E-05,
     $      2.441E-05,  2.524E-05,  2.599E-05,  2.759E-05,  2.874E-05/
      DATA (Y(I),I= 1901, 1950)/
     1      2.954E-05,  3.028E-05,  3.064E-05,  3.063E-05,  3.021E-05,
     2      2.866E-05,  2.788E-05,  2.788E-05,  2.749E-05,  2.711E-05,
     3      2.866E-05,  3.016E-05,  3.250E-05,  3.323E-05,  3.366E-05,
     4      3.328E-05,  3.367E-05,  3.366E-05,  3.483E-05,  3.521E-05,
     5      3.599E-05,  3.603E-05,  3.606E-05,  3.685E-05,  3.764E-05,
     6      3.960E-05,  4.076E-05,  4.079E-05,  4.041E-05,  3.926E-05,
     7      3.848E-05,  4.004E-05,  4.161E-05,  4.350E-05,  4.426E-05,
     8      4.386E-05,  4.305E-05,  4.344E-05,  4.457E-05,  4.570E-05,
     9      4.685E-05,  4.567E-05,  4.450E-05,  4.337E-05,  4.297E-05,
     $      4.414E-05,  4.527E-05,  4.644E-05,  4.721E-05,  4.798E-05/
      DATA (Y(I),I= 1951, 2000)/
     1      4.873E-05,  4.952E-05,  4.797E-05,  4.798E-05,  4.681E-05,
     2      4.721E-05,  4.760E-05,  4.724E-05,  4.802E-05,  4.805E-05,
     3      4.766E-05,  4.690E-05,  4.727E-05,  4.800E-05,  4.911E-05,
     4      4.947E-05,  4.636E-05,  4.401E-05,  4.164E-05,  4.089E-05,
     5      4.167E-05,  4.245E-05,  4.248E-05,  4.136E-05,  3.944E-05,
     6      3.832E-05,  3.753E-05,  3.867E-05,  3.982E-05,  4.022E-05,
     7      3.829E-05,  3.714E-05,  3.481E-05,  3.365E-05,  3.438E-05,
     8      3.478E-05,  3.474E-05,  3.395E-05,  3.277E-05,  3.085E-05,
     9      2.966E-05,  3.045E-05,  3.044E-05,  3.123E-05,  3.042E-05,
     $      2.965E-05,  2.969E-05,  2.966E-05,  3.082E-05,  3.273E-05/
      DATA (Y(I),I= 2001, 2050)/
     1      3.388E-05,  3.310E-05,  3.195E-05,  3.196E-05,  3.080E-05,
     2      3.004E-05,  3.004E-05,  2.927E-05,  2.927E-05,  3.044E-05,
     3      3.081E-05,  3.198E-05,  3.197E-05,  3.237E-05,  3.200E-05,
     4      3.236E-05,  3.277E-05,  3.239E-05,  3.238E-05,  3.276E-05,
     5      3.199E-05,  3.161E-05,  3.200E-05,  3.199E-05,  3.236E-05,
     6      3.352E-05,  3.351E-05,  3.429E-05,  3.388E-05,  3.350E-05,
     7      3.426E-05,  3.541E-05,  3.618E-05,  3.697E-05,  3.697E-05,
     8      3.696E-05,  3.737E-05,  3.737E-05,  3.736E-05,  3.854E-05,
     9      3.893E-05,  3.738E-05,  3.624E-05,  3.395E-05,  3.318E-05,
     $      3.550E-05,  3.819E-05,  4.092E-05,  4.204E-05,  4.245E-05/
      DATA (Y(I),I= 2051, 2100)/
     1      4.284E-05,  4.325E-05,  4.324E-05,  4.363E-05,  4.403E-05,
     2      4.442E-05,  4.365E-05,  4.284E-05,  4.246E-05,  4.244E-05,
     3      4.360E-05,  4.552E-05,  4.629E-05,  4.548E-05,  4.431E-05,
     4      4.316E-05,  4.238E-05,  4.353E-05,  4.428E-05,  4.504E-05,
     5      4.581E-05,  4.621E-05,  4.543E-05,  4.583E-05,  4.620E-05,
     6      4.659E-05,  4.622E-05,  4.661E-05,  4.621E-05,  4.583E-05,
     7      4.583E-05,  4.503E-05,  4.582E-05,  4.541E-05,  4.578E-05,
     8      4.577E-05,  4.655E-05,  4.691E-05,  4.769E-05,  4.693E-05,
     9      4.693E-05,  4.539E-05,  4.462E-05,  4.462E-05,  4.425E-05,
     $      4.351E-05,  4.391E-05,  4.468E-05,  4.588E-05,  4.741E-05/
      DATA (Y(I),I= 2101, 2150)/
     1      4.818E-05,  5.011E-05,  5.046E-05,  5.198E-05,  5.124E-05,
     2      5.086E-05,  4.931E-05,  4.776E-05,  4.700E-05,  4.623E-05,
     3      4.469E-05,  4.472E-05,  4.625E-05,  4.856E-05,  5.087E-05,
     4      5.239E-05,  5.318E-05,  5.474E-05,  5.547E-05,  5.626E-05,
     5      5.706E-05,  5.668E-05,  5.704E-05,  5.629E-05,  5.361E-05,
     6      5.208E-05,  5.017E-05,  5.209E-05,  5.554E-05,  5.864E-05,
     7      6.132E-05,  6.251E-05,  6.253E-05,  6.293E-05,  6.371E-05,
     8      6.489E-05,  6.606E-05,  6.608E-05,  6.644E-05,  6.453E-05,
     9      6.336E-05,  6.178E-05,  6.177E-05,  6.327E-05,  6.557E-05,
     $      6.711E-05,  6.864E-05,  6.901E-05,  6.937E-05,  6.975E-05/
      DATA (Y(I),I= 2151, 2200)/
     1      6.821E-05,  6.704E-05,  6.554E-05,  6.438E-05,  6.398E-05,
     2      6.318E-05,  6.319E-05,  6.280E-05,  6.048E-05,  5.933E-05,
     3      5.817E-05,  5.660E-05,  5.620E-05,  5.619E-05,  5.616E-05,
     4      5.501E-05,  5.226E-05,  4.917E-05,  4.647E-05,  4.528E-05,
     5      4.565E-05,  4.602E-05,  4.641E-05,  4.524E-05,  4.213E-05,
     6      3.943E-05,  3.713E-05,  3.711E-05,  3.866E-05,  3.982E-05,
     7      4.137E-05,  3.944E-05,  3.674E-05,  3.368E-05,  3.139E-05,
     8      3.100E-05,  3.139E-05,  3.212E-05,  3.251E-05,  3.135E-05,
     9      2.904E-05,  2.633E-05,  2.441E-05,  2.478E-05,  2.480E-05,
     $      2.557E-05,  2.594E-05,  2.477E-05,  2.404E-05,  2.404E-05/
      DATA (Y(I),I= 2201, 2250)/
     1      2.325E-05,  2.365E-05,  2.286E-05,  2.365E-05,  2.442E-05,
     2      2.404E-05,  2.406E-05,  2.401E-05,  2.403E-05,  2.480E-05,
     3      2.480E-05,  2.598E-05,  2.675E-05,  2.675E-05,  2.793E-05,
     4      2.793E-05,  2.793E-05,  2.796E-05,  2.758E-05,  2.761E-05,
     5      2.721E-05,  2.875E-05,  2.991E-05,  3.147E-05,  3.301E-05,
     6      3.266E-05,  3.187E-05,  3.148E-05,  3.110E-05,  3.261E-05,
     7      3.417E-05,  3.532E-05,  3.686E-05,  3.684E-05,  3.606E-05,
     8      3.568E-05,  3.493E-05,  3.645E-05,  3.797E-05,  4.030E-05,
     9      4.144E-05,  4.221E-05,  4.181E-05,  4.184E-05,  4.142E-05,
     $      4.181E-05,  4.257E-05,  4.372E-05,  4.412E-05,  4.334E-05/
      DATA (Y(I),I= 2251, 2300)/
     1      4.220E-05,  4.146E-05,  3.991E-05,  4.067E-05,  4.223E-05,
     2      4.377E-05,  4.608E-05,  4.609E-05,  4.608E-05,  4.571E-05,
     3      4.569E-05,  4.532E-05,  4.608E-05,  4.648E-05,  4.684E-05,
     4      4.764E-05,  4.725E-05,  4.803E-05,  4.805E-05,  4.726E-05,
     5      4.806E-05,  4.806E-05,  4.810E-05,  4.770E-05,  4.731E-05,
     6      4.617E-05,  4.577E-05,  4.424E-05,  4.464E-05,  4.502E-05,
     7      4.577E-05,  4.538E-05,  4.619E-05,  4.696E-05,  4.696E-05,
     8      4.696E-05,  4.655E-05,  4.657E-05,  4.540E-05,  4.537E-05,
     9      4.618E-05,  4.770E-05,  4.847E-05,  5.041E-05,  5.041E-05,
     $      4.928E-05,  4.773E-05,  4.619E-05,  4.502E-05,  4.582E-05/
      DATA (Y(I),I= 2301, 2350)/
     1      4.658E-05,  4.581E-05,  4.658E-05,  4.734E-05,  4.813E-05,
     2      4.809E-05,  4.927E-05,  4.848E-05,  4.814E-05,  4.813E-05,
     3      4.736E-05,  4.931E-05,  4.968E-05,  5.124E-05,  5.317E-05,
     4      5.242E-05,  5.087E-05,  5.013E-05,  4.899E-05,  4.819E-05,
     5      5.129E-05,  5.321E-05,  5.514E-05,  5.667E-05,  5.666E-05,
     6      5.669E-05,  5.745E-05,  5.708E-05,  5.744E-05,  5.746E-05,
     7      5.823E-05,  5.862E-05,  5.745E-05,  5.706E-05,  5.593E-05,
     8      5.513E-05,  5.512E-05,  5.512E-05,  5.586E-05,  5.663E-05,
     9      5.701E-05,  5.662E-05,  5.619E-05,  5.657E-05,  5.618E-05,
     $      5.502E-05,  5.385E-05,  5.307E-05,  5.191E-05,  5.113E-05/
      DATA (Y(I),I= 2351, 2400)/
     1      5.074E-05,  5.033E-05,  5.036E-05,  4.994E-05,  4.801E-05,
     2      4.723E-05,  4.609E-05,  4.531E-05,  4.724E-05,  4.992E-05,
     3      5.185E-05,  5.298E-05,  5.028E-05,  4.607E-05,  4.225E-05,
     4      3.800E-05,  3.608E-05,  3.610E-05,  3.454E-05,  3.455E-05,
     5      3.339E-05,  3.302E-05,  3.263E-05,  3.264E-05,  3.186E-05,
     6      3.262E-05,  3.417E-05,  3.493E-05,  3.647E-05,  3.530E-05,
     7      3.375E-05,  3.221E-05,  3.068E-05,  3.031E-05,  2.992E-05,
     8      3.029E-05,  3.029E-05,  3.108E-05,  3.108E-05,  3.031E-05,
     9      3.031E-05,  3.031E-05,  3.031E-05,  3.031E-05,  3.070E-05,
     $      2.994E-05,  3.111E-05,  3.109E-05,  3.184E-05,  3.300E-05/
      DATA (Y(I),I= 2401, 2450)/
     1      3.377E-05,  3.379E-05,  3.378E-05,  3.378E-05,  3.380E-05,
     2      3.456E-05,  3.535E-05,  3.572E-05,  3.651E-05,  3.650E-05,
     3      3.766E-05,  3.766E-05,  3.843E-05,  3.959E-05,  3.961E-05,
     4      4.037E-05,  4.039E-05,  4.038E-05,  4.077E-05,  3.923E-05,
     5      3.806E-05,  3.654E-05,  3.537E-05,  3.690E-05,  3.804E-05,
     6      3.956E-05,  4.189E-05,  4.264E-05,  4.228E-05,  4.150E-05,
     7      4.113E-05,  4.152E-05,  4.189E-05,  4.189E-05,  4.265E-05,
     8      4.304E-05,  4.227E-05,  4.189E-05,  4.035E-05,  3.921E-05,
     9      3.844E-05,  3.920E-05,  3.996E-05,  4.072E-05,  4.265E-05,
     $      4.265E-05,  4.265E-05,  4.189E-05,  4.148E-05,  4.151E-05/
      DATA (Y(I),I= 2451, 2500)/
     1      4.151E-05,  4.188E-05,  4.188E-05,  4.188E-05,  4.306E-05,
     2      4.226E-05,  4.304E-05,  4.385E-05,  4.357E-05,  4.342E-05,
     3      4.250E-05,  4.239E-05,  4.224E-05,  4.160E-05,  4.217E-05,
     4      4.281E-05,  4.298E-05,  4.359E-05,  4.302E-05,  4.287E-05,
     5      4.195E-05,  4.139E-05,  4.195E-05,  4.233E-05,  4.270E-05,
     6      4.229E-05,  4.268E-05,  4.234E-05,  4.191E-05,  4.138E-05,
     7      4.084E-05,  4.065E-05,  4.062E-05,  4.052E-05,  4.045E-05,
     8      4.053E-05,  4.100E-05,  4.151E-05,  4.206E-05,  4.260E-05,
     9      4.292E-05,  4.304E-05,  4.320E-05,  4.333E-05,  4.345E-05,
     $      4.377E-05,  4.404E-05,  4.432E-05,  4.460E-05,  4.453E-05/
      DATA (Y(I),I= 2501, 2550)/
     1      4.380E-05,  4.316E-05,  4.239E-05,  4.170E-05,  4.214E-05,
     2      4.302E-05,  4.387E-05,  4.472E-05,  4.523E-05,  4.454E-05,
     3      4.369E-05,  4.300E-05,  4.215E-05,  4.239E-05,  4.346E-05,
     4      4.450E-05,  4.554E-05,  4.650E-05,  4.638E-05,  4.595E-05,
     5      4.568E-05,  4.529E-05,  4.501E-05,  4.520E-05,  4.527E-05,
     6      4.534E-05,  4.541E-05,  4.559E-05,  4.584E-05,  4.613E-05,
     7      4.626E-05,  4.628E-05,  4.530E-05,  4.413E-05,  4.292E-05,
     8      4.171E-05,  4.101E-05,  4.133E-05,  4.153E-05,  4.180E-05,
     9      4.208E-05,  4.158E-05,  4.078E-05,  4.003E-05,  3.923E-05,
     $      3.854E-05,  3.888E-05,  3.936E-05,  3.993E-05,  4.038E-05/
      DATA (Y(I),I= 2551, 2600)/
     1      4.018E-05,  3.852E-05,  3.686E-05,  3.515E-05,  3.349E-05,
     2      3.280E-05,  3.284E-05,  3.285E-05,  3.277E-05,  3.273E-05,
     3      3.205E-05,  3.128E-05,  3.044E-05,  2.968E-05,  2.902E-05,
     4      2.895E-05,  2.891E-05,  2.895E-05,  2.895E-05,  2.876E-05,
     5      2.837E-05,  2.803E-05,  2.760E-05,  2.726E-05,  2.707E-05,
     6      2.707E-05,  2.700E-05,  2.700E-05,  2.689E-05,  2.685E-05,
     7      2.681E-05,  2.674E-05,  2.667E-05,  2.651E-05,  2.628E-05,
     8      2.602E-05,  2.575E-05,  2.552E-05,  2.544E-05,  2.567E-05,
     9      2.591E-05,  2.614E-05,  2.633E-05,  2.610E-05,  2.557E-05,
     $      2.511E-05,  2.457E-05,  2.416E-05,  2.450E-05,  2.512E-05/
      DATA (Y(I),I= 2601, 2650)/
     1      2.578E-05,  2.635E-05,  2.678E-05,  2.655E-05,  2.628E-05,
     2      2.594E-05,  2.563E-05,  2.548E-05,  2.606E-05,  2.667E-05,
     3      2.725E-05,  2.779E-05,  2.802E-05,  2.768E-05,  2.737E-05,
     4      2.699E-05,  2.661E-05,  2.681E-05,  2.746E-05,  2.816E-05,
     5      2.877E-05,  2.939E-05,  2.951E-05,  2.935E-05,  2.924E-05,
     6      2.913E-05,  2.897E-05,  2.916E-05,  2.936E-05,  2.951E-05,
     7      2.978E-05,  2.990E-05,  2.997E-05,  3.001E-05,  3.004E-05,
     8      3.004E-05,  3.000E-05,  2.981E-05,  2.957E-05,  2.938E-05,
     9      2.918E-05,  2.903E-05,  2.918E-05,  2.930E-05,  2.937E-05,
     $      2.948E-05,  2.944E-05,  2.898E-05,  2.848E-05,  2.806E-05/
      DATA (Y(I),I= 2651, 2700)/
     1      2.756E-05,  2.733E-05,  2.744E-05,  2.751E-05,  2.763E-05,
     2      2.767E-05,  2.751E-05,  2.716E-05,  2.682E-05,  2.635E-05,
     3      2.605E-05,  2.585E-05,  2.581E-05,  2.578E-05,  2.570E-05,
     4      2.574E-05,  2.566E-05,  2.554E-05,  2.551E-05,  2.543E-05,
     5      2.535E-05,  2.535E-05,  2.544E-05,  2.548E-05,  2.560E-05,
     6      2.564E-05,  2.569E-05,  2.585E-05,  2.589E-05,  2.597E-05,
     7      2.601E-05,  2.617E-05,  2.641E-05,  2.656E-05,  2.671E-05,
     8      2.695E-05,  2.710E-05,  2.730E-05,  2.746E-05,  2.769E-05,
     9      2.785E-05,  2.800E-05,  2.812E-05,  2.831E-05,  2.842E-05,
     $      2.858E-05,  2.878E-05,  2.889E-05,  2.908E-05,  2.920E-05/
      DATA (Y(I),I= 2701, 2750)/
     1      2.939E-05,  2.935E-05,  2.935E-05,  2.931E-05,  2.923E-05,
     2      2.927E-05,  2.919E-05,  2.915E-05,  2.915E-05,  2.911E-05,
     3      2.903E-05,  2.895E-05,  2.864E-05,  2.837E-05,  2.814E-05,
     4      2.790E-05,  2.763E-05,  2.740E-05,  2.716E-05,  2.693E-05,
     5      2.666E-05,  2.639E-05,  2.612E-05,  2.584E-05,  2.550E-05,
     6      2.526E-05,  2.499E-05,  2.464E-05,  2.437E-05,  2.414E-05,
     7      2.391E-05,  2.360E-05,  2.336E-05,  2.313E-05,  2.293E-05,
     8      2.266E-05,  2.247E-05,  2.227E-05,  2.196E-05,  2.177E-05,
     9      2.157E-05,  2.130E-05,  2.118E-05,  2.107E-05,  2.095E-05,
     $      2.075E-05,  2.063E-05,  2.044E-05,  2.036E-05,  2.017E-05/
      DATA (Y(I),I= 2751, 2800)/
     1      2.005E-05,  1.994E-05,  1.982E-05,  1.971E-05,  1.959E-05,
     2      1.944E-05,  1.933E-05,  1.925E-05,  1.906E-05,  1.898E-05,
     3      1.887E-05,  1.876E-05,  1.864E-05,  1.856E-05,  1.849E-05,
     4      1.830E-05,  1.827E-05,  1.815E-05,  1.807E-05,  1.800E-05,
     5      1.788E-05,  1.778E-05,  1.770E-05,  1.766E-05,  1.762E-05,
     6      1.755E-05,  1.755E-05,  1.751E-05,  1.740E-05,  1.736E-05,
     7      1.736E-05,  1.728E-05,  1.725E-05,  1.721E-05,  1.717E-05,
     8      1.712E-05,  1.712E-05,  1.709E-05,  1.704E-05,  1.693E-05,
     9      1.693E-05,  1.685E-05,  1.685E-05,  1.680E-05,  1.680E-05,
     $      1.672E-05,  1.672E-05,  1.672E-05,  1.672E-05,  1.668E-05/
      DATA (Y(I),I= 2801, 2850)/
     1      1.660E-05,  1.659E-05,  1.659E-05,  1.651E-05,  1.655E-05,
     2      1.655E-05,  1.655E-05,  1.651E-05,  1.647E-05,  1.647E-05,
     3      1.647E-05,  1.643E-05,  1.643E-05,  1.636E-05,  1.636E-05,
     4      1.640E-05,  1.640E-05,  1.636E-05,  1.640E-05,  1.640E-05,
     5      1.640E-05,  1.637E-05,  1.637E-05,  1.641E-05,  1.637E-05,
     6      1.641E-05,  1.638E-05,  1.626E-05,  1.622E-05,  1.619E-05,
     7      1.611E-05,  1.611E-05,  1.607E-05,  1.604E-05,  1.600E-05,
     8      1.596E-05,  1.589E-05,  1.573E-05,  1.562E-05,  1.550E-05,
     9      1.538E-05,  1.531E-05,  1.519E-05,  1.499E-05,  1.496E-05,
     $      1.480E-05,  1.469E-05,  1.461E-05,  1.445E-05,  1.437E-05/
      DATA (Y(I),I= 2851, 2900)/
     1      1.422E-05,  1.418E-05,  1.403E-05,  1.399E-05,  1.383E-05,
     2      1.375E-05,  1.360E-05,  1.356E-05,  1.356E-05,  1.360E-05,
     3      1.360E-05,  1.364E-05,  1.360E-05,  1.360E-05,  1.360E-05,
     4      1.360E-05,  1.356E-05,  1.364E-05,  1.364E-05,  1.368E-05,
     5      1.380E-05,  1.384E-05,  1.388E-05,  1.392E-05,  1.403E-05,
     6      1.403E-05,  1.411E-05,  1.411E-05,  1.419E-05,  1.423E-05,
     7      1.431E-05,  1.427E-05,  1.435E-05,  1.439E-05,  1.439E-05,
     8      1.443E-05,  1.451E-05,  1.454E-05,  1.454E-05,  1.458E-05,
     9      1.458E-05,  1.450E-05,  1.450E-05,  1.446E-05,  1.446E-05,
     $      1.446E-05,  1.438E-05,  1.442E-05,  1.441E-05,  1.433E-05/
      DATA (Y(I),I= 2901, 2950)/
     1      1.429E-05,  1.422E-05,  1.418E-05,  1.406E-05,  1.398E-05,
     2      1.387E-05,  1.379E-05,  1.367E-05,  1.360E-05,  1.352E-05,
     3      1.340E-05,  1.329E-05,  1.313E-05,  1.310E-05,  1.290E-05,
     4      1.275E-05,  1.263E-05,  1.248E-05,  1.233E-05,  1.225E-05,
     5      1.210E-05,  1.198E-05,  1.187E-05,  1.171E-05,  1.156E-05,
     6      1.144E-05,  1.133E-05,  1.118E-05,  1.106E-05,  1.098E-05,
     7      1.081E-05,  1.073E-05,  1.055E-05,  1.046E-05,  1.035E-05,
     8      1.025E-05,  1.013E-05,  1.003E-05,  9.912E-06,  9.800E-06,
     9      9.691E-06,  9.582E-06,  9.470E-06,  9.362E-06,  9.253E-06,
     $      9.175E-06,  9.097E-06,  9.027E-06,  8.950E-06,  8.883E-06/
      DATA (Y(I),I= 2951, 3000)/
     1      8.802E-06,  8.732E-06,  8.658E-06,  8.588E-06,  8.510E-06,
     2      8.440E-06,  8.385E-06,  8.342E-06,  8.307E-06,  8.267E-06,
     3      8.228E-06,  8.189E-06,  8.149E-06,  8.106E-06,  8.063E-06,
     4      8.024E-06,  7.984E-06,  7.957E-06,  7.949E-06,  7.944E-06,
     5      7.928E-06,  7.920E-06,  7.908E-06,  7.907E-06,  7.895E-06,
     6      7.883E-06,  7.874E-06,  7.870E-06,  7.854E-06,  7.858E-06,
     7      7.853E-06,  7.850E-06,  7.853E-06,  7.849E-06,  7.853E-06,
     8      7.849E-06,  7.845E-06,  7.852E-06,  7.848E-06,  7.844E-06,
     9      7.836E-06,  7.825E-06,  7.810E-06,  7.798E-06,  7.787E-06,
     $      7.776E-06,  7.753E-06,  7.741E-06,  7.726E-06,  7.715E-06/
      DATA (Y(I),I= 3001, 3050)/
     1      7.703E-06,  7.684E-06,  7.642E-06,  7.604E-06,  7.565E-06,
     2      7.527E-06,  7.492E-06,  7.454E-06,  7.416E-06,  7.370E-06,
     3      7.332E-06,  7.297E-06,  7.251E-06,  7.197E-06,  7.132E-06,
     4      7.070E-06,  7.005E-06,  6.935E-06,  6.870E-06,  6.805E-06,
     5      6.743E-06,  6.674E-06,  6.608E-06,  6.543E-06,  6.481E-06,
     6      6.393E-06,  6.312E-06,  6.223E-06,  6.143E-06,  6.058E-06,
     7      5.977E-06,  5.897E-06,  5.816E-06,  5.723E-06,  5.643E-06,
     8      5.558E-06,  5.493E-06,  5.420E-06,  5.358E-06,  5.285E-06,
     9      5.212E-06,  5.150E-06,  5.073E-06,  5.004E-06,  4.939E-06,
     $      4.866E-06,  4.804E-06,  4.743E-06,  4.731E-06,  4.727E-06/
      DATA (Y(I),I= 3051, 3100)/
     1      4.723E-06,  4.719E-06,  4.715E-06,  4.707E-06,  4.699E-06,
     2      4.699E-06,  4.695E-06,  4.691E-06,  4.680E-06,  4.695E-06,
     3      4.745E-06,  4.798E-06,  4.844E-06,  4.890E-06,  4.944E-06,
     4      4.986E-06,  5.040E-06,  5.086E-06,  5.136E-06,  5.189E-06,
     5      5.236E-06,  5.262E-06,  5.270E-06,  5.289E-06,  5.297E-06,
     6      5.305E-06,  5.324E-06,  5.332E-06,  5.347E-06,  5.355E-06,
     7      5.370E-06,  5.381E-06,  5.389E-06,  5.347E-06,  5.297E-06,
     8      5.231E-06,  5.181E-06,  5.127E-06,  5.065E-06,  5.011E-06,
     9      4.950E-06,  4.896E-06,  4.838E-06,  4.780E-06,  4.726E-06,
     $      4.661E-06,  4.595E-06,  4.525E-06,  4.456E-06,  4.390E-06/
      DATA (Y(I),I= 3101, 3150)/
     1      4.328E-06,  4.255E-06,  4.193E-06,  4.127E-06,  4.058E-06,
     2      3.988E-06,  3.926E-06,  3.849E-06,  3.760E-06,  3.679E-06,
     3      3.598E-06,  3.520E-06,  3.439E-06,  3.357E-06,  3.272E-06,
     4      3.191E-06,  3.118E-06,  3.033E-06,  2.944E-06,  2.824E-06,
     5      2.693E-06,  2.561E-06,  2.434E-06,  2.306E-06,  2.179E-06,
     6      2.043E-06,  1.917E-06,  1.790E-06,  1.663E-06,  1.533E-06,
     7      1.300E-06,  1.100E-06,  0.910E-06,  0.725E-06,  0.575E-06,
     8      0.450E-06,  0.370E-06,  0.295E-06,  0.225E-06,  0.160E-06,
     9      0.100E-06,  0.055E-06,  0.035E-06,  0.025E-06,  0.020E-06,
     *      0.015E-06,  0.010E-06,  0.006E-06,  0.003E-06,  0.00000  /
C
      DATA (Z(I),I=    1,   50)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 51, 100)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  101,  150)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 151, 200)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  201,  250)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 251, 300)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  301,  350)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 351, 400)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  401,  450)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 451, 500)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  501,  550)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 551, 600)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  601,  650)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 651, 700)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  701,  750)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I= 751, 800)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     6      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     7      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     8      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     9      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     $      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
      DATA (Z(I),I=  801,  850)/
     1      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     2      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     3      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     4      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,
     5     -5.000e-10, -1.000e-09, -1.500e-09, -2.500e-09, -5.000e-09,
     6     -7.500e-09, -1.000e-08, -1.375e-08, -1.750e-08, -2.325e-08,
     7     -3.200e-08, -4.100e-08, -5.000e-08, -6.000e-08, -7.000e-08,
     8     -8.000e-08, -9.500e-08, -1.200e-07, -2.400e-07, -3.300e-07,
     9     -3.400e-07, -3.476E-07, -3.507E-07, -3.498E-07, -3.441E-07,
     $     -3.405E-07, -3.384E-07, -3.439E-07, -3.538E-07, -3.699E-07/
      DATA (Z(I),I= 851, 900)/
     1     -3.651E-07, -3.462E-07, -3.292E-07, -3.177E-07, -3.239E-07,
     2     -3.321E-07, -3.447E-07, -3.404E-07, -3.284E-07, -3.301E-07,
     3     -3.295E-07, -3.214E-07, -3.194E-07, -3.230E-07, -3.202E-07,
     4     -3.182E-07, -3.245E-07, -3.272E-07, -3.181E-07, -3.133E-07,
     5     -3.115E-07, -3.134E-07, -3.144E-07, -3.226E-07, -3.292E-07,
     6     -3.239E-07, -3.248E-07, -3.088E-07, -3.269E-07, -3.166E-07,
     7     -3.283E-07, -3.143E-07, -3.064E-07, -3.140E-07, -2.938E-07,
     8     -2.938E-07, -3.014E-07, -2.873E-07, -2.809E-07, -2.809E-07,
     9     -3.092E-07, -2.748E-07, -2.809E-07, -2.744E-07, -2.809E-07,
     $     -2.669E-07, -2.666E-07, -2.870E-07, -2.806E-07, -2.806E-07/
      DATA (Z(I),I=  901,  950)/
     1     -2.806E-07, -2.870E-07, -2.730E-07, -2.590E-07, -2.512E-07,
     2     -2.792E-07, -2.856E-07, -2.860E-07, -2.781E-07, -2.565E-07,
     3     -2.425E-07, -2.626E-07, -2.767E-07, -2.624E-07, -2.484E-07,
     4     -2.343E-07, -2.484E-07, -2.689E-07, -2.268E-07, -2.329E-07,
     5     -2.470E-07, -2.193E-07, -2.254E-07, -2.318E-07, -2.318E-07,
     6     -2.459E-07, -2.240E-07, -2.021E-07, -2.161E-07, -2.021E-07,
     7     -2.162E-07, -2.226E-07, -2.022E-07, -2.086E-07, -2.083E-07,
     8     -2.148E-07, -2.288E-07, -2.285E-07, -2.286E-07, -2.288E-07,
     9     -2.209E-07, -2.069E-07, -2.069E-07, -2.131E-07, -2.131E-07,
     $     -1.991E-07, -1.991E-07, -1.991E-07, -1.707E-07, -1.851E-07/
      DATA (Z(I),I= 951, 1000)/
     1     -1.772E-07, -1.775E-07, -1.912E-07, -1.976E-07, -2.117E-07,
     2     -1.974E-07, -1.772E-07, -1.837E-07, -1.693E-07, -1.693E-07,
     3     -1.769E-07, -1.693E-07, -1.632E-07, -1.489E-07, -1.352E-07,
     4     -1.492E-07, -1.495E-07, -1.556E-07, -1.492E-07, -1.352E-07,
     5     -1.352E-07, -1.492E-07, -1.492E-07, -1.349E-07, -1.285E-07,
     6     -1.147E-07, -1.287E-07, -1.150E-07, -1.147E-07, -1.007E-07,
     7     -9.455E-08, -8.022E-08, -6.618E-08, -7.990E-08, -8.016E-08,
     8     -7.992E-08, -8.779E-08, -7.986E-08, -1.015E-07, -9.567E-08,
     9     -7.373E-08, -8.135E-08, -9.567E-08, -7.520E-08, -8.955E-08,
     $     -7.521E-08, -1.033E-07, -8.133E-08, -8.895E-08, -7.493E-08/
      DATA (Z(I),I= 1001, 1050)/
     1     -5.446E-08, -6.850E-08, -9.655E-08, -6.846E-08, -9.040E-08,
     2     -7.640E-08, -7.669E-08, -6.239E-08, -5.600E-08, -7.612E-08,
     3     -6.998E-08, -4.985E-08, -5.597E-08, -7.032E-08, -7.032E-08,
     4     -7.031E-08, -6.243E-08, -5.453E-08, -5.458E-08, -4.817E-08,
     5     -5.461E-08, -4.059E-08, -6.867E-08, -4.671E-08, -5.316E-08,
     6     -3.889E-08, -3.915E-08, -3.888E-08, -3.131E-08, -1.732E-08,
     7     -1.729E-08, -3.134E-08, -1.735E-08, -9.496E-09, -9.835E-09,
     8     -1.596E-08, -3.031E-08, -3.644E-08, -4.288E-08, -2.243E-08,
     9     -1.488E-08, -8.478E-10,  2.750E-08,  3.508E-08,  3.537E-08,
     $      3.540E-08,  2.780E-08, -6.711E-09, -6.391E-10,  2.739E-08/
      DATA (Z(I),I= 1051, 1100)/
     1      3.530E-08,  4.290E-08,  4.289E-08,  5.107E-08,  4.316E-08,
     2      4.345E-08,  5.102E-08,  2.937E-08,  3.696E-08,  3.699E-08,
     3      3.695E-08,  5.099E-08,  4.486E-08,  4.480E-08,  4.479E-08,
     4      5.240E-08,  4.448E-08,  6.496E-08,  7.107E-08,  8.540E-08,
     5      5.706E-08,  2.900E-08,  4.300E-08,  7.128E-08,  9.295E-08,
     6      9.322E-08,  9.936E-08,  9.935E-08,  9.934E-08,  1.134E-07,
     7      9.146E-08,  1.195E-07,  1.055E-07,  9.177E-08,  9.207E-08,
     8      9.848E-08,  9.876E-08,  9.846E-08,  9.846E-08,  1.125E-07,
     9      1.122E-07,  1.122E-07,  1.183E-07,  1.324E-07,  1.321E-07,
     $      1.245E-07,  1.245E-07,  1.108E-07,  1.251E-07,  1.315E-07/
      DATA (Z(I),I= 1101, 1150)/
     1      1.380E-07,  1.303E-07,  1.227E-07,  1.148E-07,  1.148E-07,
     2      1.289E-07,  1.350E-07,  1.210E-07,  1.351E-07,  1.556E-07,
     3      1.559E-07,  1.479E-07,  1.681E-07,  1.465E-07,  1.591E-07,
     4      1.591E-07,  1.667E-07,  1.808E-07,  1.670E-07,  1.667E-07,
     5      1.592E-07,  1.376E-07,  1.498E-07,  1.562E-07,  1.986E-07,
     6      2.206E-07,  2.209E-07,  2.147E-07,  2.069E-07,  2.069E-07,
     7      2.130E-07,  1.990E-07,  2.051E-07,  2.268E-07,  2.552E-07,
     8      2.689E-07,  2.973E-07,  3.034E-07,  3.177E-07,  3.101E-07,
     9      3.242E-07,  3.102E-07,  3.038E-07,  3.053E-07,  2.909E-07,
     $      2.909E-07,  2.988E-07,  2.988E-07,  3.128E-07,  3.128E-07/
      DATA (Z(I),I= 1151, 1200)/
     1      3.193E-07,  3.257E-07,  3.257E-07,  3.461E-07,  3.321E-07,
     2      3.321E-07,  3.257E-07,  3.256E-07,  3.257E-07,  3.394E-07,
     3      3.391E-07,  3.248E-07,  3.248E-07,  3.251E-07,  3.172E-07,
     4      3.175E-07,  3.239E-07,  3.236E-07,  3.028E-07,  3.229E-07,
     5      3.306E-07,  3.309E-07,  3.528E-07,  3.531E-07,  3.394E-07,
     6      3.315E-07,  3.175E-07,  3.236E-07,  3.017E-07,  3.157E-07,
     7      3.157E-07,  3.295E-07,  3.154E-07,  3.154E-07,  3.093E-07,
     8      2.952E-07,  2.812E-07,  2.671E-07,  2.672E-07,  2.734E-07,
     9      2.731E-07,  2.792E-07,  2.997E-07,  2.933E-07,  2.792E-07,
     $      2.792E-07,  2.653E-07,  2.372E-07,  2.289E-07,  2.290E-07/
      DATA (Z(I),I= 1201, 1250)/
     1      2.276E-07,  2.132E-07,  2.211E-07,  2.431E-07,  2.369E-07,
     2      2.448E-07,  2.388E-07,  2.247E-07,  2.043E-07,  2.043E-07,
     3      1.903E-07,  1.967E-07,  2.251E-07,  2.235E-07,  2.300E-07,
     4      2.301E-07,  2.234E-07,  2.295E-07,  2.368E-07,  2.365E-07,
     5      2.368E-07,  2.089E-07,  2.093E-07,  1.894E-07,  1.757E-07,
     6      1.895E-07,  2.035E-07,  2.319E-07,  2.394E-07,  2.236E-07,
     7      2.376E-07,  1.937E-07,  1.658E-07,  1.374E-07,  1.456E-07,
     8      1.395E-07,  1.690E-07,  1.769E-07,  1.849E-07,  1.988E-07,
     9      1.988E-07,  1.769E-07,  1.708E-07,  1.976E-07,  2.042E-07,
     $      2.105E-07,  2.171E-07,  2.310E-07,  2.249E-07,  2.389E-07/
      DATA (Z(I),I= 1251, 1300)/
     1      2.246E-07,  2.109E-07,  1.822E-07,  1.682E-07,  1.536E-07,
     2      1.467E-07,  1.529E-07,  1.593E-07,  1.521E-07,  1.866E-07,
     3      1.930E-07,  1.729E-07,  1.947E-07,  1.745E-07,  1.623E-07,
     4      1.562E-07,  1.636E-07,  1.629E-07,  1.626E-07,  1.623E-07,
     5      1.688E-07,  1.750E-07,  1.676E-07,  1.878E-07,  1.946E-07,
     6      2.151E-07,  2.086E-07,  2.086E-07,  1.946E-07,  2.087E-07,
     7      2.007E-07,  1.784E-07,  1.630E-07,  1.613E-07,  1.551E-07,
     8      1.425E-07,  1.507E-07,  1.242E-07,  1.184E-07,  1.259E-07,
     9      1.605E-07,  1.540E-07,  1.962E-07,  2.102E-07,  1.804E-07,
     $      1.585E-07,  1.225E-07,  9.271E-08,  7.834E-08,  1.067E-07/
      DATA (Z(I),I= 1301, 1350)/
     1      1.144E-07,  1.205E-07,  1.489E-07,  1.568E-07,  1.162E-07,
     2      1.039E-07,  1.399E-07,  1.196E-07,  1.263E-07,  1.267E-07,
     3      9.704E-08,  8.956E-08,  6.793E-08,  9.629E-08,  1.241E-07,
     4      1.383E-07,  1.662E-07,  1.742E-07,  1.256E-07,  7.070E-08,
     5      6.420E-08,  3.083E-08,  5.131E-08,  5.781E-08,  7.816E-08,
     6      6.300E-08,  8.356E-08,  7.602E-08,  9.796E-08,  9.801E-08,
     7      1.402E-07,  1.758E-07,  1.478E-07,  1.261E-07,  6.339E-08,
     8      2.789E-08, -2.071E-08, -2.791E-08, -4.312E-08, -3.657E-08,
     9     -9.756E-09,  3.890E-08,  2.843E-09, -3.339E-09, -3.133E-08,
     $     -5.142E-08, -7.951E-08, -5.144E-08, -4.518E-08, -9.156E-09/
      DATA (Z(I),I= 1351, 1400)/
     1      1.096E-08,  1.085E-08,  3.110E-08,  5.755E-08,  7.769E-08,
     2      6.960E-08,  6.819E-08,  4.017E-08,  2.609E-08, -8.452E-09,
     3     -2.233E-08, -4.403E-08, -1.706E-08,  1.748E-08,  2.296E-09,
     4      1.513E-08, -6.365E-09, -3.449E-08, -6.250E-08, -5.452E-08,
     5     -8.259E-08, -1.186E-07, -1.342E-07, -1.281E-07, -1.579E-07,
     6     -1.736E-07, -1.812E-07, -1.666E-07, -1.662E-07, -1.655E-07,
     7     -1.787E-07, -1.864E-07, -2.023E-07, -2.179E-07, -2.260E-07,
     8     -2.417E-07, -2.499E-07, -2.504E-07, -2.510E-07, -2.516E-07,
     9     -2.657E-07, -2.738E-07, -2.889E-07, -3.182E-07, -3.259E-07,
     $     -3.346E-07, -3.486E-07, -3.754E-07, -3.819E-07, -3.947E-07/
      DATA (Z(I),I= 1401, 1450)/
     1     -3.998E-07, -4.139E-07, -4.214E-07, -4.162E-07, -4.378E-07,
     2     -4.672E-07, -4.888E-07, -4.964E-07, -4.825E-07, -4.964E-07,
     3     -4.964E-07, -5.040E-07, -5.524E-07, -4.964E-07, -4.963E-07,
     4     -4.963E-07, -4.963E-07, -4.963E-07, -4.964E-07, -4.963E-07,
     5     -4.964E-07, -4.963E-07, -4.963E-07, -4.963E-07, -4.963E-07,
     6     -5.724E-07, -4.964E-07, -4.963E-07, -6.365E-07, -4.963E-07,
     7     -5.577E-07, -6.978E-07, -5.575E-07, -5.576E-07, -6.950E-07,
     8     -5.577E-07, -5.576E-07, -6.951E-07, -5.577E-07, -5.577E-07,
     9     -5.577E-07, -5.576E-07, -5.576E-07, -6.950E-07, -5.577E-07,
     $     -6.978E-07, -5.576E-07, -5.577E-07, -5.576E-07, -5.577E-07/
      DATA (Z(I),I= 1451, 1500)/
     1     -5.576E-07, -5.577E-07, -5.576E-07, -7.010E-07, -8.381E-07,
     2     -8.383E-07, -6.977E-07, -8.411E-07, -8.411E-07, -5.604E-07,
     3     -8.411E-07, -7.005E-07, -8.410E-07, -7.006E-07, -6.978E-07,
     4     -8.379E-07, -7.006E-07, -8.381E-07, -8.378E-07, -8.349E-07,
     5     -8.380E-07, -8.380E-07, -8.351E-07, -8.379E-07, -8.379E-07,
     6     -8.380E-07, -8.379E-07, -8.350E-07, -8.350E-07, -8.350E-07,
     7     -8.351E-07, -9.751E-07, -9.751E-07, -8.350E-07, -8.351E-07,
     8     -8.380E-07, -8.379E-07, -9.751E-07, -8.377E-07, -8.379E-07,
     9     -8.379E-07, -9.140E-07, -9.139E-07, -8.527E-07, -9.286E-07,
     $     -8.526E-07, -9.286E-07, -7.881E-07, -8.526E-07, -7.882E-07/
      DATA (Z(I),I= 1501, 1550)/
     1     -6.477E-07, -7.122E-07, -7.121E-07, -7.883E-07, -6.475E-07,
     2     -8.524E-07, -8.522E-07, -6.478E-07, -7.881E-07, -7.883E-07,
     3     -7.124E-07, -6.475E-07, -6.476E-07, -7.238E-07, -6.478E-07,
     4     -5.834E-07, -8.637E-07, -5.832E-07, -8.025E-07, -5.832E-07,
     5     -5.219E-07, -5.218E-07, -7.412E-07, -7.413E-07, -6.010E-07,
     6     -6.011E-07, -7.413E-07, -7.410E-07, -5.221E-07, -5.219E-07,
     7     -4.576E-07, -6.007E-07, -4.604E-07, -4.603E-07, -3.814E-07,
     8     -3.169E-07, -3.199E-07, -3.813E-07, -3.813E-07, -4.601E-07,
     9     -3.201E-07, -3.203E-07, -3.201E-07, -3.200E-07, -3.202E-07,
     $     -3.991E-07, -2.587E-07, -2.560E-07, -2.589E-07, -2.587E-07/
      DATA (Z(I),I= 1551, 1600)/
     1     -2.591E-07, -3.202E-07, -3.950E-08, -3.960E-08,  3.970E-08,
     2      3.965E-08, -1.009E-07, -3.950E-08,  3.965E-08, -1.009E-07,
     3      1.038E-07, -3.955E-08, -3.950E-08, -3.955E-08, -3.960E-08,
     4     -3.950E-08, -3.955E-08,  2.176E-08,  9.783E-08, -4.263E-08,
     5      3.668E-08, -1.039E-07,  3.329E-08,  3.605E-08, -2.791E-08,
     6     -2.515E-08,  1.121E-07,  1.149E-07,  1.150E-07,  1.150E-07,
     7      1.148E-07,  1.150E-07,  1.149E-07,  1.150E-07, -2.838E-08,
     8     -2.823E-08, -2.833E-08,  1.118E-07, -3.157E-08,  4.727E-08,
     9     -9.266E-08, -9.276E-08,  5.061E-08,  2.666E-07,  1.263E-07,
     $      1.906E-07,  1.262E-07, -1.403E-08, -1.544E-07, -9.308E-08/
      DATA (Z(I),I= 1601, 1650)/
     1      5.024E-08, -9.318E-08, -1.543E-07, -1.545E-07, -1.435E-08,
     2     -1.440E-08, -7.873E-08, -7.863E-08, -7.878E-08, -7.883E-08,
     3     -7.878E-08, -7.878E-08, -2.221E-07, -7.873E-08, -2.833E-07,
     4     -7.878E-08, -1.432E-07, -2.220E-07, -2.865E-07, -3.478E-07,
     5     -3.478E-07, -3.478E-07, -3.478E-07, -4.091E-07, -3.448E-07,
     6     -3.448E-07, -3.448E-07, -2.046E-07, -2.046E-07, -4.090E-07,
     7     -2.044E-07, -3.446E-07, -4.092E-07, -2.045E-07, -2.046E-07,
     8     -2.045E-07, -2.046E-07, -2.046E-07, -2.044E-07, -2.047E-07,
     9     -2.045E-07, -2.046E-07, -2.046E-07, -4.852E-07, -3.446E-07,
     $     -2.806E-07, -2.806E-07, -2.805E-07, -5.611E-07, -3.418E-07/
      DATA (Z(I),I= 1651, 1700)/
     1     -3.416E-07, -4.820E-07, -4.031E-07, -3.418E-07, -3.418E-07,
     2     -4.819E-07, -4.059E-07, -3.447E-07, -3.447E-07, -3.447E-07,
     3     -3.447E-07, -4.848E-07, -3.447E-07, -3.446E-07, -3.448E-07,
     4     -3.447E-07, -3.446E-07, -3.447E-07, -3.447E-07, -4.202E-07,
     5     -4.203E-07, -4.205E-07, -4.204E-07, -2.801E-07, -2.802E-07,
     6     -2.801E-07, -4.203E-07, -4.203E-07, -2.830E-07, -2.830E-07,
     7     -4.230E-07, -4.233E-07, -2.831E-07, -4.230E-07, -2.829E-07,
     8     -3.619E-07, -3.590E-07, -3.618E-07, -4.992E-07, -4.991E-07,
     9     -4.992E-07, -4.378E-07, -4.378E-07, -2.973E-07, -7.816E-08,
     $     -2.974E-07, -2.186E-07, -2.186E-07, -2.186E-07, -7.800E-08/
      DATA (Z(I),I= 1701, 1750)/
     1     -7.800E-08, -2.181E-07, -7.816E-08, -2.182E-07, -7.810E-08,
     2     -1.568E-07, -7.758E-08, -7.774E-08,  6.250E-08, -7.810E-08,
     3     -2.183E-07, -1.568E-07, -1.569E-07, -2.179E-07, -7.810E-08,
     4     -7.800E-08, -7.795E-08, -7.800E-08, -2.973E-07, -2.973E-07,
     5     -7.805E-08, -2.185E-07, -3.586E-07, -3.588E-07, -2.973E-07,
     6     -1.541E-07, -3.586E-07, -2.152E-07, -7.805E-08, -2.154E-07,
     7     -1.541E-07, -7.518E-08, -2.154E-07, -7.523E-08, -1.393E-07,
     8     -7.518E-08, -1.396E-07,  7.826E-10, -1.396E-07, -2.798E-07,
     9      9.391E-10, -2.798E-07, -2.799E-07, -1.396E-07, -1.396E-07,
     $     -1.425E-07, -1.397E-07, -1.397E-07, -1.430E-07, -1.427E-07/
      DATA (Z(I),I= 1751, 1800)/
     1     -1.430E-07, -2.609E-09,  4.174E-10,  5.217E-10,  3.652E-10,
     2     -2.348E-09, -1.397E-07,  3.652E-10,  5.739E-10,  4.696E-10,
     3      4.696E-10,  4.696E-10,  2.609E-10,  4.696E-10, -6.407E-08,
     4      1.406E-07,  7.654E-08, -6.402E-08, -6.402E-08, -6.412E-08,
     5     -6.402E-08,  1.205E-08, -6.396E-08, -6.407E-08, -6.402E-08,
     6     -6.407E-08, -6.402E-08, -6.402E-08, -1.252E-07,  1.487E-08,
     7     -6.130E-08, -6.130E-08, -1.256E-07, -2.046E-07, -2.045E-07,
     8     -6.130E-08, -6.130E-08, -6.438E-08, -6.449E-08, -1.288E-07,
     9     -6.438E-08, -6.449E-08, -6.449E-08, -6.443E-08, -6.443E-08,
     $     -6.438E-08, -6.438E-08, -1.258E-07,  1.435E-08,  1.440E-08/
      DATA (Z(I),I= 1801, 1850)/
     1      1.440E-08, -2.645E-08, -4.043E-08, -8.129E-08, -8.791E-08,
     2     -7.268E-08, -5.113E-08, -5.770E-08, -7.821E-08, -9.219E-08,
     3     -1.203E-07, -1.483E-07, -1.639E-07, -1.921E-07, -1.559E-07,
     4     -1.422E-07, -1.278E-07, -7.962E-08, -1.560E-08,  4.664E-08,
     5      6.057E-08,  4.649E-08,  1.883E-08,  4.696E-09, -2.332E-08,
     6     -3.751E-08, -3.751E-08, -5.937E-08, -6.553E-08, -7.983E-08,
     7     -7.983E-08, -6.584E-08, -5.932E-08, -6.548E-08, -1.079E-07,
     8     -1.436E-07, -1.717E-07, -1.576E-07, -1.782E-07, -1.719E-07,
     9     -1.502E-07, -1.489E-07, -9.939E-08, -9.814E-08, -9.819E-08,
     $     -1.110E-07, -1.174E-07, -1.313E-07, -1.454E-07, -1.594E-07/
      DATA (Z(I),I= 1851, 1900)/
     1     -1.454E-07, -1.314E-07, -1.173E-07, -1.044E-07, -1.121E-07,
     2     -1.198E-07, -1.554E-07, -1.630E-07, -1.502E-07, -1.798E-07,
     3     -1.797E-07, -1.578E-07, -1.653E-07, -1.654E-07, -1.657E-07,
     4     -1.578E-07, -1.517E-07, -1.514E-07, -1.374E-07, -1.233E-07,
     5     -1.297E-07, -1.516E-07, -1.513E-07, -1.798E-07, -1.793E-07,
     6     -1.933E-07, -1.934E-07, -2.013E-07, -1.796E-07, -1.513E-07,
     7     -1.437E-07, -1.145E-07, -1.007E-07, -9.318E-08, -9.198E-08,
     8     -9.835E-08, -1.049E-07, -9.078E-08, -9.083E-08, -9.845E-08,
     9     -1.266E-07, -1.545E-07, -2.108E-07, -2.169E-07, -1.949E-07,
     $     -2.011E-07, -1.853E-07, -1.853E-07, -1.555E-07, -1.413E-07/
      DATA (Z(I),I= 1901, 1950)/
     1     -1.333E-07, -1.474E-07, -1.478E-07, -1.619E-07, -1.834E-07,
     2     -2.255E-07, -2.536E-07, -2.536E-07, -2.472E-07, -2.677E-07,
     3     -2.255E-07, -2.055E-07, -1.555E-07, -1.495E-07, -1.480E-07,
     4     -1.684E-07, -1.748E-07, -1.955E-07, -1.880E-07, -1.945E-07,
     5     -2.007E-07, -1.926E-07, -1.989E-07, -1.910E-07, -1.831E-07,
     6     -1.469E-07, -1.393E-07, -1.454E-07, -1.658E-07, -2.003E-07,
     7     -2.352E-07, -2.133E-07, -1.981E-07, -1.844E-07, -1.705E-07,
     8     -1.780E-07, -1.932E-07, -1.930E-07, -1.724E-07, -1.521E-07,
     9     -1.176E-07, -1.393E-07, -1.676E-07, -1.812E-07, -1.955E-07,
     $     -1.672E-07, -1.536E-07, -1.252E-07, -1.113E-07, -9.715E-08/
      DATA (Z(I),I= 1951, 2000)/
     1     -9.725E-08, -6.918E-08, -1.112E-07, -9.720E-08, -1.256E-07,
     2     -1.111E-07, -8.343E-08, -8.309E-08, -5.502E-08, -2.687E-08,
     3     -2.721E-08, -2.723E-08, -2.079E-08, -1.463E-08,  5.113E-09,
     4      1.158E-08, -5.241E-08, -1.165E-07, -1.602E-07, -1.805E-07,
     5     -1.456E-07, -1.176E-07, -1.169E-07, -1.166E-07, -1.583E-07,
     6     -1.580E-07, -1.658E-07, -1.382E-07, -1.103E-07, -8.259E-08,
     7     -1.109E-07, -1.252E-07, -1.611E-07, -1.754E-07, -1.626E-07,
     8     -1.483E-07, -1.561E-07, -1.574E-07, -1.790E-07, -2.005E-07,
     9     -2.363E-07, -2.081E-07, -2.082E-07, -1.801E-07, -2.020E-07,
     $     -2.160E-07, -2.082E-07, -2.020E-07, -1.741E-07, -1.325E-07/
      DATA (Z(I),I= 2001, 2050)/
     1     -1.047E-07, -9.861E-08, -1.127E-07, -9.884E-08, -1.131E-07,
     2     -1.270E-07, -1.406E-07, -1.546E-07, -1.478E-07, -1.194E-07,
     3     -1.333E-07, -1.050E-07, -1.190E-07, -1.114E-07, -1.319E-07,
     4     -1.321E-07, -1.246E-07, -1.450E-07, -1.590E-07, -1.656E-07,
     5     -1.935E-07, -2.140E-07, -2.205E-07, -2.345E-07, -2.550E-07,
     6     -2.475E-07, -2.615E-07, -2.335E-07, -2.411E-07, -2.414E-07,
     7     -2.273E-07, -1.928E-07, -1.788E-07, -1.709E-07, -1.709E-07,
     8     -1.850E-07, -1.909E-07, -1.908E-07, -2.049E-07, -1.967E-07,
     9     -2.032E-07, -2.452E-07, -2.589E-07, -3.211E-07, -3.352E-07,
     $     -3.133E-07, -2.710E-07, -2.275E-07, -2.278E-07, -2.203E-07/
      DATA (Z(I),I= 2051, 2100)/
     1     -2.267E-07, -2.392E-07, -2.600E-07, -2.665E-07, -2.729E-07,
     2     -2.793E-07, -2.934E-07, -3.153E-07, -3.291E-07, -3.229E-07,
     3     -3.153E-07, -2.871E-07, -2.730E-07, -2.950E-07, -3.166E-07,
     4     -3.510E-07, -3.791E-07, -3.446E-07, -3.446E-07, -3.447E-07,
     5     -3.307E-07, -3.231E-07, -3.512E-07, -3.436E-07, -3.439E-07,
     6     -3.503E-07, -3.500E-07, -3.565E-07, -3.708E-07, -3.845E-07,
     7     -3.845E-07, -4.065E-07, -3.986E-07, -4.061E-07, -4.065E-07,
     8     -4.205E-07, -3.925E-07, -3.861E-07, -3.580E-07, -3.720E-07,
     9     -3.720E-07, -4.001E-07, -4.141E-07, -4.140E-07, -4.138E-07,
     $     -4.338E-07, -4.263E-07, -4.123E-07, -3.900E-07, -3.620E-07/
      DATA (Z(I),I= 2101, 2150)/
     1     -3.480E-07, -3.265E-07, -3.408E-07, -3.268E-07, -3.470E-07,
     2     -3.674E-07, -4.095E-07, -4.515E-07, -4.723E-07, -4.863E-07,
     3     -5.143E-07, -5.205E-07, -4.925E-07, -4.571E-07, -4.151E-07,
     4     -4.011E-07, -3.932E-07, -3.713E-07, -3.652E-07, -3.574E-07,
     5     -3.495E-07, -3.632E-07, -3.635E-07, -3.977E-07, -4.602E-07,
     6     -5.084E-07, -5.569E-07, -5.286E-07, -4.930E-07, -4.431E-07,
     7     -4.216E-07, -4.061E-07, -4.262E-07, -4.187E-07, -4.248E-07,
     8     -4.234E-07, -4.018E-07, -4.219E-07, -4.223E-07, -4.640E-07,
     9     -4.924E-07, -5.283E-07, -5.424E-07, -5.223E-07, -4.802E-07,
     $     -4.522E-07, -4.241E-07, -4.178E-07, -4.181E-07, -3.976E-07/
      DATA (Z(I),I= 2151, 2200)/
     1     -4.256E-07, -4.472E-07, -4.674E-07, -4.750E-07, -4.825E-07,
     2     -4.702E-07, -4.562E-07, -4.498E-07, -4.716E-07, -4.792E-07,
     3     -4.867E-07, -4.884E-07, -4.820E-07, -4.618E-07, -4.355E-07,
     4     -4.290E-07, -4.664E-07, -5.023E-07, -5.177E-07, -5.191E-07,
     5     -4.852E-07, -4.446E-07, -4.101E-07, -3.975E-07, -4.272E-07,
     6     -4.426E-07, -4.505E-07, -4.235E-07, -3.815E-07, -3.330E-07,
     7     -2.910E-07, -3.124E-07, -3.279E-07, -3.699E-07, -3.978E-07,
     8     -3.914E-07, -3.569E-07, -3.166E-07, -2.821E-07, -2.896E-07,
     9     -3.115E-07, -3.269E-07, -3.485E-07, -3.280E-07, -3.140E-07,
     $     -2.999E-07, -2.863E-07, -3.079E-07, -3.139E-07, -3.139E-07/
      DATA (Z(I),I= 2201, 2250)/
     1     -3.218E-07, -3.075E-07, -3.496E-07, -3.417E-07, -3.277E-07,
     2     -3.415E-07, -3.476E-07, -3.695E-07, -3.757E-07, -3.617E-07,
     3     -3.617E-07, -3.603E-07, -3.463E-07, -3.463E-07, -3.448E-07,
     4     -3.448E-07, -3.448E-07, -3.509E-07, -3.714E-07, -3.775E-07,
     5     -3.918E-07, -3.840E-07, -3.764E-07, -3.545E-07, -3.265E-07,
     6     -3.531E-07, -3.812E-07, -4.156E-07, -4.361E-07, -4.362E-07,
     7     -4.143E-07, -4.208E-07, -3.928E-07, -4.068E-07, -4.349E-07,
     8     -4.553E-07, -4.895E-07, -4.755E-07, -4.616E-07, -4.256E-07,
     9     -4.119E-07, -4.181E-07, -4.325E-07, -4.386E-07, -4.670E-07,
     $     -4.734E-07, -4.734E-07, -4.659E-07, -4.723E-07, -5.004E-07/
      DATA (Z(I),I= 2251, 2300)/
     1     -5.141E-07, -5.342E-07, -5.763E-07, -5.623E-07, -5.404E-07,
     2     -5.124E-07, -4.905E-07, -4.905E-07, -4.905E-07, -5.109E-07,
     3     -5.250E-07, -5.455E-07, -5.315E-07, -5.379E-07, -5.382E-07,
     4     -5.303E-07, -5.239E-07, -4.958E-07, -5.020E-07, -5.098E-07,
     5     -4.879E-07, -4.879E-07, -4.800E-07, -4.876E-07, -4.879E-07,
     6     -5.016E-07, -5.159E-07, -5.439E-07, -5.296E-07, -5.158E-07,
     7     -5.159E-07, -5.094E-07, -4.875E-07, -4.735E-07, -4.735E-07,
     8     -4.735E-07, -4.811E-07, -4.671E-07, -4.954E-07, -4.893E-07,
     9     -4.673E-07, -4.534E-07, -4.394E-07, -4.037E-07, -4.037E-07,
     $     -4.241E-07, -4.662E-07, -4.942E-07, -5.226E-07, -5.147E-07/
      DATA (Z(I),I= 2301, 2350)/
     1     -5.007E-07, -5.147E-07, -5.007E-07, -5.007E-07, -4.928E-07,
     2     -5.008E-07, -4.791E-07, -5.072E-07, -5.198E-07, -5.338E-07,
     3     -5.478E-07, -5.324E-07, -5.327E-07, -5.108E-07, -4.892E-07,
     4     -5.093E-07, -5.514E-07, -5.715E-07, -6.060E-07, -6.481E-07,
     5     -5.982E-07, -5.766E-07, -5.483E-07, -5.203E-07, -5.344E-07,
     6     -5.405E-07, -5.265E-07, -5.469E-07, -5.405E-07, -5.466E-07,
     7     -5.326E-07, -5.324E-07, -5.607E-07, -5.542E-07, -5.747E-07,
     8     -5.826E-07, -5.966E-07, -5.966E-07, -5.765E-07, -5.625E-07,
     9     -5.487E-07, -5.423E-07, -5.438E-07, -5.300E-07, -5.236E-07,
     $     -5.311E-07, -5.527E-07, -5.466E-07, -5.541E-07, -5.480E-07/
      DATA (Z(I),I= 2351, 2400)/
     1     -5.415E-07, -5.289E-07, -5.009E-07, -4.883E-07, -5.099E-07,
     2     -5.037E-07, -4.973E-07, -4.911E-07, -4.286E-07, -3.660E-07,
     3     -3.035E-07, -2.489E-07, -2.844E-07, -3.340E-07, -3.760E-07,
     4     -4.335E-07, -4.410E-07, -4.269E-07, -4.488E-07, -4.348E-07,
     5     -4.423E-07, -4.286E-07, -4.221E-07, -4.081E-07, -4.019E-07,
     6     -4.019E-07, -3.599E-07, -3.459E-07, -3.179E-07, -3.395E-07,
     7     -3.473E-07, -3.754E-07, -4.034E-07, -4.031E-07, -4.033E-07,
     8     -3.969E-07, -3.969E-07, -3.890E-07, -3.891E-07, -4.031E-07,
     9     -4.031E-07, -4.031E-07, -4.031E-07, -4.031E-07, -4.095E-07,
     $     -4.235E-07, -4.019E-07, -4.160E-07, -3.958E-07, -3.882E-07/
      DATA (Z(I),I= 2401, 2450)/
     1     -3.742E-07, -3.804E-07, -3.944E-07, -4.011E-07, -4.073E-07,
     2     -4.073E-07, -3.994E-07, -3.997E-07, -3.919E-07, -4.059E-07,
     3     -3.983E-07, -3.983E-07, -4.045E-07, -3.902E-07, -3.963E-07,
     4     -3.963E-07, -4.025E-07, -4.165E-07, -4.230E-07, -4.510E-07,
     5     -4.793E-07, -5.276E-07, -5.492E-07, -5.212E-07, -5.074E-07,
     6     -4.935E-07, -4.576E-07, -4.576E-07, -4.640E-07, -4.921E-07,
     7     -4.985E-07, -4.982E-07, -4.985E-07, -4.985E-07, -4.944E-07,
     8     -4.952E-07, -5.176E-07, -5.243E-07, -5.594E-07, -5.868E-07,
     9     -6.078E-07, -6.009E-07, -5.939E-07, -5.869E-07, -5.502E-07,
     $     -5.515E-07, -5.487E-07, -5.599E-07, -5.661E-07, -5.695E-07/
      DATA (Z(I),I= 2451, 2500)/
     1     -5.667E-07, -5.574E-07, -5.560E-07, -5.532E-07, -5.221E-07,
     2     -5.314E-07, -5.061E-07, -4.870E-07, -4.811E-07, -4.715E-07,
     3     -4.773E-07, -4.683E-07, -4.601E-07, -4.733E-07, -4.546E-07,
     4     -4.414E-07, -4.385E-07, -4.191E-07, -4.363E-07, -4.309E-07,
     5     -4.395E-07, -4.499E-07, -4.353E-07, -4.183E-07, -4.019E-07,
     6     -4.052E-07, -3.913E-07, -3.900E-07, -3.902E-07, -3.931E-07,
     7     -3.953E-07, -3.960E-07, -3.944E-07, -3.943E-07, -3.948E-07,
     8     -3.911E-07, -3.817E-07, -3.704E-07, -3.583E-07, -3.475E-07,
     9     -3.403E-07, -3.394E-07, -3.356E-07, -3.333E-07, -3.309E-07,
     $     -3.257E-07, -3.206E-07, -3.148E-07, -3.096E-07, -3.108E-07/
      DATA (Z(I),I= 2501, 2550)/
     1     -3.247E-07, -3.364E-07, -3.517E-07, -3.635E-07, -3.555E-07,
     2     -3.400E-07, -3.238E-07, -3.082E-07, -3.003E-07, -3.170E-07,
     3     -3.359E-07, -3.518E-07, -3.707E-07, -3.698E-07, -3.543E-07,
     4     -3.388E-07, -3.232E-07, -3.092E-07, -3.154E-07, -3.275E-07,
     5     -3.359E-07, -3.486E-07, -3.578E-07, -3.578E-07, -3.607E-07,
     6     -3.643E-07, -3.672E-07, -3.668E-07, -3.645E-07, -3.628E-07,
     7     -3.626E-07, -3.652E-07, -3.859E-07, -4.088E-07, -4.338E-07,
     8     -4.595E-07, -4.744E-07, -4.748E-07, -4.758E-07, -4.761E-07,
     9     -4.772E-07, -4.909E-07, -5.102E-07, -5.308E-07, -5.501E-07,
     $     -5.653E-07, -5.538E-07, -5.375E-07, -5.190E-07, -5.028E-07/
      DATA (Z(I),I= 2551, 2600)/
     1     -4.997E-07, -5.220E-07, -5.441E-07, -5.685E-07, -5.907E-07,
     2     -5.939E-07, -5.815E-07, -5.711E-07, -5.615E-07, -5.499E-07,
     3     -5.515E-07, -5.545E-07, -5.589E-07, -5.625E-07, -5.621E-07,
     4     -5.540E-07, -5.438E-07, -5.328E-07, -5.212E-07, -5.144E-07,
     5     -5.105E-07, -5.065E-07, -5.034E-07, -4.994E-07, -4.926E-07,
     6     -4.837E-07, -4.761E-07, -4.672E-07, -4.611E-07, -4.521E-07,
     7     -4.445E-07, -4.370E-07, -4.294E-07, -4.233E-07, -4.186E-07,
     8     -4.159E-07, -4.126E-07, -4.093E-07, -4.018E-07, -3.901E-07,
     9     -3.776E-07, -3.659E-07, -3.542E-07, -3.515E-07, -3.558E-07,
     $     -3.600E-07, -3.643E-07, -3.651E-07, -3.546E-07, -3.379E-07/
      DATA (Z(I),I= 2601, 2650)/
     1     -3.204E-07, -3.058E-07, -2.933E-07, -2.947E-07, -2.975E-07,
     2     -3.010E-07, -3.052E-07, -3.066E-07, -2.933E-07, -2.793E-07,
     3     -2.667E-07, -2.563E-07, -2.499E-07, -2.561E-07, -2.603E-07,
     4     -2.686E-07, -2.762E-07, -2.712E-07, -2.599E-07, -2.459E-07,
     5     -2.360E-07, -2.240E-07, -2.218E-07, -2.260E-07, -2.294E-07,
     6     -2.328E-07, -2.370E-07, -2.342E-07, -2.306E-07, -2.292E-07,
     7     -2.250E-07, -2.243E-07, -2.236E-07, -2.216E-07, -2.216E-07,
     8     -2.210E-07, -2.225E-07, -2.246E-07, -2.290E-07, -2.318E-07,
     9     -2.354E-07, -2.376E-07, -2.313E-07, -2.245E-07, -2.198E-07,
     $     -2.129E-07, -2.102E-07, -2.160E-07, -2.210E-07, -2.253E-07/
      DATA (Z(I),I= 2651, 2700)/
     1     -2.297E-07, -2.306E-07, -2.237E-07, -2.176E-07, -2.121E-07,
     2     -2.066E-07, -2.062E-07, -2.077E-07, -2.092E-07, -2.143E-07,
     3     -2.144E-07, -2.140E-07, -2.119E-07, -2.085E-07, -2.058E-07,
     4     -2.010E-07, -1.997E-07, -1.970E-07, -1.936E-07, -1.915E-07,
     5     -1.896E-07, -1.874E-07, -1.838E-07, -1.816E-07, -1.794E-07,
     6     -1.778E-07, -1.750E-07, -1.706E-07, -1.698E-07, -1.662E-07,
     7     -1.648E-07, -1.632E-07, -1.602E-07, -1.580E-07, -1.586E-07,
     8     -1.557E-07, -1.549E-07, -1.527E-07, -1.525E-07, -1.489E-07,
     9     -1.487E-07, -1.479E-07, -1.486E-07, -1.484E-07, -1.491E-07,
     $     -1.489E-07, -1.488E-07, -1.494E-07, -1.493E-07, -1.485E-07/
      DATA (Z(I),I= 2701, 2750)/
     1     -1.484E-07, -1.518E-07, -1.526E-07, -1.547E-07, -1.569E-07,
     2     -1.575E-07, -1.611E-07, -1.632E-07, -1.640E-07, -1.654E-07,
     3     -1.682E-07, -1.696E-07, -1.740E-07, -1.783E-07, -1.811E-07,
     4     -1.855E-07, -1.885E-07, -1.927E-07, -1.956E-07, -1.992E-07,
     5     -2.028E-07, -2.051E-07, -2.067E-07, -2.069E-07, -2.099E-07,
     6     -2.094E-07, -2.109E-07, -2.126E-07, -2.141E-07, -2.136E-07,
     7     -2.131E-07, -2.140E-07, -2.135E-07, -2.124E-07, -2.105E-07,
     8     -2.093E-07, -2.081E-07, -2.055E-07, -2.058E-07, -2.039E-07,
     9     -2.028E-07, -2.023E-07, -1.996E-07, -1.984E-07, -1.951E-07,
     $     -1.959E-07, -1.940E-07, -1.928E-07, -1.901E-07, -1.889E-07/
      DATA (Z(I),I= 2751, 2800)/
     1     -1.870E-07, -1.843E-07, -1.823E-07, -1.782E-07, -1.755E-07,
     2     -1.728E-07, -1.702E-07, -1.647E-07, -1.634E-07, -1.593E-07,
     3     -1.566E-07, -1.518E-07, -1.491E-07, -1.464E-07, -1.436E-07,
     4     -1.423E-07, -1.381E-07, -1.368E-07, -1.334E-07, -1.299E-07,
     5     -1.287E-07, -1.258E-07, -1.232E-07, -1.203E-07, -1.163E-07,
     6     -1.157E-07, -1.129E-07, -1.088E-07, -1.081E-07, -1.047E-07,
     7     -1.019E-07, -1.006E-07, -9.780E-08, -9.514E-08, -9.455E-08,
     8     -9.330E-08, -9.128E-08, -8.923E-08, -8.864E-08, -8.879E-08,
     9     -8.536E-08, -8.615E-08, -8.414E-08, -8.288E-08, -8.226E-08,
     $     -8.305E-08, -8.305E-08, -8.244E-08, -8.183E-08, -8.198E-08/
      DATA (Z(I),I= 2801, 2850)/
     1     -8.478E-08, -8.417E-08, -8.356E-08, -8.435E-08, -8.370E-08,
     2     -8.370E-08, -8.432E-08, -8.575E-08, -8.651E-08, -8.651E-08,
     3     -8.510E-08, -8.648E-08, -8.648E-08, -8.788E-08, -8.928E-08,
     4     -8.785E-08, -8.645E-08, -8.770E-08, -8.627E-08, -8.627E-08,
     5     -8.627E-08, -8.685E-08, -8.685E-08, -8.609E-08, -8.606E-08,
     6     -8.462E-08, -8.257E-08, -8.394E-08, -8.330E-08, -8.186E-08,
     7     -8.265E-08, -8.186E-08, -8.122E-08, -7.978E-08, -7.913E-08,
     8     -7.910E-08, -7.641E-08, -7.517E-08, -7.453E-08, -7.326E-08,
     9     -7.060E-08, -6.998E-08, -6.732E-08, -6.886E-08, -6.541E-08,
     $     -6.417E-08, -6.353E-08, -6.230E-08, -6.308E-08, -6.246E-08/
      DATA (Z(I),I= 2851, 2900)/
     1     -6.325E-08, -6.126E-08, -6.205E-08, -6.000E-08, -6.078E-08,
     2     -5.815E-08, -5.893E-08, -5.896E-08, -5.896E-08, -5.960E-08,
     3     -5.960E-08, -5.823E-08, -6.027E-08, -6.013E-08, -6.048E-08,
     4     -6.090E-08, -6.196E-08, -6.042E-08, -6.145E-08, -6.167E-08,
     5     -5.982E-08, -6.005E-08, -6.021E-08, -6.029E-08, -5.844E-08,
     6     -6.004E-08, -5.883E-08, -5.981E-08, -5.916E-08, -5.857E-08,
     7     -5.779E-08, -5.930E-08, -5.790E-08, -5.801E-08, -5.863E-08,
     8     -5.798E-08, -5.719E-08, -5.717E-08, -5.717E-08, -5.714E-08,
     9     -5.653E-08, -5.933E-08, -5.872E-08, -5.808E-08, -5.948E-08,
     $     -5.887E-08, -5.965E-08, -5.963E-08, -5.901E-08, -6.121E-08/
      DATA (Z(I),I= 2901, 2950)/
     1     -6.056E-08, -6.135E-08, -5.991E-08, -6.146E-08, -6.146E-08,
     2     -6.222E-08, -6.300E-08, -6.376E-08, -6.315E-08, -6.455E-08,
     3     -6.530E-08, -6.465E-08, -6.544E-08, -6.339E-08, -6.414E-08,
     4     -6.493E-08, -6.367E-08, -6.506E-08, -6.377E-08, -6.377E-08,
     5     -6.315E-08, -6.390E-08, -6.283E-08, -6.306E-08, -6.258E-08,
     6     -6.160E-08, -6.109E-08, -6.061E-08, -6.010E-08, -5.830E-08,
     7     -5.820E-08, -5.667E-08, -5.671E-08, -5.521E-08, -5.468E-08,
     8     -5.347E-08, -5.241E-08, -5.126E-08, -5.032E-08, -4.923E-08,
     9     -4.822E-08, -4.715E-08, -4.619E-08, -4.512E-08, -4.411E-08,
     $     -4.294E-08, -4.198E-08, -4.082E-08, -3.972E-08, -3.848E-08/
      DATA (Z(I),I= 2951, 3000)/
     1     -3.765E-08, -3.656E-08, -3.531E-08, -3.421E-08, -3.319E-08,
     2     -3.209E-08, -3.131E-08, -3.082E-08, -3.010E-08, -2.947E-08,
     3     -2.890E-08, -2.833E-08, -2.776E-08, -2.718E-08, -2.669E-08,
     4     -2.626E-08, -2.569E-08, -2.517E-08, -2.505E-08, -2.494E-08,
     5     -2.504E-08, -2.485E-08, -2.481E-08, -2.454E-08, -2.444E-08,
     6     -2.453E-08, -2.435E-08, -2.416E-08, -2.425E-08, -2.419E-08,
     7     -2.441E-08, -2.454E-08, -2.448E-08, -2.470E-08, -2.463E-08,
     8     -2.491E-08, -2.499E-08, -2.485E-08, -2.506E-08, -2.514E-08,
     9     -2.521E-08, -2.494E-08, -2.501E-08, -2.488E-08, -2.488E-08,
     $     -2.460E-08, -2.475E-08, -2.468E-08, -2.470E-08, -2.462E-08/
      DATA (Z(I),I= 3001, 3050)/
     1     -2.442E-08, -2.429E-08, -2.417E-08, -2.398E-08, -2.373E-08,
     2     -2.353E-08, -2.314E-08, -2.295E-08, -2.276E-08, -2.264E-08,
     3     -2.245E-08, -2.212E-08, -2.200E-08, -2.163E-08, -2.146E-08,
     4     -2.114E-08, -2.091E-08, -2.081E-08, -2.064E-08, -2.047E-08,
     5     -2.009E-08, -1.999E-08, -1.982E-08, -1.959E-08, -1.907E-08,
     6     -1.885E-08, -1.841E-08, -1.804E-08, -1.760E-08, -1.724E-08,
     7     -1.674E-08, -1.630E-08, -1.579E-08, -1.564E-08, -1.513E-08,
     8     -1.483E-08, -1.431E-08, -1.401E-08, -1.350E-08, -1.320E-08,
     9     -1.289E-08, -1.238E-08, -1.228E-08, -1.190E-08, -1.146E-08,
     $     -1.116E-08, -1.064E-08, -1.033E-08, -1.040E-08, -1.028E-08/
      DATA (Z(I),I= 3051, 3100)/
     1     -1.042E-08, -1.036E-08, -1.023E-08, -1.037E-08, -1.045E-08,
     2     -1.025E-08, -1.033E-08, -1.026E-08, -1.040E-08, -1.054E-08,
     3     -1.080E-08, -1.111E-08, -1.150E-08, -1.190E-08, -1.207E-08,
     4     -1.260E-08, -1.272E-08, -1.317E-08, -1.342E-08, -1.368E-08,
     5     -1.399E-08, -1.426E-08, -1.466E-08, -1.493E-08, -1.533E-08,
     6     -1.568E-08, -1.586E-08, -1.621E-08, -1.654E-08, -1.688E-08,
     7     -1.714E-08, -1.741E-08, -1.775E-08, -1.805E-08, -1.794E-08,
     8     -1.819E-08, -1.814E-08, -1.810E-08, -1.848E-08, -1.844E-08,
     9     -1.868E-08, -1.864E-08, -1.874E-08, -1.885E-08, -1.887E-08,
     $     -1.864E-08, -1.835E-08, -1.819E-08, -1.789E-08, -1.760E-08/
      DATA (Z(I),I= 3101, 3150)/
     1     -1.722E-08, -1.706E-08, -1.663E-08, -1.639E-08, -1.610E-08,
     2     -1.580E-08, -1.536E-08, -1.480E-08, -1.459E-08, -1.403E-08,
     3     -1.346E-08, -1.304E-08, -1.248E-08, -1.199E-08, -1.157E-08,
     4     -1.115E-08, -1.045E-08, -1.002E-08, -9.538E-09, -8.929E-09,
     5     -8.535E-09, -8.001E-09, -7.532E-09, -6.920E-09, -6.422E-09,
     6     -6.094E-09, -5.472E-09, -4.954E-09, -4.367E-09, -3.874E-09,
     7     -3.500E-09, -3.300E-09, -3.000E-09, -2.700E-09, -2.500E-09,
     8     -2.300E-09, -2.100E-09, -1.900E-09, -1.700E-09, -1.500E-09,
     9     -1.300E-09, -1.125E-09, -0.950E-09, -0.775E-09, -0.600E-09,
     *     -0.425E-09, -0.250E-09, -0.100E-09, -0.005E-09,  0.00000  /
C
      END   
C
C     --------------------------------------------------------------
C
      SUBROUTINE O3HHT0 (V1C,V2C,DVC,NPTC,C)                              F21120
C                                                                         F21130
      IMPLICIT REAL*8           (V)                                     ! F21140
C                                                                         F21150
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F21160
      COMMON /O3HH0/ V1S,V2S,DVS,NPTS,S(2687)                             F21170
      DIMENSION C(*)                                                      F21180
C                                                                         F21190
      DVC = DVS                                                           F21200
      V1C = V1ABS-DVC                                                     F21210
      V2C = V2ABS+DVC                                                     F21220
C                                                                         F21230
      I1 = (V1C-V1S)/DVS                                                  F21240
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F21310
         I = I1+J                                                         F21320
         C(J) = 0.                                                        F21330
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F21340
         VJ = V1C+DVC* REAL(J-1)                                          F21350
         C(J) = S(I)/VJ                                                   F21360
C                                                                         F21370
C     RADIATION FLD REMOVED FROM DIFFUSE OZONE                            F21380
C                                                                         F21390
   10 CONTINUE                                                            F21400
C                                                                         F21410
      RETURN                                                              F21420
C                                                                         F21430
      END                                                                 F21440
C
C     --------------------------------------------------------------
C
      BLOCK DATA BO3HH0                                                   F21450
C                                                                         F21460
      IMPLICIT REAL*8           (V)                                     ! F21470
C                                                                         F21480
C     O3HH0 CONTAINS O3 HARTLEY HUGGINS CROSS SECTIONS FOR 273K           F21490
C               UNITS OF (CM**2/MOL)*1.E-20                               F21500
C                                                                         F21510
C     NOW INCLUDES MOLINA & MOLINA AT 273K WITH THE TEMPERATURE           F21520
C     DEPENDENCE DETERMINED FROM THE 195K HARVARD MEASUREMENTS,           F21530
C     EMPLOYING THE BASS ALGORITHM
C
C              (CO(1+C1*(T-273.15)+C2*(T-273.15)**2);                     F21540
C
C     THIS IS ONLY FOR THE WAVELENGTH RANGE FROM .34 TO .35 MICRONS;      F21550
C     OTHERWISE, THE BASS DATA ALONE HAVE BEEN EMPLOYED BETWEEN           F21560
C     .34 AND .245 MICRONS.                                               F21570
C                                                                         F21580
C     NEW T-DEPENDENT X-SECTIONS BETWEEN .345 AND .36 MICRONS             F21590
C     HAVE NOW BEEN ADDED, BASED ON WORK BY CACCIANI, DISARRA             F21600
C     AND FIOCCO, UNIVERSITY OF ROME, 1987.  QUADRATIC TEMP               F21610
C     HAS BEEN DERIVED, AS ABOVE.                                         F21620
C                                                                         F21630
C     MOLINA & MOLINA HAVE AGAIN BEEN USED BETWEEN .245 AND .185          F21640
C     MICRONS (NO TEMPERATURE DEPENDENCE)                                 F21650
C                                                                         F21660
C     AGREEMENT AMONGST THE FOUR DATA SETS IS REASONABLE (<10%)           F21670
C     AND OFTEN EXCELLENT (0-3%)                                          F21680
C                                                                         F21690
C                                                                         F21700
      COMMON /O3HH0/ V1C,V2C,DVC,NC,                                      F21710
     *               O30001(80),O30081(80),O30161(80),O30241(80),         F21720
     *               O30321(80),O30401( 7),                               F21730
     *               C00001(80),C00081(80),C00161(80),C00241(80),         F21740
     *               C00321(80),C00401(80),C00481(80),C00561(80),         F21750
     *               C00641(80),C00721(80),C00801(80),C00881(80),         F21760
     *               C00961(80),C01041(80),C01121(80),C01201(80),         F21770
     *               C01281(80),C01361(80),C01441(80),C01521(80),         F21780
     *               C01601(80),C01681(80),C01761(80),C01841(80),         F21790
     *               C01921(80),C02001(80),C02081(80),C02161(80),         F21800
     *               C02241(40)                                           F21810
C                                                                         F21820
C     DATA V1C  /27370./,V2C  /29400./,DVC  /5./,NC  /407/ INN & TANAKA   F21830
C         DATA FROM INN & TANAKA, HANDBOOK OF GEOPHYSICS, 1957, P 16-24   F21840
C                LINEARLY INTERPOLATED BY SAC, JUNE 1985                  F21850
C                CONVERSION: (I&T)/(LOSCHMIDT 1 1987*1.2)                 F21860
C                                                                         F21870
C     DATA V1C /29405./, V2C /40800./ ,DVC /5./, NC /2280/  BASS          F21880
C         DATA FROM BASS, JUNE 1985                                       F21890
C                                                                         F21900
      DATA V1C /27370./, V2C /40800./ ,DVC /5./, NC /2687/                F21910
C                                                                         F21920
C                                                                         F21930
C    X 2.08858E-03, 1.98947E-03, 1.89037E-03, 1.79126E-03, 1.69215E-03,   F21940
C     THIS LINE OF DATA HAS BEEN REPLACED BY MONOTONICALLY INCREASING     F21950
C     VALUES                                                              F21960
C                                                                         F21970
      DATA O30001/                                                        F21980
     * 1.00000E-03, 1.15000E-03, 1.25000E-03, 1.40000E-03, 1.50000E-03,   F21990
     * 1.59304E-03, 1.62396E-03, 1.76216E-03, 1.90036E-03, 2.03856E-03,   F22000
     * 2.16538E-03, 2.02324E-03, 1.88110E-03, 1.73896E-03, 1.59682E-03,   F22010
     * 1.45468E-03, 1.31253E-03, 1.17039E-03, 1.02825E-03, 8.86108E-04,   F22020
     * 7.43963E-04, 6.01821E-04, 4.59679E-04, 5.14820E-04, 5.73044E-04,   F22030
     * 6.31269E-04, 6.89493E-04, 7.47718E-04, 8.05942E-04, 8.64167E-04,   F22040
     * 9.22392E-04, 9.80617E-04, 1.03884E-03, 1.09707E-03, 1.15528E-03,   F22050
     * 1.21351E-03, 1.27173E-03, 1.32996E-03, 1.38818E-03, 1.44641E-03,   F22060
     * 1.50463E-03, 1.56286E-03, 1.62108E-03, 1.67931E-03, 1.73753E-03,   F22070
     * 1.79575E-03, 1.85398E-03, 1.91220E-03, 1.97043E-03, 2.02865E-03,   F22080
     * 2.08688E-03, 2.14510E-03, 2.20333E-03, 2.26155E-03, 2.31978E-03,   F22090
     * 2.37800E-03, 2.43623E-03, 2.49444E-03, 2.55267E-03, 2.61089E-03,   F22100
     * 2.66912E-03, 2.72734E-03, 2.78557E-03, 2.84379E-03, 2.90202E-03,   F22110
     * 2.96024E-03, 3.01847E-03, 3.07669E-03, 3.13491E-03, 3.19313E-03,   F22120
     * 3.25136E-03, 3.30958E-03, 3.36781E-03, 3.31660E-03, 3.21583E-03,   F22130
     * 3.11505E-03, 3.22165E-03, 3.46058E-03, 3.69953E-03, 3.93846E-03/   F22140
      DATA O30081/                                                        F22150
     * 4.17739E-03, 4.41633E-03, 4.42256E-03, 4.13791E-03, 4.17894E-03,   F22160
     * 4.25583E-03, 4.33273E-03, 4.40963E-03, 4.49259E-03, 4.44532E-03,   F22170
     * 4.17540E-03, 3.84814E-03, 3.41823E-03, 3.11003E-03, 2.86548E-03,   F22180
     * 2.73912E-03, 2.70800E-03, 2.70882E-03, 2.70866E-03, 2.70816E-03,   F22190
     * 2.71228E-03, 2.78044E-03, 2.86135E-03, 3.00163E-03, 3.15222E-03,   F22200
     * 3.33394E-03, 3.48231E-03, 3.64966E-03, 3.83242E-03, 3.97733E-03,   F22210
     * 4.10299E-03, 4.26332E-03, 4.41165E-03, 4.54040E-03, 4.65544E-03,   F22220
     * 4.91897E-03, 5.23429E-03, 5.45390E-03, 5.74420E-03, 5.96314E-03,   F22230
     * 6.07198E-03, 6.07338E-03, 5.99162E-03, 5.95079E-03, 6.04655E-03,   F22240
     * 6.18239E-03, 6.56998E-03, 6.93885E-03, 7.38561E-03, 7.73029E-03,   F22250
     * 7.90493E-03, 7.72072E-03, 7.40226E-03, 6.53860E-03, 5.30328E-03,   F22260
     * 4.23000E-03, 3.45735E-03, 3.21167E-03, 3.16694E-03, 3.30966E-03,   F22270
     * 3.47431E-03, 3.68089E-03, 3.92006E-03, 4.05246E-03, 4.16408E-03,   F22280
     * 4.08710E-03, 3.98224E-03, 4.07316E-03, 4.19498E-03, 4.44990E-03,   F22290
     * 4.77881E-03, 5.08270E-03, 5.37384E-03, 5.70240E-03, 5.91906E-03,   F22300
     * 5.96745E-03, 5.92363E-03, 5.80363E-03, 5.60812E-03, 5.37450E-03/   F22310
      DATA O30161/                                                        F22320
     * 5.16202E-03, 4.98389E-03, 4.95294E-03, 5.04930E-03, 5.17576E-03,   F22330
     * 5.26042E-03, 5.22957E-03, 5.32404E-03, 5.39630E-03, 5.53353E-03,   F22340
     * 5.68057E-03, 5.78679E-03, 5.83795E-03, 5.93810E-03, 6.09330E-03,   F22350
     * 6.40001E-03, 6.69056E-03, 7.04863E-03, 7.41339E-03, 7.87421E-03,   F22360
     * 8.35570E-03, 8.97672E-03, 9.58486E-03, 1.01972E-02, 1.08463E-02,   F22370
     * 1.14105E-02, 1.18935E-02, 1.22404E-02, 1.25053E-02, 1.28759E-02,   F22380
     * 1.32169E-02, 1.37796E-02, 1.46488E-02, 1.57324E-02, 1.68897E-02,   F22390
     * 1.78560E-02, 1.87101E-02, 1.92197E-02, 1.94106E-02, 1.90711E-02,   F22400
     * 1.86585E-02, 1.82149E-02, 1.82219E-02, 1.85639E-02, 1.91924E-02,   F22410
     * 2.01342E-02, 2.12312E-02, 2.26362E-02, 2.39610E-02, 2.55156E-02,   F22420
     * 2.71338E-02, 2.87904E-02, 3.04268E-02, 3.17055E-02, 3.28248E-02,   F22430
     * 3.36026E-02, 3.36867E-02, 3.26393E-02, 2.99356E-02, 2.56607E-02,   F22440
     * 2.11545E-02, 1.79508E-02, 1.59757E-02, 1.49569E-02, 1.46214E-02,   F22450
     * 1.46214E-02, 1.48217E-02, 1.51379E-02, 1.53816E-02, 1.58087E-02,   F22460
     * 1.62186E-02, 1.66627E-02, 1.70961E-02, 1.76101E-02, 1.81759E-02,   F22470
     * 1.86154E-02, 1.88889E-02, 1.89577E-02, 1.89316E-02, 1.88826E-02/   F22480
      DATA O30241/                                                        F22490
     * 1.90915E-02, 1.95550E-02, 2.02707E-02, 2.11620E-02, 2.21844E-02,   F22500
     * 2.30920E-02, 2.37270E-02, 2.37422E-02, 2.33578E-02, 2.20358E-02,   F22510
     * 1.96239E-02, 1.73329E-02, 1.57013E-02, 1.50566E-02, 1.49248E-02,   F22520
     * 1.52044E-02, 1.57658E-02, 1.63436E-02, 1.68986E-02, 1.74180E-02,   F22530
     * 1.78192E-02, 1.80677E-02, 1.79927E-02, 1.77900E-02, 1.75599E-02,   F22540
     * 1.74982E-02, 1.76674E-02, 1.81633E-02, 1.87826E-02, 1.96898E-02,   F22550
     * 2.06898E-02, 2.17167E-02, 2.28231E-02, 2.40702E-02, 2.55084E-02,   F22560
     * 2.69701E-02, 2.86915E-02, 3.05796E-02, 3.22328E-02, 3.42637E-02,   F22570
     * 3.61708E-02, 3.79118E-02, 3.94418E-02, 4.07333E-02, 4.17158E-02,   F22580
     * 4.17081E-02, 4.01127E-02, 3.65411E-02, 3.25123E-02, 2.98737E-02,   F22590
     * 2.83616E-02, 2.79907E-02, 2.80571E-02, 2.84778E-02, 2.91698E-02,   F22600
     * 2.99500E-02, 3.07468E-02, 3.13903E-02, 3.19811E-02, 3.24616E-02,   F22610
     * 3.26503E-02, 3.26829E-02, 3.27688E-02, 3.36446E-02, 3.55133E-02,   F22620
     * 3.88447E-02, 4.28854E-02, 4.55381E-02, 4.77161E-02, 4.93567E-02,   F22630
     * 4.95127E-02, 5.00492E-02, 5.06233E-02, 5.12739E-02, 5.20327E-02,   F22640
     * 5.29001E-02, 5.38677E-02, 5.49272E-02, 5.60703E-02, 5.72886E-02/   F22650
      DATA O30321/                                                        F22660
     * 5.85739E-02, 5.99178E-02, 6.13170E-02, 6.28474E-02, 6.46499E-02,   F22670
     * 6.68672E-02, 6.96421E-02, 7.31174E-02, 7.74361E-02, 8.27413E-02,   F22680
     * 8.91756E-02, 9.67018E-02, 1.04844E-01, 1.13063E-01, 1.20818E-01,   F22690
     * 1.27567E-01, 1.32771E-01, 1.35888E-01, 1.36377E-01, 1.33780E-01,   F22700
     * 1.28385E-01, 1.20887E-01, 1.11978E-01, 1.02354E-01, 9.27108E-02,   F22710
     * 8.37418E-02, 7.61423E-02, 7.06032E-02, 6.74255E-02, 6.62092E-02,   F22720
     * 6.64813E-02, 6.77689E-02, 6.95995E-02, 7.15004E-02, 7.29991E-02,   F22730
     * 7.36229E-02, 7.29641E-02, 7.11015E-02, 6.83345E-02, 6.49638E-02,   F22740
     * 6.12897E-02, 5.76125E-02, 5.42326E-02, 5.14504E-02, 4.95645E-02,   F22750
     * 4.87078E-02, 4.87234E-02, 4.94254E-02, 5.06280E-02, 5.21454E-02,   F22760
     * 5.37919E-02, 5.53818E-02, 5.67293E-02, 5.76709E-02, 5.82319E-02,   F22770
     * 5.85334E-02, 5.86968E-02, 5.88439E-02, 5.90963E-02, 5.95756E-02,   F22780
     * 6.04035E-02, 6.17016E-02, 6.35548E-02, 6.59664E-02, 6.89282E-02,   F22790
     * 7.24326E-02, 7.64718E-02, 8.10380E-02, 8.61236E-02, 9.17211E-02,   F22800
     * 9.78192E-02, 1.04353E-01, 1.11218E-01, 1.18308E-01, 1.25519E-01,   F22810
     * 1.32745E-01, 1.39881E-01, 1.46821E-01, 1.53461E-01, 1.59687E-01/   F22820
C                                                                         F22830
C    X 1.64187E-01, 1.69368E-01, 1.74549E-01, 1.79731E-01, 1.84912E-01,   F22840
C      1.90094E-01, 1.95275E-01/                                          F22850
C   THE VALUE AT 29400. HAS BEEN CHANGED TO PROVIDE A SMOOTH TRANSITION   F22860
C    X 1.90094E-01, 1.85275E-01/                                          F22870
C                                                                         F22880
      DATA O30401/                                                        F22890
     * 1.65365E-01, 1.70353E-01, 1.74507E-01, 1.77686E-01, 1.79748E-01,   F22900
     * 1.80549E-01, 1.79948E-01/                                          F22910
C                                                                         F22920
C                                                                         F22930
C    FOLLOWING DATA ARE FROM BASS JUNE 1985                               F22940
C                                                                         F22950
      DATA C00001 /                                                       F22960
     * 1.81094E-01, 1.57760E-01, 1.37336E-01, 1.19475E-01, 1.17191E-01,   F22970
     * 1.14331E-01, 1.15984E-01, 1.10412E-01, 1.12660E-01, 1.16014E-01,   F22980
     * 1.15060E-01, 1.12041E-01, 1.11611E-01, 1.00378E-01, 9.54850E-02,   F22990
     * 9.87528E-02, 9.46153E-02, 9.53093E-02, 9.72653E-02, 9.66468E-02,   F23000
     * 9.39750E-02, 1.03552E-01, 1.01361E-01, 1.04315E-01, 1.12842E-01,   F23010
     * 1.02800E-01, 1.09576E-01, 1.05577E-01, 1.17334E-01, 1.25763E-01,   F23020
     * 1.27597E-01, 1.34267E-01, 1.44799E-01, 1.57366E-01, 1.67369E-01,   F23030
     * 1.81778E-01, 1.89207E-01, 2.01376E-01, 2.10310E-01, 2.21721E-01,   F23040
     * 2.43162E-01, 2.55542E-01, 2.75312E-01, 2.88576E-01, 3.02505E-01,   F23050
     * 3.15141E-01, 3.28908E-01, 3.49000E-01, 3.56620E-01, 3.59852E-01,   F23060
     * 3.57517E-01, 3.12924E-01, 2.63610E-01, 2.50854E-01, 2.25642E-01,   F23070
     * 2.15954E-01, 2.12099E-01, 2.13039E-01, 2.12286E-01, 2.17214E-01,   F23080
     * 2.28784E-01, 2.28276E-01, 2.34677E-01, 2.30730E-01, 2.16107E-01,   F23090
     * 1.99471E-01, 1.85629E-01, 1.72730E-01, 1.56229E-01, 1.38156E-01,   F23100
     * 1.37641E-01, 1.33169E-01, 1.32759E-01, 1.30102E-01, 1.35396E-01,   F23110
     * 1.37976E-01, 1.41571E-01, 1.46448E-01, 1.44508E-01, 1.47612E-01/   F23120
      DATA C00081 /                                                       F23130
     * 1.47424E-01, 1.48173E-01, 1.52936E-01, 1.58908E-01, 1.58808E-01,   F23140
     * 1.59860E-01, 1.73936E-01, 1.84109E-01, 1.95143E-01, 2.08267E-01,   F23150
     * 2.19256E-01, 2.31653E-01, 2.46400E-01, 2.60437E-01, 2.70792E-01,   F23160
     * 2.79749E-01, 2.91068E-01, 2.98080E-01, 3.10421E-01, 3.24540E-01,   F23170
     * 3.39003E-01, 3.58322E-01, 3.81520E-01, 4.02798E-01, 4.35972E-01,   F23180
     * 4.56220E-01, 4.79037E-01, 5.02597E-01, 5.24648E-01, 5.33964E-01,   F23190
     * 5.39211E-01, 5.43613E-01, 5.28793E-01, 4.94103E-01, 4.34481E-01,   F23200
     * 3.76792E-01, 3.37161E-01, 3.15750E-01, 3.11042E-01, 3.08745E-01,   F23210
     * 3.09195E-01, 3.05859E-01, 3.01443E-01, 2.88111E-01, 2.81303E-01,   F23220
     * 2.75329E-01, 2.60812E-01, 2.59337E-01, 2.45576E-01, 2.40470E-01,   F23230
     * 2.39705E-01, 2.45389E-01, 2.49801E-01, 2.53235E-01, 2.54387E-01,   F23240
     * 2.64311E-01, 2.74146E-01, 2.89737E-01, 2.96673E-01, 3.07337E-01,   F23250
     * 3.24380E-01, 3.42266E-01, 3.59522E-01, 3.78005E-01, 3.97178E-01,   F23260
     * 4.23351E-01, 4.45925E-01, 4.63029E-01, 4.94843E-01, 5.19418E-01,   F23270
     * 5.49928E-01, 5.69115E-01, 6.02396E-01, 6.43471E-01, 6.76401E-01,   F23280
     * 7.14024E-01, 7.42425E-01, 7.60916E-01, 7.83319E-01, 7.98299E-01/   F23290
      DATA C00161 /                                                       F23300
     * 7.76672E-01, 7.22769E-01, 6.45967E-01, 5.80850E-01, 5.76514E-01,   F23310
     * 5.79380E-01, 5.90359E-01, 6.21721E-01, 6.37540E-01, 6.52572E-01,   F23320
     * 6.63442E-01, 6.69026E-01, 6.69038E-01, 6.53319E-01, 6.21950E-01,   F23330
     * 5.47619E-01, 4.58994E-01, 4.14888E-01, 3.97736E-01, 3.88775E-01,   F23340
     * 3.87424E-01, 3.93567E-01, 4.03442E-01, 4.05217E-01, 4.12848E-01,   F23350
     * 4.12246E-01, 4.16620E-01, 4.13195E-01, 4.08467E-01, 4.13104E-01,   F23360
     * 4.24498E-01, 4.32002E-01, 4.46361E-01, 4.61131E-01, 4.77228E-01,   F23370
     * 4.96519E-01, 5.16764E-01, 5.38966E-01, 5.54187E-01, 5.73748E-01,   F23380
     * 6.07260E-01, 6.34358E-01, 6.60286E-01, 6.95533E-01, 7.37090E-01,   F23390
     * 7.83894E-01, 8.19557E-01, 8.49244E-01, 8.91832E-01, 9.44885E-01,   F23400
     * 9.86271E-01, 1.02262E+00, 1.07242E+00, 1.12162E+00, 1.18287E+00,   F23410
     * 1.22402E+00, 1.24978E+00, 1.24392E+00, 1.19668E+00, 1.11562E+00,   F23420
     * 1.03983E+00, 9.31884E-01, 8.35307E-01, 7.92620E-01, 7.81980E-01,   F23430
     * 7.89623E-01, 8.05987E-01, 8.27344E-01, 8.57514E-01, 8.66302E-01,   F23440
     * 8.72092E-01, 8.66840E-01, 8.40536E-01, 7.87360E-01, 7.35743E-01,   F23450
     * 6.92039E-01, 6.64032E-01, 6.48360E-01, 6.46288E-01, 6.49505E-01/   F23460
      DATA C00241 /                                                       F23470
     * 6.69937E-01, 6.81006E-01, 7.00969E-01, 7.19834E-01, 7.26964E-01,   F23480
     * 7.50591E-01, 7.73600E-01, 8.00673E-01, 8.20347E-01, 8.37855E-01,   F23490
     * 8.66780E-01, 9.04297E-01, 9.46300E-01, 9.69134E-01, 9.97928E-01,   F23500
     * 1.06388E+00, 1.11032E+00, 1.15221E+00, 1.21324E+00, 1.24462E+00,   F23510
     * 1.31978E+00, 1.35617E+00, 1.38792E+00, 1.39196E+00, 1.35161E+00,   F23520
     * 1.29381E+00, 1.30295E+00, 1.32965E+00, 1.37024E+00, 1.44064E+00,   F23530
     * 1.50484E+00, 1.57200E+00, 1.62097E+00, 1.67874E+00, 1.72676E+00,   F23540
     * 1.73383E+00, 1.66091E+00, 1.54936E+00, 1.35454E+00, 1.20070E+00,   F23550
     * 1.14609E+00, 1.13642E+00, 1.13784E+00, 1.14609E+00, 1.14531E+00,   F23560
     * 1.16024E+00, 1.16891E+00, 1.16111E+00, 1.14192E+00, 1.09903E+00,   F23570
     * 1.05745E+00, 1.02341E+00, 1.00121E+00, 1.00036E+00, 1.00576E+00,   F23580
     * 1.02405E+00, 1.04379E+00, 1.07623E+00, 1.11347E+00, 1.17305E+00,   F23590
     * 1.20016E+00, 1.22697E+00, 1.27479E+00, 1.32572E+00, 1.38690E+00,   F23600
     * 1.43768E+00, 1.48379E+00, 1.55317E+00, 1.64020E+00, 1.71268E+00,   F23610
     * 1.77183E+00, 1.85824E+00, 1.95131E+00, 2.04609E+00, 2.13151E+00,   F23620
     * 2.17777E+00, 2.22832E+00, 2.26886E+00, 2.19775E+00, 2.05087E+00/   F23630
      DATA C00321 /                                                       F23640
     * 1.96103E+00, 1.95554E+00, 1.98037E+00, 2.05440E+00, 2.11629E+00,   F23650
     * 2.17893E+00, 2.24384E+00, 2.30464E+00, 2.32525E+00, 2.29945E+00,   F23660
     * 2.21712E+00, 2.03430E+00, 1.82139E+00, 1.70354E+00, 1.64631E+00,   F23670
     * 1.62164E+00, 1.61356E+00, 1.63900E+00, 1.66313E+00, 1.67409E+00,   F23680
     * 1.69143E+00, 1.70181E+00, 1.69165E+00, 1.67699E+00, 1.67879E+00,   F23690
     * 1.67312E+00, 1.68133E+00, 1.70002E+00, 1.72500E+00, 1.76308E+00,   F23700
     * 1.80634E+00, 1.87548E+00, 1.94924E+00, 1.99812E+00, 2.05333E+00,   F23710
     * 2.14035E+00, 2.21847E+00, 2.27412E+00, 2.29752E+00, 2.30750E+00,   F23720
     * 2.36165E+00, 2.44394E+00, 2.52782E+00, 2.61343E+00, 2.71640E+00,   F23730
     * 2.81613E+00, 2.93679E+00, 3.01577E+00, 3.15995E+00, 3.15931E+00,   F23740
     * 2.96658E+00, 2.73295E+00, 2.67480E+00, 2.66652E+00, 2.69393E+00,   F23750
     * 2.75102E+00, 2.86503E+00, 2.99163E+00, 2.99576E+00, 3.02603E+00,   F23760
     * 2.98415E+00, 2.79309E+00, 2.65337E+00, 2.50962E+00, 2.43207E+00,   F23770
     * 2.34812E+00, 2.34872E+00, 2.35186E+00, 2.39477E+00, 2.42629E+00,   F23780
     * 2.48068E+00, 2.55087E+00, 2.55952E+00, 2.56497E+00, 2.64323E+00,   F23790
     * 2.67961E+00, 2.66263E+00, 2.70243E+00, 2.74911E+00, 2.81786E+00/   F23800
      DATA C00401 /                                                       F23810
     * 2.88684E+00, 2.97790E+00, 3.04305E+00, 3.13053E+00, 3.23857E+00,   F23820
     * 3.35582E+00, 3.40654E+00, 3.38117E+00, 3.36296E+00, 3.39480E+00,   F23830
     * 3.49066E+00, 3.60832E+00, 3.71817E+00, 3.83924E+00, 3.96355E+00,   F23840
     * 4.03656E+00, 4.00518E+00, 3.90389E+00, 3.74790E+00, 3.61385E+00,   F23850
     * 3.57066E+00, 3.59438E+00, 3.66182E+00, 3.71176E+00, 3.75255E+00,   F23860
     * 3.79101E+00, 3.85278E+00, 3.85027E+00, 3.81112E+00, 3.72553E+00,   F23870
     * 3.61017E+00, 3.54384E+00, 3.52406E+00, 3.54097E+00, 3.59375E+00,   F23880
     * 3.66312E+00, 3.72632E+00, 3.76825E+00, 3.86798E+00, 3.92916E+00,   F23890
     * 3.95610E+00, 4.00120E+00, 4.05865E+00, 4.11981E+00, 4.14634E+00,   F23900
     * 4.19109E+00, 4.20317E+00, 4.25754E+00, 4.35131E+00, 4.48573E+00,   F23910
     * 4.58716E+00, 4.67462E+00, 4.78228E+00, 4.91196E+00, 5.01871E+00,   F23920
     * 5.10663E+00, 5.17780E+00, 5.21393E+00, 5.18144E+00, 5.04379E+00,   F23930
     * 4.86504E+00, 4.78569E+00, 4.72717E+00, 4.69132E+00, 4.65797E+00,   F23940
     * 4.60305E+00, 4.59798E+00, 4.65300E+00, 4.69707E+00, 4.74790E+00,   F23950
     * 4.82581E+00, 4.80953E+00, 4.80517E+00, 4.82685E+00, 4.82321E+00,   F23960
     * 4.84806E+00, 4.88591E+00, 4.91759E+00, 4.98074E+00, 5.07071E+00/   F23970
      DATA C00481 /                                                       F23980
     * 5.18733E+00, 5.30567E+00, 5.38670E+00, 5.43942E+00, 5.51797E+00,   F23990
     * 5.62652E+00, 5.71228E+00, 5.82347E+00, 5.91434E+00, 6.00171E+00,   F24000
     * 6.06977E+00, 6.13040E+00, 6.21990E+00, 6.29980E+00, 6.37206E+00,   F24010
     * 6.48233E+00, 6.53068E+00, 6.53275E+00, 6.56858E+00, 6.54577E+00,   F24020
     * 6.50472E+00, 6.41504E+00, 6.33853E+00, 6.31184E+00, 6.21253E+00,   F24030
     * 6.22034E+00, 6.26918E+00, 6.28982E+00, 6.29461E+00, 6.35418E+00,   F24040
     * 6.40956E+00, 6.38020E+00, 6.39784E+00, 6.45383E+00, 6.50134E+00,   F24050
     * 6.56808E+00, 6.58850E+00, 6.58882E+00, 6.65097E+00, 6.75259E+00,   F24060
     * 6.83256E+00, 6.92593E+00, 6.98083E+00, 7.03632E+00, 7.11147E+00,   F24070
     * 7.15622E+00, 7.21106E+00, 7.27319E+00, 7.33382E+00, 7.38601E+00,   F24080
     * 7.48971E+00, 7.61459E+00, 7.70134E+00, 7.76194E+00, 7.85534E+00,   F24090
     * 7.99519E+00, 8.12227E+00, 8.25461E+00, 8.34670E+00, 8.42733E+00,   F24100
     * 8.51806E+00, 8.57638E+00, 8.56481E+00, 8.55461E+00, 8.55593E+00,   F24110
     * 8.58756E+00, 8.50070E+00, 8.54400E+00, 8.57575E+00, 8.62083E+00,   F24120
     * 8.60684E+00, 8.67824E+00, 8.72069E+00, 8.79127E+00, 8.85479E+00,   F24130
     * 8.86770E+00, 8.90574E+00, 8.91531E+00, 8.94800E+00, 9.00167E+00/   F24140
      DATA C00561 /                                                       F24150
     * 9.14051E+00, 9.25421E+00, 9.39694E+00, 9.50896E+00, 9.53190E+00,   F24160
     * 9.55977E+00, 9.53482E+00, 9.49662E+00, 9.53359E+00, 9.54007E+00,   F24170
     * 9.49809E+00, 9.49373E+00, 9.53282E+00, 9.63757E+00, 9.67855E+00,   F24180
     * 9.67633E+00, 9.67045E+00, 9.79481E+00, 9.93420E+00, 1.00234E+01,   F24190
     * 1.01372E+01, 1.02577E+01, 1.05056E+01, 1.07873E+01, 1.09967E+01,   F24200
     * 1.10873E+01, 1.11624E+01, 1.13006E+01, 1.14875E+01, 1.16106E+01,   F24210
     * 1.16744E+01, 1.17582E+01, 1.17709E+01, 1.18537E+01, 1.19623E+01,   F24220
     * 1.19763E+01, 1.19879E+01, 1.20384E+01, 1.20763E+01, 1.20826E+01,   F24230
     * 1.20449E+01, 1.19747E+01, 1.20227E+01, 1.21805E+01, 1.23134E+01,   F24240
     * 1.24042E+01, 1.25614E+01, 1.26828E+01, 1.26645E+01, 1.26963E+01,   F24250
     * 1.28226E+01, 1.28720E+01, 1.28981E+01, 1.29462E+01, 1.29363E+01,   F24260
     * 1.29199E+01, 1.29797E+01, 1.28860E+01, 1.29126E+01, 1.30205E+01,   F24270
     * 1.31327E+01, 1.31722E+01, 1.31901E+01, 1.33189E+01, 1.34833E+01,   F24280
     * 1.36228E+01, 1.37474E+01, 1.38548E+01, 1.39450E+01, 1.40926E+01,   F24290
     * 1.43099E+01, 1.44836E+01, 1.46257E+01, 1.47755E+01, 1.49163E+01,   F24300
     * 1.51038E+01, 1.53308E+01, 1.54194E+01, 1.54852E+01, 1.55968E+01/   F24310
      DATA C00641 /                                                       F24320
     * 1.57025E+01, 1.58667E+01, 1.60365E+01, 1.61427E+01, 1.62967E+01,   F24330
     * 1.64735E+01, 1.66123E+01, 1.67268E+01, 1.67673E+01, 1.67825E+01,   F24340
     * 1.68898E+01, 1.68178E+01, 1.68216E+01, 1.68574E+01, 1.68799E+01,   F24350
     * 1.70317E+01, 1.70767E+01, 1.71508E+01, 1.72965E+01, 1.73421E+01,   F24360
     * 1.73937E+01, 1.74420E+01, 1.74535E+01, 1.75110E+01, 1.75497E+01,   F24370
     * 1.75149E+01, 1.75955E+01, 1.78260E+01, 1.78271E+01, 1.79750E+01,   F24380
     * 1.80600E+01, 1.81597E+01, 1.83454E+01, 1.85243E+01, 1.87382E+01,   F24390
     * 1.88904E+01, 1.90395E+01, 1.92759E+01, 1.95398E+01, 1.97712E+01,   F24400
     * 1.98487E+01, 1.99522E+01, 2.02363E+01, 2.03271E+01, 2.07090E+01,   F24410
     * 2.09195E+01, 2.10974E+01, 2.11702E+01, 2.12964E+01, 2.14339E+01,   F24420
     * 2.15764E+01, 2.17351E+01, 2.18486E+01, 2.19700E+01, 2.21663E+01,   F24430
     * 2.24244E+01, 2.24813E+01, 2.25248E+01, 2.26357E+01, 2.26457E+01,   F24440
     * 2.27249E+01, 2.27172E+01, 2.27123E+01, 2.26859E+01, 2.27216E+01,   F24450
     * 2.29306E+01, 2.30711E+01, 2.31374E+01, 2.31815E+01, 2.33423E+01,   F24460
     * 2.33810E+01, 2.36430E+01, 2.36807E+01, 2.36676E+01, 2.38607E+01,   F24470
     * 2.41559E+01, 2.43413E+01, 2.44401E+01, 2.45968E+01, 2.47927E+01/   F24480
      DATA C00721 /                                                       F24490
     * 2.50743E+01, 2.53667E+01, 2.55749E+01, 2.57357E+01, 2.58927E+01,   F24500
     * 2.61523E+01, 2.64110E+01, 2.66650E+01, 2.68829E+01, 2.70635E+01,   F24510
     * 2.72797E+01, 2.75064E+01, 2.77229E+01, 2.80341E+01, 2.82003E+01,   F24520
     * 2.83346E+01, 2.83909E+01, 2.86212E+01, 2.88006E+01, 2.89577E+01,   F24530
     * 2.90965E+01, 2.91834E+01, 2.93224E+01, 2.94094E+01, 2.94848E+01,   F24540
     * 2.96584E+01, 2.96749E+01, 2.97760E+01, 2.99163E+01, 3.00238E+01,   F24550
     * 3.01290E+01, 3.02307E+01, 3.03663E+01, 3.05897E+01, 3.07937E+01,   F24560
     * 3.10403E+01, 3.11778E+01, 3.13271E+01, 3.15799E+01, 3.18435E+01,   F24570
     * 3.21614E+01, 3.25097E+01, 3.27701E+01, 3.29600E+01, 3.32583E+01,   F24580
     * 3.36348E+01, 3.40282E+01, 3.41751E+01, 3.44128E+01, 3.46199E+01,   F24590
     * 3.49363E+01, 3.52087E+01, 3.54056E+01, 3.55596E+01, 3.56694E+01,   F24600
     * 3.58104E+01, 3.60276E+01, 3.62818E+01, 3.63505E+01, 3.66069E+01,   F24610
     * 3.67544E+01, 3.70664E+01, 3.72525E+01, 3.73491E+01, 3.76006E+01,   F24620
     * 3.77102E+01, 3.78970E+01, 3.81254E+01, 3.82728E+01, 3.81720E+01,   F24630
     * 3.82781E+01, 3.84982E+01, 3.87202E+01, 3.89958E+01, 3.94148E+01,   F24640
     * 3.98434E+01, 3.98952E+01, 4.01573E+01, 4.06014E+01, 4.09651E+01/   F24650
      DATA C00801 /                                                       F24660
     * 4.12821E+01, 4.16849E+01, 4.19899E+01, 4.22719E+01, 4.27736E+01,   F24670
     * 4.32254E+01, 4.33883E+01, 4.39831E+01, 4.39414E+01, 4.42613E+01,   F24680
     * 4.46503E+01, 4.49027E+01, 4.50384E+01, 4.52929E+01, 4.57269E+01,   F24690
     * 4.56433E+01, 4.57350E+01, 4.60128E+01, 4.60487E+01, 4.61183E+01,   F24700
     * 4.64397E+01, 4.68211E+01, 4.70706E+01, 4.72821E+01, 4.74972E+01,   F24710
     * 4.78253E+01, 4.81615E+01, 4.84480E+01, 4.85703E+01, 4.87397E+01,   F24720
     * 4.90015E+01, 4.93673E+01, 4.97291E+01, 4.99836E+01, 5.02975E+01,   F24730
     * 5.05572E+01, 5.08226E+01, 5.13433E+01, 5.17112E+01, 5.19703E+01,   F24740
     * 5.23128E+01, 5.27305E+01, 5.30599E+01, 5.34555E+01, 5.39625E+01,   F24750
     * 5.43627E+01, 5.45446E+01, 5.49263E+01, 5.53511E+01, 5.57270E+01,   F24760
     * 5.60904E+01, 5.63875E+01, 5.68475E+01, 5.73172E+01, 5.81134E+01,   F24770
     * 5.86399E+01, 5.90384E+01, 5.91417E+01, 5.90883E+01, 5.93610E+01,   F24780
     * 5.95794E+01, 5.99600E+01, 5.98493E+01, 5.99441E+01, 6.02748E+01,   F24790
     * 6.04778E+01, 6.05233E+01, 6.07194E+01, 6.11589E+01, 6.13324E+01,   F24800
     * 6.17685E+01, 6.23166E+01, 6.31055E+01, 6.38211E+01, 6.42320E+01,   F24810
     * 6.45195E+01, 6.51125E+01, 6.56765E+01, 6.59286E+01, 6.62716E+01/   F24820
      DATA C00881 /                                                       F24830
     * 6.65693E+01, 6.68906E+01, 6.72246E+01, 6.75177E+01, 6.78476E+01,   F24840
     * 6.82599E+01, 6.84400E+01, 6.89072E+01, 6.95720E+01, 7.01410E+01,   F24850
     * 7.05519E+01, 7.09367E+01, 7.13975E+01, 7.22128E+01, 7.28222E+01,   F24860
     * 7.33808E+01, 7.38828E+01, 7.44496E+01, 7.49983E+01, 7.54178E+01,   F24870
     * 7.60554E+01, 7.62484E+01, 7.67892E+01, 7.71262E+01, 7.76235E+01,   F24880
     * 7.81413E+01, 7.85694E+01, 7.91248E+01, 7.94715E+01, 7.96200E+01,   F24890
     * 8.00270E+01, 8.03783E+01, 8.07100E+01, 8.11929E+01, 8.17375E+01,   F24900
     * 8.18410E+01, 8.23341E+01, 8.26754E+01, 8.30893E+01, 8.34232E+01,   F24910
     * 8.35533E+01, 8.36017E+01, 8.38589E+01, 8.43366E+01, 8.47593E+01,   F24920
     * 8.51614E+01, 8.55271E+01, 8.58979E+01, 8.64892E+01, 8.74367E+01,   F24930
     * 8.82440E+01, 8.89105E+01, 8.90980E+01, 8.97266E+01, 9.04886E+01,   F24940
     * 9.12709E+01, 9.21243E+01, 9.26673E+01, 9.31331E+01, 9.38190E+01,   F24950
     * 9.44877E+01, 9.50636E+01, 9.57445E+01, 9.65211E+01, 9.68623E+01,   F24960
     * 9.75356E+01, 9.81991E+01, 9.88881E+01, 9.94554E+01, 9.99292E+01,   F24970
     * 1.00357E+02, 1.00670E+02, 1.01227E+02, 1.01529E+02, 1.01889E+02,   F24980
     * 1.02033E+02, 1.02254E+02, 1.02731E+02, 1.02914E+02, 1.03120E+02/   F24990
      DATA C00961 /                                                       F25000
     * 1.03674E+02, 1.03768E+02, 1.04146E+02, 1.04850E+02, 1.05525E+02,   F25010
     * 1.06263E+02, 1.06653E+02, 1.07084E+02, 1.07461E+02, 1.08052E+02,   F25020
     * 1.08793E+02, 1.09395E+02, 1.09811E+02, 1.10079E+02, 1.10656E+02,   F25030
     * 1.11575E+02, 1.12544E+02, 1.13453E+02, 1.14440E+02, 1.15292E+02,   F25040
     * 1.15869E+02, 1.16925E+02, 1.17854E+02, 1.18723E+02, 1.19574E+02,   F25050
     * 1.19940E+02, 1.21108E+02, 1.21807E+02, 1.22490E+02, 1.23278E+02,   F25060
     * 1.24094E+02, 1.24816E+02, 1.25469E+02, 1.26217E+02, 1.26878E+02,   F25070
     * 1.27536E+02, 1.28168E+02, 1.28682E+02, 1.29076E+02, 1.30171E+02,   F25080
     * 1.30667E+02, 1.31242E+02, 1.31665E+02, 1.31961E+02, 1.32347E+02,   F25090
     * 1.32805E+02, 1.33152E+02, 1.33869E+02, 1.34261E+02, 1.34498E+02,   F25100
     * 1.35028E+02, 1.36049E+02, 1.36577E+02, 1.37491E+02, 1.38078E+02,   F25110
     * 1.38389E+02, 1.38819E+02, 1.39653E+02, 1.39770E+02, 1.40812E+02,   F25120
     * 1.40926E+02, 1.41267E+02, 1.41872E+02, 1.42233E+02, 1.43447E+02,   F25130
     * 1.44641E+02, 1.45500E+02, 1.45996E+02, 1.47040E+02, 1.48767E+02,   F25140
     * 1.48785E+02, 1.49525E+02, 1.50266E+02, 1.50814E+02, 1.51443E+02,   F25150
     * 1.52272E+02, 1.52846E+02, 1.54000E+02, 1.54629E+02, 1.54907E+02/   F25160
      DATA C01041 /                                                       F25170
     * 1.55527E+02, 1.56642E+02, 1.57436E+02, 1.59036E+02, 1.59336E+02,   F25180
     * 1.59661E+02, 1.60287E+02, 1.61202E+02, 1.62410E+02, 1.63040E+02,   F25190
     * 1.62872E+02, 1.63248E+02, 1.63776E+02, 1.64313E+02, 1.65782E+02,   F25200
     * 1.65692E+02, 1.66049E+02, 1.66701E+02, 1.67786E+02, 1.69150E+02,   F25210
     * 1.69996E+02, 1.71634E+02, 1.71137E+02, 1.71372E+02, 1.72525E+02,   F25220
     * 1.73816E+02, 1.75219E+02, 1.76091E+02, 1.78260E+02, 1.79299E+02,   F25230
     * 1.79904E+02, 1.81718E+02, 1.83807E+02, 1.85488E+02, 1.85929E+02,   F25240
     * 1.86787E+02, 1.88282E+02, 1.89546E+02, 1.91489E+02, 1.92646E+02,   F25250
     * 1.93399E+02, 1.93838E+02, 1.94406E+02, 1.95829E+02, 1.96745E+02,   F25260
     * 1.96978E+02, 1.97243E+02, 1.97636E+02, 1.98025E+02, 1.98227E+02,   F25270
     * 1.99552E+02, 2.00304E+02, 2.01031E+02, 2.01788E+02, 2.02432E+02,   F25280
     * 2.03817E+02, 2.04866E+02, 2.05561E+02, 2.06180E+02, 2.07024E+02,   F25290
     * 2.08303E+02, 2.09426E+02, 2.10575E+02, 2.11637E+02, 2.12559E+02,   F25300
     * 2.13361E+02, 2.14191E+02, 2.15264E+02, 2.16366E+02, 2.17316E+02,   F25310
     * 2.17717E+02, 2.17154E+02, 2.19172E+02, 2.20346E+02, 2.20849E+02,   F25320
     * 2.21539E+02, 2.22810E+02, 2.22740E+02, 2.22824E+02, 2.23285E+02/   F25330
      DATA C01121 /                                                       F25340
     * 2.23696E+02, 2.23864E+02, 2.23968E+02, 2.23544E+02, 2.24804E+02,   F25350
     * 2.25953E+02, 2.26753E+02, 2.27732E+02, 2.29505E+02, 2.30108E+02,   F25360
     * 2.31232E+02, 2.32552E+02, 2.33979E+02, 2.36677E+02, 2.38481E+02,   F25370
     * 2.41797E+02, 2.44025E+02, 2.45113E+02, 2.47373E+02, 2.47258E+02,   F25380
     * 2.48617E+02, 2.49790E+02, 2.50562E+02, 2.51198E+02, 2.51289E+02,   F25390
     * 2.52509E+02, 2.54136E+02, 2.55335E+02, 2.55808E+02, 2.56567E+02,   F25400
     * 2.57977E+02, 2.58987E+02, 2.59622E+02, 2.60170E+02, 2.61127E+02,   F25410
     * 2.60655E+02, 2.62129E+02, 2.64020E+02, 2.65659E+02, 2.67086E+02,   F25420
     * 2.67615E+02, 2.69800E+02, 2.71452E+02, 2.73314E+02, 2.76972E+02,   F25430
     * 2.78005E+02, 2.79815E+02, 2.81709E+02, 2.84043E+02, 2.87070E+02,   F25440
     * 2.88842E+02, 2.90555E+02, 2.92401E+02, 2.94314E+02, 2.96074E+02,   F25450
     * 2.97103E+02, 2.98037E+02, 2.98113E+02, 2.97705E+02, 2.97350E+02,   F25460
     * 2.97329E+02, 2.97016E+02, 2.96752E+02, 2.96599E+02, 2.96637E+02,   F25470
     * 2.97057E+02, 2.97585E+02, 2.98179E+02, 2.98997E+02, 3.00012E+02,   F25480
     * 3.00806E+02, 3.00908E+02, 3.02369E+02, 3.04063E+02, 3.05325E+02,   F25490
     * 3.06737E+02, 3.08066E+02, 3.09694E+02, 3.11530E+02, 3.13132E+02/   F25500
      DATA C01201 /                                                       F25510
     * 3.13296E+02, 3.15513E+02, 3.16887E+02, 3.17682E+02, 3.18296E+02,   F25520
     * 3.18654E+02, 3.18912E+02, 3.19236E+02, 3.19626E+02, 3.20020E+02,   F25530
     * 3.20186E+02, 3.20709E+02, 3.21628E+02, 3.22625E+02, 3.23504E+02,   F25540
     * 3.25479E+02, 3.26825E+02, 3.28146E+02, 3.29404E+02, 3.30512E+02,   F25550
     * 3.32634E+02, 3.34422E+02, 3.35602E+02, 3.36833E+02, 3.39372E+02,   F25560
     * 3.43446E+02, 3.46374E+02, 3.48719E+02, 3.50881E+02, 3.53160E+02,   F25570
     * 3.54890E+02, 3.57162E+02, 3.59284E+02, 3.60876E+02, 3.62295E+02,   F25580
     * 3.63987E+02, 3.64835E+02, 3.65257E+02, 3.65738E+02, 3.65904E+02,   F25590
     * 3.65976E+02, 3.66460E+02, 3.67087E+02, 3.67377E+02, 3.69079E+02,   F25600
     * 3.70694E+02, 3.70940E+02, 3.70557E+02, 3.72693E+02, 3.73852E+02,   F25610
     * 3.75679E+02, 3.77863E+02, 3.79964E+02, 3.81368E+02, 3.82716E+02,   F25620
     * 3.85556E+02, 3.89072E+02, 3.91796E+02, 3.92766E+02, 3.96551E+02,   F25630
     * 3.97833E+02, 3.97285E+02, 4.01929E+02, 4.02158E+02, 4.04553E+02,   F25640
     * 4.06451E+02, 4.06236E+02, 4.08135E+02, 4.07797E+02, 4.08415E+02,   F25650
     * 4.10111E+02, 4.11781E+02, 4.12735E+02, 4.11547E+02, 4.11606E+02,   F25660
     * 4.13548E+02, 4.12557E+02, 4.12923E+02, 4.12866E+02, 4.13009E+02/   F25670
      DATA C01281 /                                                       F25680
     * 4.14447E+02, 4.16032E+02, 4.17032E+02, 4.19064E+02, 4.22458E+02,   F25690
     * 4.26021E+02, 4.25192E+02, 4.25684E+02, 4.27536E+02, 4.29972E+02,   F25700
     * 4.31994E+02, 4.36037E+02, 4.39132E+02, 4.40363E+02, 4.40716E+02,   F25710
     * 4.40342E+02, 4.42063E+02, 4.44408E+02, 4.45454E+02, 4.47835E+02,   F25720
     * 4.48256E+02, 4.48831E+02, 4.50257E+02, 4.51427E+02, 4.52373E+02,   F25730
     * 4.53899E+02, 4.55496E+02, 4.56311E+02, 4.57314E+02, 4.59922E+02,   F25740
     * 4.61048E+02, 4.59840E+02, 4.62144E+02, 4.63152E+02, 4.64565E+02,   F25750
     * 4.66715E+02, 4.69380E+02, 4.70751E+02, 4.72012E+02, 4.73482E+02,   F25760
     * 4.75524E+02, 4.79307E+02, 4.82035E+02, 4.84423E+02, 4.86712E+02,   F25770
     * 4.88754E+02, 4.90102E+02, 4.92047E+02, 4.94150E+02, 4.95375E+02,   F25780
     * 4.95828E+02, 4.97555E+02, 4.98559E+02, 4.97618E+02, 4.99265E+02,   F25790
     * 4.99979E+02, 5.00681E+02, 5.01386E+02, 5.00868E+02, 5.01935E+02,   F25800
     * 5.03151E+02, 5.04329E+02, 5.05546E+02, 5.08259E+02, 5.09222E+02,   F25810
     * 5.09818E+02, 5.11397E+02, 5.12391E+02, 5.13326E+02, 5.14329E+02,   F25820
     * 5.15443E+02, 5.16533E+02, 5.21417E+02, 5.25071E+02, 5.26581E+02,   F25830
     * 5.27762E+02, 5.29274E+02, 5.31704E+02, 5.34310E+02, 5.35727E+02/   F25840
      DATA C01361 /                                                       F25850
     * 5.36838E+02, 5.37082E+02, 5.36733E+02, 5.36170E+02, 5.36063E+02,   F25860
     * 5.36451E+02, 5.37870E+02, 5.40475E+02, 5.42268E+02, 5.41972E+02,   F25870
     * 5.42532E+02, 5.44764E+02, 5.46844E+02, 5.47525E+02, 5.49150E+02,   F25880
     * 5.52049E+02, 5.55423E+02, 5.56259E+02, 5.57424E+02, 5.59189E+02,   F25890
     * 5.61167E+02, 5.64512E+02, 5.66753E+02, 5.68183E+02, 5.69628E+02,   F25900
     * 5.73474E+02, 5.76192E+02, 5.78058E+02, 5.79588E+02, 5.81619E+02,   F25910
     * 5.83530E+02, 5.84852E+02, 5.85326E+02, 5.88130E+02, 5.90570E+02,   F25920
     * 5.91785E+02, 5.91371E+02, 5.90931E+02, 5.90942E+02, 5.91168E+02,   F25930
     * 5.91291E+02, 5.89791E+02, 5.91146E+02, 5.90804E+02, 5.87847E+02,   F25940
     * 5.89067E+02, 5.91027E+02, 5.90951E+02, 5.89227E+02, 5.93389E+02,   F25950
     * 5.92921E+02, 5.92739E+02, 5.94544E+02, 5.98941E+02, 6.02302E+02,   F25960
     * 6.03908E+02, 6.04265E+02, 6.06737E+02, 6.08560E+02, 6.11272E+02,   F25970
     * 6.14992E+02, 6.18595E+02, 6.20930E+02, 6.22107E+02, 6.22957E+02,   F25980
     * 6.26710E+02, 6.28657E+02, 6.30132E+02, 6.31543E+02, 6.33043E+02,   F25990
     * 6.36932E+02, 6.38248E+02, 6.37126E+02, 6.41648E+02, 6.48274E+02,   F26000
     * 6.52638E+02, 6.53922E+02, 6.56647E+02, 6.59351E+02, 6.60525E+02/   F26010
      DATA C01441 /                                                       F26020
     * 6.60130E+02, 6.61375E+02, 6.62660E+02, 6.63976E+02, 6.65181E+02,   F26030
     * 6.64820E+02, 6.64458E+02, 6.64927E+02, 6.66555E+02, 6.66759E+02,   F26040
     * 6.68218E+02, 6.70323E+02, 6.72703E+02, 6.76085E+02, 6.79180E+02,   F26050
     * 6.80850E+02, 6.80017E+02, 6.79928E+02, 6.80886E+02, 6.82038E+02,   F26060
     * 6.82271E+02, 6.84057E+02, 6.85309E+02, 6.86816E+02, 6.90180E+02,   F26070
     * 6.93205E+02, 6.95870E+02, 6.98794E+02, 7.03776E+02, 7.04010E+02,   F26080
     * 7.05041E+02, 7.07254E+02, 7.07432E+02, 7.10736E+02, 7.13791E+02,   F26090
     * 7.15542E+02, 7.16468E+02, 7.17412E+02, 7.17783E+02, 7.17340E+02,   F26100
     * 7.18184E+02, 7.18716E+02, 7.18809E+02, 7.18282E+02, 7.20317E+02,   F26110
     * 7.18568E+02, 7.16274E+02, 7.19119E+02, 7.20852E+02, 7.21727E+02,   F26120
     * 7.22607E+02, 7.26369E+02, 7.26412E+02, 7.27101E+02, 7.29404E+02,   F26130
     * 7.30786E+02, 7.30910E+02, 7.30656E+02, 7.30566E+02, 7.33408E+02,   F26140
     * 7.37064E+02, 7.39178E+02, 7.36713E+02, 7.37365E+02, 7.40861E+02,   F26150
     * 7.45281E+02, 7.46178E+02, 7.46991E+02, 7.48035E+02, 7.49777E+02,   F26160
     * 7.54665E+02, 7.56585E+02, 7.57408E+02, 7.58131E+02, 7.58155E+02,   F26170
     * 7.60838E+02, 7.64792E+02, 7.68161E+02, 7.69263E+02, 7.73166E+02/   F26180
      DATA C01521 /                                                       F26190
     * 7.79006E+02, 7.82037E+02, 7.83109E+02, 7.84674E+02, 7.87444E+02,   F26200
     * 7.89510E+02, 7.90130E+02, 7.91364E+02, 7.95225E+02, 8.03599E+02,   F26210
     * 8.06340E+02, 8.05105E+02, 8.05120E+02, 8.08515E+02, 8.10907E+02,   F26220
     * 8.11388E+02, 8.13432E+02, 8.12579E+02, 8.10564E+02, 8.08719E+02,   F26230
     * 8.07682E+02, 8.05009E+02, 8.01754E+02, 8.01013E+02, 7.99926E+02,   F26240
     * 7.99067E+02, 7.98369E+02, 7.94090E+02, 7.92883E+02, 7.94244E+02,   F26250
     * 7.98220E+02, 7.98201E+02, 7.98332E+02, 7.99289E+02, 8.02355E+02,   F26260
     * 8.03621E+02, 8.05302E+02, 8.08368E+02, 8.09983E+02, 8.11529E+02,   F26270
     * 8.13068E+02, 8.14717E+02, 8.16441E+02, 8.19241E+02, 8.22944E+02,   F26280
     * 8.23768E+02, 8.25030E+02, 8.26103E+02, 8.26374E+02, 8.28331E+02,   F26290
     * 8.32620E+02, 8.38618E+02, 8.43666E+02, 8.45212E+02, 8.46324E+02,   F26300
     * 8.48536E+02, 8.50192E+02, 8.53083E+02, 8.56653E+02, 8.59614E+02,   F26310
     * 8.62000E+02, 8.64593E+02, 8.67678E+02, 8.70908E+02, 8.73408E+02,   F26320
     * 8.74779E+02, 8.74005E+02, 8.76718E+02, 8.80445E+02, 8.84365E+02,   F26330
     * 8.83806E+02, 8.84292E+02, 8.85539E+02, 8.87474E+02, 8.84905E+02,   F26340
     * 8.84039E+02, 8.85105E+02, 8.83733E+02, 8.82224E+02, 8.79865E+02/   F26350
      DATA C01601 /                                                       F26360
     * 8.75663E+02, 8.75575E+02, 8.73144E+02, 8.68602E+02, 8.70278E+02,   F26370
     * 8.69659E+02, 8.68701E+02, 8.69250E+02, 8.71057E+02, 8.72860E+02,   F26380
     * 8.74361E+02, 8.74458E+02, 8.77576E+02, 8.81613E+02, 8.84358E+02,   F26390
     * 8.87440E+02, 8.91549E+02, 8.96568E+02, 8.99836E+02, 9.02880E+02,   F26400
     * 9.05428E+02, 9.06891E+02, 9.07349E+02, 9.10151E+02, 9.15917E+02,   F26410
     * 9.16197E+02, 9.18571E+02, 9.21219E+02, 9.20292E+02, 9.21949E+02,   F26420
     * 9.24509E+02, 9.27454E+02, 9.29474E+02, 9.31348E+02, 9.32818E+02,   F26430
     * 9.32658E+02, 9.36280E+02, 9.39512E+02, 9.39667E+02, 9.44078E+02,   F26440
     * 9.47196E+02, 9.48291E+02, 9.46150E+02, 9.46918E+02, 9.49093E+02,   F26450
     * 9.51372E+02, 9.53109E+02, 9.56308E+02, 9.61335E+02, 9.58214E+02,   F26460
     * 9.56188E+02, 9.55660E+02, 9.58633E+02, 9.57541E+02, 9.54879E+02,   F26470
     * 9.51663E+02, 9.52839E+02, 9.52055E+02, 9.49253E+02, 9.50187E+02,   F26480
     * 9.50323E+02, 9.50937E+02, 9.54362E+02, 9.55855E+02, 9.56350E+02,   F26490
     * 9.55908E+02, 9.57963E+02, 9.61866E+02, 9.66948E+02, 9.69786E+02,   F26500
     * 9.74302E+02, 9.79061E+02, 9.82465E+02, 9.86019E+02, 9.89930E+02,   F26510
     * 9.94294E+02, 9.97011E+02, 9.98207E+02, 9.98607E+02, 1.00175E+03/   F26520
      DATA C01681 /                                                       F26530
     * 1.00275E+03, 1.00284E+03, 1.00294E+03, 1.00485E+03, 1.00593E+03,   F26540
     * 1.00524E+03, 1.00415E+03, 1.00335E+03, 1.00278E+03, 1.00185E+03,   F26550
     * 9.99982E+02, 9.98177E+02, 9.97959E+02, 9.99161E+02, 9.98810E+02,   F26560
     * 9.95415E+02, 9.94342E+02, 9.92998E+02, 9.91340E+02, 9.90900E+02,   F26570
     * 9.90407E+02, 9.89232E+02, 9.85447E+02, 9.86312E+02, 9.87461E+02,   F26580
     * 9.86090E+02, 9.86670E+02, 9.85534E+02, 9.81877E+02, 9.84946E+02,   F26590
     * 9.86392E+02, 9.86709E+02, 9.88086E+02, 9.90269E+02, 9.92566E+02,   F26600
     * 9.94029E+02, 9.95795E+02, 9.97788E+02, 1.00005E+03, 1.00287E+03,   F26610
     * 1.00566E+03, 1.00833E+03, 1.00982E+03, 1.01348E+03, 1.01862E+03,   F26620
     * 1.02322E+03, 1.02786E+03, 1.03179E+03, 1.03339E+03, 1.03833E+03,   F26630
     * 1.04317E+03, 1.04598E+03, 1.04753E+03, 1.04981E+03, 1.05321E+03,   F26640
     * 1.05492E+03, 1.05721E+03, 1.05978E+03, 1.06033E+03, 1.06107E+03,   F26650
     * 1.06155E+03, 1.06035E+03, 1.05838E+03, 1.05649E+03, 1.05553E+03,   F26660
     * 1.05498E+03, 1.05387E+03, 1.05171E+03, 1.04877E+03, 1.04725E+03,   F26670
     * 1.04748E+03, 1.04733E+03, 1.04704E+03, 1.04643E+03, 1.04411E+03,   F26680
     * 1.04435E+03, 1.04520E+03, 1.04233E+03, 1.04047E+03, 1.03992E+03/   F26690
      DATA C01761 /                                                       F26700
     * 1.04192E+03, 1.04171E+03, 1.04140E+03, 1.04197E+03, 1.04415E+03,   F26710
     * 1.04548E+03, 1.04533E+03, 1.04616E+03, 1.04705E+03, 1.04800E+03,   F26720
     * 1.05025E+03, 1.05219E+03, 1.05412E+03, 1.05808E+03, 1.06062E+03,   F26730
     * 1.06292E+03, 1.06780E+03, 1.07219E+03, 1.07610E+03, 1.07913E+03,   F26740
     * 1.08405E+03, 1.08798E+03, 1.08835E+03, 1.09140E+03, 1.09447E+03,   F26750
     * 1.09676E+03, 1.10015E+03, 1.10272E+03, 1.10410E+03, 1.10749E+03,   F26760
     * 1.10991E+03, 1.11121E+03, 1.10981E+03, 1.10981E+03, 1.11063E+03,   F26770
     * 1.10714E+03, 1.10500E+03, 1.10357E+03, 1.10093E+03, 1.09898E+03,   F26780
     * 1.09679E+03, 1.09188E+03, 1.09088E+03, 1.09040E+03, 1.08586E+03,   F26790
     * 1.08178E+03, 1.07752E+03, 1.07243E+03, 1.07178E+03, 1.07084E+03,   F26800
     * 1.06693E+03, 1.06527E+03, 1.06405E+03, 1.06285E+03, 1.06287E+03,   F26810
     * 1.06276E+03, 1.06221E+03, 1.06464E+03, 1.06579E+03, 1.06498E+03,   F26820
     * 1.06596E+03, 1.06812E+03, 1.07159E+03, 1.07361E+03, 1.07556E+03,   F26830
     * 1.07751E+03, 1.08128E+03, 1.08523E+03, 1.08927E+03, 1.09193E+03,   F26840
     * 1.09612E+03, 1.10133E+03, 1.10435E+03, 1.10781E+03, 1.11168E+03,   F26850
     * 1.11641E+03, 1.12217E+03, 1.12839E+03, 1.13298E+03, 1.13575E+03/   F26860
      DATA C01841 /                                                       F26870
     * 1.13742E+03, 1.13929E+03, 1.14132E+03, 1.14340E+03, 1.14518E+03,   F26880
     * 1.14742E+03, 1.14943E+03, 1.14935E+03, 1.14975E+03, 1.15086E+03,   F26890
     * 1.15420E+03, 1.15267E+03, 1.15007E+03, 1.15155E+03, 1.14982E+03,   F26900
     * 1.14663E+03, 1.14301E+03, 1.13986E+03, 1.13676E+03, 1.13307E+03,   F26910
     * 1.12898E+03, 1.12516E+03, 1.12284E+03, 1.12068E+03, 1.11855E+03,   F26920
     * 1.11632E+03, 1.11464E+03, 1.11318E+03, 1.11180E+03, 1.11163E+03,   F26930
     * 1.11160E+03, 1.11035E+03, 1.11178E+03, 1.11395E+03, 1.11447E+03,   F26940
     * 1.11439E+03, 1.11440E+03, 1.11582E+03, 1.11560E+03, 1.11478E+03,   F26950
     * 1.11448E+03, 1.11454E+03, 1.11494E+03, 1.11607E+03, 1.11736E+03,   F26960
     * 1.11854E+03, 1.11875E+03, 1.11989E+03, 1.12165E+03, 1.12427E+03,   F26970
     * 1.12620E+03, 1.12758E+03, 1.12774E+03, 1.12870E+03, 1.13001E+03,   F26980
     * 1.13006E+03, 1.13078E+03, 1.13172E+03, 1.12971E+03, 1.12857E+03,   F26990
     * 1.12810E+03, 1.12740E+03, 1.12659E+03, 1.12564E+03, 1.12338E+03,   F27000
     * 1.12117E+03, 1.11902E+03, 1.11878E+03, 1.11855E+03, 1.11828E+03,   F27010
     * 1.11791E+03, 1.11784E+03, 1.11815E+03, 1.11957E+03, 1.12046E+03,   F27020
     * 1.12042E+03, 1.11929E+03, 1.12074E+03, 1.12708E+03, 1.12600E+03/   F27030
      DATA C01921 /                                                       F27040
     * 1.12538E+03, 1.12871E+03, 1.13167E+03, 1.13388E+03, 1.13444E+03,   F27050
     * 1.13595E+03, 1.13801E+03, 1.14096E+03, 1.14230E+03, 1.14304E+03,   F27060
     * 1.14421E+03, 1.14580E+03, 1.14767E+03, 1.15000E+03, 1.15126E+03,   F27070
     * 1.15181E+03, 1.15197E+03, 1.15364E+03, 1.15626E+03, 1.15538E+03,   F27080
     * 1.15636E+03, 1.15908E+03, 1.16024E+03, 1.16188E+03, 1.16411E+03,   F27090
     * 1.16310E+03, 1.16430E+03, 1.16927E+03, 1.17035E+03, 1.17052E+03,   F27100
     * 1.17013E+03, 1.16968E+03, 1.16969E+03, 1.17106E+03, 1.17123E+03,   F27110
     * 1.17006E+03, 1.16536E+03, 1.16087E+03, 1.15691E+03, 1.15608E+03,   F27120
     * 1.15388E+03, 1.15077E+03, 1.14967E+03, 1.14793E+03, 1.14554E+03,   F27130
     * 1.14212E+03, 1.13908E+03, 1.13654E+03, 1.13499E+03, 1.13308E+03,   F27140
     * 1.13033E+03, 1.13051E+03, 1.13073E+03, 1.12898E+03, 1.12941E+03,   F27150
     * 1.13051E+03, 1.13086E+03, 1.13189E+03, 1.13304E+03, 1.13192E+03,   F27160
     * 1.13131E+03, 1.13110E+03, 1.13499E+03, 1.13914E+03, 1.14359E+03,   F27170
     * 1.14383E+03, 1.14390E+03, 1.14435E+03, 1.14540E+03, 1.14646E+03,   F27180
     * 1.14716E+03, 1.14880E+03, 1.15062E+03, 1.15170E+03, 1.15093E+03,   F27190
     * 1.14926E+03, 1.15133E+03, 1.15167E+03, 1.15043E+03, 1.15134E+03/   F27200
      DATA C02001 /                                                       F27210
     * 1.15135E+03, 1.15000E+03, 1.15087E+03, 1.15118E+03, 1.14935E+03,   F27220
     * 1.14780E+03, 1.14647E+03, 1.14560E+03, 1.14404E+03, 1.14238E+03,   F27230
     * 1.14406E+03, 1.14245E+03, 1.13781E+03, 1.13664E+03, 1.13653E+03,   F27240
     * 1.13778E+03, 1.13813E+03, 1.13794E+03, 1.13681E+03, 1.13515E+03,   F27250
     * 1.13328E+03, 1.13132E+03, 1.13080E+03, 1.13130E+03, 1.13400E+03,   F27260
     * 1.13526E+03, 1.13494E+03, 1.13193E+03, 1.12898E+03, 1.12654E+03,   F27270
     * 1.12739E+03, 1.12849E+03, 1.12774E+03, 1.12733E+03, 1.12733E+03,   F27280
     * 1.12943E+03, 1.13014E+03, 1.12967E+03, 1.12731E+03, 1.12671E+03,   F27290
     * 1.12885E+03, 1.13050E+03, 1.13201E+03, 1.13345E+03, 1.13488E+03,   F27300
     * 1.13605E+03, 1.13530E+03, 1.13737E+03, 1.14186E+03, 1.14250E+03,   F27310
     * 1.14305E+03, 1.14383E+03, 1.14510E+03, 1.14659E+03, 1.14848E+03,   F27320
     * 1.14949E+03, 1.14995E+03, 1.14934E+03, 1.15058E+03, 1.15368E+03,   F27330
     * 1.15435E+03, 1.15422E+03, 1.15296E+03, 1.15228E+03, 1.15189E+03,   F27340
     * 1.15198E+03, 1.15081E+03, 1.14881E+03, 1.14562E+03, 1.14276E+03,   F27350
     * 1.14030E+03, 1.13637E+03, 1.13254E+03, 1.12942E+03, 1.12653E+03,   F27360
     * 1.12362E+03, 1.11987E+03, 1.11712E+03, 1.11522E+03, 1.11403E+03/   F27370
      DATA C02081 /                                                       F27380
     * 1.11226E+03, 1.10947E+03, 1.10956E+03, 1.10976E+03, 1.10748E+03,   F27390
     * 1.10673E+03, 1.10688E+03, 1.10675E+03, 1.10533E+03, 1.10230E+03,   F27400
     * 1.10384E+03, 1.10496E+03, 1.10274E+03, 1.10197E+03, 1.10196E+03,   F27410
     * 1.10278E+03, 1.10257E+03, 1.10147E+03, 1.10205E+03, 1.10308E+03,   F27420
     * 1.10478E+03, 1.10358E+03, 1.10197E+03, 1.10305E+03, 1.10390E+03,   F27430
     * 1.10456E+03, 1.10526E+03, 1.10588E+03, 1.10640E+03, 1.10747E+03,   F27440
     * 1.10904E+03, 1.11214E+03, 1.11350E+03, 1.11359E+03, 1.11604E+03,   F27450
     * 1.11706E+03, 1.11594E+03, 1.11600E+03, 1.11616E+03, 1.11561E+03,   F27460
     * 1.11556E+03, 1.11547E+03, 1.11370E+03, 1.11289E+03, 1.11276E+03,   F27470
     * 1.11338E+03, 1.11437E+03, 1.11595E+03, 1.11309E+03, 1.10958E+03,   F27480
     * 1.10887E+03, 1.10573E+03, 1.10068E+03, 1.10194E+03, 1.10165E+03,   F27490
     * 1.09813E+03, 1.09973E+03, 1.10233E+03, 1.10121E+03, 1.10097E+03,   F27500
     * 1.10149E+03, 1.10162E+03, 1.10222E+03, 1.10389E+03, 1.10315E+03,   F27510
     * 1.10158E+03, 1.10193E+03, 1.10186E+03, 1.10135E+03, 1.10336E+03,   F27520
     * 1.10500E+03, 1.10459E+03, 1.10592E+03, 1.10784E+03, 1.10076E+03,   F27530
     * 1.09615E+03, 1.09496E+03, 1.09422E+03, 1.09350E+03, 1.09244E+03/   F27540
      DATA C02161 /                                                       F27550
     * 1.08955E+03, 1.08535E+03, 1.08379E+03, 1.08184E+03, 1.07889E+03,   F27560
     * 1.07563E+03, 1.07238E+03, 1.07042E+03, 1.06882E+03, 1.06761E+03,   F27570
     * 1.06816E+03, 1.06772E+03, 1.06327E+03, 1.06313E+03, 1.06563E+03,   F27580
     * 1.06254E+03, 1.06072E+03, 1.06095E+03, 1.06173E+03, 1.06269E+03,   F27590
     * 1.06361E+03, 1.06438E+03, 1.06501E+03, 1.06465E+03, 1.06481E+03,   F27600
     * 1.06685E+03, 1.06642E+03, 1.06447E+03, 1.06701E+03, 1.06791E+03,   F27610
     * 1.06612E+03, 1.06471E+03, 1.06403E+03, 1.06774E+03, 1.06823E+03,   F27620
     * 1.06524E+03, 1.06479E+03, 1.06453E+03, 1.06346E+03, 1.06175E+03,   F27630
     * 1.05958E+03, 1.05941E+03, 1.05936E+03, 1.05938E+03, 1.05736E+03,   F27640
     * 1.05449E+03, 1.05307E+03, 1.05180E+03, 1.05074E+03, 1.04810E+03,   F27650
     * 1.04536E+03, 1.04477E+03, 1.04389E+03, 1.04272E+03, 1.04006E+03,   F27660
     * 1.03739E+03, 1.03533E+03, 1.03476E+03, 1.03516E+03, 1.03275E+03,   F27670
     * 1.03093E+03, 1.03062E+03, 1.02997E+03, 1.02919E+03, 1.02993E+03,   F27680
     * 1.02983E+03, 1.02837E+03, 1.02611E+03, 1.02386E+03, 1.02426E+03,   F27690
     * 1.02542E+03, 1.02750E+03, 1.02638E+03, 1.02496E+03, 1.02608E+03,   F27700
     * 1.02568E+03, 1.02388E+03, 1.02522E+03, 1.02692E+03, 1.02834E+03/   F27710
      DATA C02241 /                                                       F27720
     * 1.02828E+03, 1.02716E+03, 1.02667E+03, 1.02607E+03, 1.02503E+03,   F27730
     * 1.02723E+03, 1.03143E+03, 1.02881E+03, 1.02646E+03, 1.02500E+03,   F27740
     * 1.02569E+03, 1.02743E+03, 1.02608E+03, 1.02548E+03, 1.02620E+03,   F27750
     * 1.02733E+03, 1.02839E+03, 1.02575E+03, 1.02432E+03, 1.02471E+03,   F27760
     * 1.02392E+03, 1.02267E+03, 1.02077E+03, 1.01964E+03, 1.01957E+03,   F27770
     * 1.01848E+03, 1.01704E+03, 1.01524E+03, 1.01352E+03, 1.01191E+03,   F27780
     * 1.01066E+03, 1.00952E+03, 1.00849E+03, 1.00660E+03, 1.00368E+03,   F27790
     * 9.99713E+02, 9.95921E+02, 9.94845E+02, 9.93286E+02, 9.91204E+02/   F27800
C                                                                         F27810
      END                                                                 F27820
C
C     --------------------------------------------------------------
C
      SUBROUTINE O3HHT1 (V1C,V2C,DVC,NPTC,C)                              F27830
C                                                                         F27840
      IMPLICIT REAL*8           (V)                                     ! F27850
C                                                                         F27860
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F27870
      COMMON /O3HH1/ V1S,V2S,DVS,NPTS,S(2687)                             F27880
      DIMENSION C(*)                                                      F27890
C                                                                         F27900
      DVC = DVS                                                           F27910
      V1C = V1ABS-DVC                                                     F27920
      V2C = V2ABS+DVC                                                     F27930
C                                                                         F27940
      I1 = (V1C-V1S)/DVS                                                  F27950
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F28020
         I = I1+J                                                         F28030
         C(J) = 0.                                                        F28040
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F28050
         C(J) = S(I)                                                      F28060
   10 CONTINUE                                                            F28070
C                                                                         F28080
      RETURN                                                              F28090
C                                                                         F28100
      END                                                                 F28110
C
C     --------------------------------------------------------------
C
      BLOCK DATA BO3HH1                                                   F28120
C                                                                         F28130
      IMPLICIT REAL*8           (V)                                     ! F28140
C                                                                         F28150
C     RATIO (C1/C0)                                                       F28160
C     DATA FROM BASS 1985                                                 F28170
C                                                                         F28180
C     NOW INCLUDES MOLINA & MOLINA AT 273K WITH THE TEMPERATURE           F28190
C     DEPENDENCE DETERMINED FROM THE 195K HARVARD MEASUREMENTS,           F28200
C     EMPLOYING THE BASS ALGORITHM
C
C              (CO(1+C1*(T-273.15)+C2*(T-273.15)**2);                     F28210
C
C     THIS IS ONLY FOR THE WAVELENGTH RANGE FROM .34 TO .35 MICRONS;      F28220
C     OTHERWISE, THE BASS DATA ALONE HAVE BEEN EMPLOYED BETWEEN           F28230
C     .34 AND .245 MICRONS.                                               F28240
C                                                                         F28250
C     NEW T-DEPENDENT X-SECTIONS BETWEEN .345 AND .36 MICRONS             F28260
C     HAVE NOW BEEN ADDED, BASED ON WORK BY CACCIANI, DISARRA             F28270
C     AND FIOCCO, UNIVERSITY OF ROME, 1987.  QUADRATIC TEMP               F28280
C     HAS BEEN DERIVED, AS ABOVE.                                         F28290
C                                                                         F28300
C     AGREEMENT AMONGST THE FOUR DATA SETS IS REASONABLE (<10%)           F28310
C     AND OFTEN EXCELLENT (0-3%)                                          F28320
C                                                                         F28330
C                                                                         F28340
      COMMON /O3HH1/ V1C,V2C,DVC,NC,                                      F28350
     *               O31001(85),C10086(80),C10166(80),C10246(65),         F28360
     *               C10311(16),C10327(80),C10407( 1),                    F28370
     *               C10001(80),C10081(80),C10161(80),C10241(80),         F28380
     *               C10321(80),C10401(80),C10481(80),C10561(80),         F28390
     *               C10641(80),C10721(80),C10801(80),C10881(80),         F28400
     *               C10961(80),C11041(80),C11121(80),C11201(80),         F28410
     *               C11281(80),C11361(80),C11441(80),C11521(80),         F28420
     *               C11601(80),C11681(80),C11761(80),C11841(80),         F28430
     *               C11921(80),C12001(80),C12081(80),C12161(80),         F28440
     *               C12241(40)                                           F28450
C                                                                         F28460
C     DATA V1C /29405./, V2C /40800./ ,DVC /5./, NC /2280/   BASS         F28470
C                                                                         F28480
      DATA V1C /27370./, V2C /40800./ ,DVC /5./, NC /2687/                F28490
C                                                                         F28500
      DATA O31001/85*1.3E-3/                                              F28510
C                                                                         F28520
      DATA C10086/                                                        F28530
     * 1.37330E-03, 1.62821E-03, 2.01703E-03, 2.54574E-03, 3.20275E-03,   F28540
     * 3.89777E-03, 4.62165E-03, 5.26292E-03, 5.86986E-03, 6.41494E-03,   F28550
     * 6.96761E-03, 7.48539E-03, 7.89600E-03, 7.87305E-03, 7.81981E-03,   F28560
     * 7.63864E-03, 7.67455E-03, 7.72586E-03, 7.69784E-03, 7.57367E-03,   F28570
     * 7.27336E-03, 7.14064E-03, 7.24207E-03, 7.09851E-03, 6.93654E-03,   F28580
     * 6.89385E-03, 7.05768E-03, 6.85578E-03, 6.58301E-03, 6.50848E-03,   F28590
     * 6.52083E-03, 6.46590E-03, 6.70692E-03, 6.92053E-03, 7.17734E-03,   F28600
     * 7.05364E-03, 6.63440E-03, 6.54702E-03, 6.27173E-03, 5.98150E-03,   F28610
     * 5.66579E-03, 5.51549E-03, 5.50291E-03, 5.93271E-03, 6.36950E-03,   F28620
     * 7.18562E-03, 7.51767E-03, 6.53815E-03, 7.22341E-03, 8.63056E-03,   F28630
     * 9.11740E-03, 8.80903E-03, 8.59902E-03, 7.74287E-03, 7.33509E-03,   F28640
     * 7.50180E-03, 7.81686E-03, 7.85635E-03, 8.08554E-03, 7.21968E-03,   F28650
     * 7.99028E-03, 9.90724E-03, 1.29121E-02, 1.54686E-02, 1.60876E-02,   F28660
     * 1.59530E-02, 1.57040E-02, 1.59499E-02, 1.63961E-02, 1.72670E-02,   F28670
     * 1.81634E-02, 1.95519E-02, 2.14181E-02, 2.28670E-02, 2.33506E-02,   F28680
     * 2.22736E-02, 2.14296E-02, 2.15271E-02, 2.30730E-02, 2.36220E-02/   F28690
      DATA C10166/                                                        F28700
     * 2.44466E-02, 2.44476E-02, 2.39223E-02, 2.41386E-02, 2.53687E-02,   F28710
     * 2.67491E-02, 2.80425E-02, 2.77558E-02, 2.82626E-02, 2.86776E-02,   F28720
     * 2.88781E-02, 2.89248E-02, 2.89983E-02, 2.85534E-02, 2.87102E-02,   F28730
     * 2.83695E-02, 2.76719E-02, 2.76091E-02, 2.90733E-02, 2.80388E-02,   F28740
     * 2.73706E-02, 2.65055E-02, 2.61268E-02, 2.45892E-02, 2.37213E-02,   F28750
     * 2.22542E-02, 2.10116E-02, 2.02852E-02, 1.97635E-02, 1.94079E-02,   F28760
     * 1.90997E-02, 1.85598E-02, 1.79221E-02, 1.77887E-02, 1.73709E-02,   F28770
     * 1.67263E-02, 1.60932E-02, 1.50775E-02, 1.39563E-02, 1.23691E-02,   F28780
     * 1.07402E-02, 9.35859E-03, 8.43786E-03, 7.92075E-03, 7.33239E-03,   F28790
     * 6.73638E-03, 6.28740E-03, 5.85640E-03, 5.85384E-03, 6.10577E-03,   F28800
     * 7.26050E-03, 9.66384E-03, 1.29629E-02, 1.69596E-02, 2.03465E-02,   F28810
     * 2.26429E-02, 2.39653E-02, 2.47970E-02, 2.51993E-02, 2.51383E-02,   F28820
     * 2.52014E-02, 2.47766E-02, 2.47171E-02, 2.47478E-02, 2.43986E-02,   F28830
     * 2.43498E-02, 2.40537E-02, 2.40574E-02, 2.40446E-02, 2.40847E-02,   F28840
     * 2.39400E-02, 2.42127E-02, 2.47123E-02, 2.52914E-02, 2.52103E-02,   F28850
     * 2.51421E-02, 2.43229E-02, 2.37902E-02, 2.30865E-02, 2.28174E-02/   F28860
      DATA C10246/                                                        F28870
     * 2.28830E-02, 2.33671E-02, 2.38274E-02, 2.46699E-02, 2.56739E-02,   F28880
     * 2.61408E-02, 2.62898E-02, 2.64228E-02, 2.55561E-02, 2.47095E-02,   F28890
     * 2.39071E-02, 2.34319E-02, 2.28738E-02, 2.23434E-02, 2.18888E-02,   F28900
     * 2.13639E-02, 2.11937E-02, 2.10110E-02, 2.07672E-02, 2.00697E-02,   F28910
     * 1.97605E-02, 1.91208E-02, 1.82056E-02, 1.73945E-02, 1.64542E-02,   F28920
     * 1.53969E-02, 1.41816E-02, 1.35665E-02, 1.27109E-02, 1.18254E-02,   F28930
     * 1.11489E-02, 1.03984E-02, 1.00760E-02, 9.86649E-03, 9.76766E-03,   F28940
     * 9.41662E-03, 9.19082E-03, 9.44272E-03, 1.04547E-02, 1.24713E-02,   F28950
     * 1.49310E-02, 1.70272E-02, 1.86057E-02, 1.93555E-02, 1.98350E-02,   F28960
     * 2.00041E-02, 2.01233E-02, 2.01917E-02, 1.98918E-02, 1.96649E-02,   F28970
     * 1.95162E-02, 2.01044E-02, 2.06711E-02, 2.08881E-02, 2.04812E-02,   F28980
     * 1.92249E-02, 1.80188E-02, 1.69496E-02, 1.60488E-02, 1.52865E-02,   F28990
     * 1.46940E-02, 1.41067E-02, 1.35675E-02, 1.31094E-02, 1.27542E-02/   F29000
      DATA C10311/                                                        F29010
     *                                                     1.3073E-02,    F29020
     * 1.2795E-02,  1.2753E-02,  1.2868E-02,  1.2885E-02,  1.2554E-02,    F29030
     * 1.2106E-02,  1.1616E-02,  1.1394E-02,  1.1092E-02,  1.0682E-02,    F29040
     * 1.0519E-02,  9.7219E-03,  9.3434E-03,  8.5260E-03,  8.3333E-03/    F29050
      DATA C10327/                                                        F29060
     * 7.8582E-03,  6.8295E-03,  6.7963E-03,  6.7516E-03,  6.2930E-03,    F29070
     * 6.1615E-03,  6.1250E-03,  5.9011E-03,  5.7823E-03,  5.4688E-03,    F29080
     * 5.0978E-03,  4.4526E-03,  3.8090E-03,  3.2310E-03,  3.0128E-03,    F29090
     * 3.9063E-03,  6.7911E-03,  9.3161E-03,  1.0256E-02,  1.0183E-02,    F29100
     * 9.8289E-03,  9.5683E-03,  9.0406E-03,  8.7148E-03,  8.5284E-03,    F29110
     * 8.6149E-03,  8.7238E-03,  9.3679E-03,  1.0683E-02,  1.2016E-02,    F29120
     * 1.3097E-02,  1.3610E-02,  1.3588E-02,  1.3805E-02,  1.3928E-02,    F29130
     * 1.3903E-02,  1.3446E-02,  1.3258E-02,  1.3194E-02,  1.2703E-02,    F29140
     * 1.2393E-02,  1.2487E-02,  1.2341E-02,  1.2388E-02,  1.2061E-02,    F29150
     * 1.2122E-02,  1.1850E-02,  1.2032E-02,  1.1806E-02,  1.1810E-02,    F29160
     * 1.1572E-02,  1.1397E-02,  1.0980E-02,  1.1012E-02,  1.0524E-02,    F29170
     * 1.0518E-02,  1.0227E-02,  9.6837E-03,  9.6425E-03,  8.9938E-03,    F29180
     * 9.1488E-03,  8.8595E-03,  8.5976E-03,  8.4447E-03,  8.0731E-03,    F29190
     * 8.0283E-03,  7.7827E-03,  7.7638E-03,  7.2438E-03,  6.8246E-03,    F29200
     * 6.3457E-03,  5.6632E-03,  5.2500E-03,  4.3593E-03,  3.9431E-03,    F29210
     * 3.1580E-03,  2.2298E-03,  1.7818E-03,  1.4513E-03,  1.3188E-03/    F29220
      DATA C10407/                                                        F29230
     * 2.1034E-03/                                                        F29240
      DATA C10001 /                                                       F29250
     * 6.45621E-03, 7.11308E-03, 1.06130E-02, 1.36338E-02, 1.27746E-02,   F29260
     * 1.42152E-02, 1.41144E-02, 1.64830E-02, 1.67110E-02, 1.57368E-02,   F29270
     * 1.54644E-02, 1.45248E-02, 1.43206E-02, 1.56946E-02, 1.54268E-02,   F29280
     * 1.37500E-02, 1.50224E-02, 1.60919E-02, 1.49099E-02, 1.53960E-02,   F29290
     * 1.61871E-02, 1.46539E-02, 1.38258E-02, 1.32571E-02, 1.21580E-02,   F29300
     * 1.39596E-02, 1.16029E-02, 1.47042E-02, 1.07441E-02, 1.08999E-02,   F29310
     * 1.05562E-02, 1.00589E-02, 9.60711E-03, 9.36950E-03, 7.65303E-03,   F29320
     * 6.86216E-03, 7.05344E-03, 6.90728E-03, 6.78627E-03, 6.97435E-03,   F29330
     * 5.75456E-03, 5.81685E-03, 5.00915E-03, 4.90259E-03, 4.42545E-03,   F29340
     * 4.14633E-03, 3.61657E-03, 3.08178E-03, 2.91680E-03, 2.94554E-03,   F29350
     * 3.35794E-03, 5.49025E-03, 7.09867E-03, 6.82592E-03, 8.84835E-03,   F29360
     * 9.15718E-03, 9.17935E-03, 8.31848E-03, 7.79481E-03, 7.75125E-03,   F29370
     * 6.95844E-03, 7.34506E-03, 7.53823E-03, 7.03272E-03, 7.57051E-03,   F29380
     * 9.20239E-03, 1.10864E-02, 1.16188E-02, 1.30029E-02, 1.44364E-02,   F29390
     * 1.29292E-02, 1.36031E-02, 1.35967E-02, 1.30412E-02, 1.29874E-02,   F29400
     * 1.14829E-02, 1.18009E-02, 1.20829E-02, 1.17831E-02, 1.21489E-02/   F29410
      DATA C10081 /                                                       F29420
     * 1.27019E-02, 1.25557E-02, 1.23812E-02, 1.20158E-02, 1.26749E-02,   F29430
     * 1.17139E-02, 1.14552E-02, 1.11268E-02, 9.79143E-03, 8.79741E-03,   F29440
     * 8.85709E-03, 8.57653E-03, 8.93908E-03, 8.46205E-03, 8.56506E-03,   F29450
     * 8.14319E-03, 8.14415E-03, 7.74205E-03, 7.80727E-03, 7.49886E-03,   F29460
     * 7.71114E-03, 6.55963E-03, 6.87550E-03, 6.39162E-03, 5.55359E-03,   F29470
     * 5.43275E-03, 4.90649E-03, 4.41165E-03, 4.21875E-03, 3.62592E-03,   F29480
     * 3.40700E-03, 2.40267E-03, 2.61479E-03, 2.75677E-03, 4.10842E-03,   F29490
     * 5.79601E-03, 7.10708E-03, 8.07826E-03, 8.16166E-03, 8.72620E-03,   F29500
     * 8.85878E-03, 8.72755E-03, 8.25811E-03, 8.12100E-03, 7.78534E-03,   F29510
     * 7.39762E-03, 8.43880E-03, 8.53789E-03, 9.90072E-03, 1.01668E-02,   F29520
     * 1.00827E-02, 9.73556E-03, 9.57462E-03, 1.01289E-02, 1.10670E-02,   F29530
     * 1.03508E-02, 1.00929E-02, 9.10236E-03, 9.39459E-03, 8.79601E-03,   F29540
     * 8.67936E-03, 8.53862E-03, 7.95459E-03, 8.04037E-03, 7.95361E-03,   F29550
     * 7.87432E-03, 6.99165E-03, 7.37107E-03, 6.09187E-03, 6.21030E-03,   F29560
     * 5.33277E-03, 5.04633E-03, 4.45811E-03, 4.34153E-03, 3.98596E-03,   F29570
     * 3.84225E-03, 3.41943E-03, 3.60535E-03, 2.81691E-03, 2.49771E-03/   F29580
      DATA C10161 /                                                       F29590
     * 2.35046E-03, 2.50947E-03, 3.75462E-03, 4.92349E-03, 5.09294E-03,   F29600
     * 4.98312E-03, 5.19325E-03, 4.41827E-03, 4.25192E-03, 4.46745E-03,   F29610
     * 4.08731E-03, 3.84776E-03, 3.67507E-03, 3.76845E-03, 3.69210E-03,   F29620
     * 4.59864E-03, 6.42677E-03, 7.83255E-03, 7.89247E-03, 8.10883E-03,   F29630
     * 8.00825E-03, 8.40322E-03, 7.97108E-03, 8.24714E-03, 8.39006E-03,   F29640
     * 8.68787E-03, 8.61108E-03, 8.81552E-03, 9.36996E-03, 9.08243E-03,   F29650
     * 9.69116E-03, 9.66185E-03, 9.22856E-03, 9.65086E-03, 9.35398E-03,   F29660
     * 9.06358E-03, 8.76851E-03, 8.43072E-03, 7.85659E-03, 7.93936E-03,   F29670
     * 7.49712E-03, 7.20199E-03, 6.94581E-03, 6.64086E-03, 6.12627E-03,   F29680
     * 6.13967E-03, 5.67310E-03, 5.09928E-03, 4.59112E-03, 3.95257E-03,   F29690
     * 3.67652E-03, 3.28781E-03, 2.77471E-03, 2.74494E-03, 2.15529E-03,   F29700
     * 1.95283E-03, 1.75043E-03, 1.60419E-03, 1.82688E-03, 2.34667E-03,   F29710
     * 2.92502E-03, 3.88322E-03, 4.39984E-03, 4.67814E-03, 4.80395E-03,   F29720
     * 4.69130E-03, 4.54564E-03, 4.46773E-03, 4.59178E-03, 4.37498E-03,   F29730
     * 4.12706E-03, 4.18299E-03, 4.57267E-03, 5.60127E-03, 6.51936E-03,   F29740
     * 7.10498E-03, 7.49870E-03, 7.89554E-03, 7.97428E-03, 8.21044E-03/   F29750
      DATA C10241 /                                                       F29760
     * 8.06324E-03, 7.76648E-03, 7.62238E-03, 7.77675E-03, 7.46905E-03,   F29770
     * 7.61115E-03, 7.42715E-03, 7.28461E-03, 7.51514E-03, 7.38782E-03,   F29780
     * 6.97206E-03, 6.52738E-03, 6.10147E-03, 5.87553E-03, 5.49218E-03,   F29790
     * 4.94873E-03, 4.47920E-03, 4.25005E-03, 3.98094E-03, 3.92084E-03,   F29800
     * 3.41707E-03, 3.30501E-03, 3.09208E-03, 3.19686E-03, 3.55283E-03,   F29810
     * 4.20775E-03, 4.11155E-03, 3.72193E-03, 3.52000E-03, 3.13572E-03,   F29820
     * 2.87629E-03, 2.64251E-03, 2.33451E-03, 2.22426E-03, 2.05800E-03,   F29830
     * 1.75214E-03, 2.32530E-03, 2.68651E-03, 3.66315E-03, 4.93904E-03,   F29840
     * 5.32850E-03, 5.43978E-03, 5.32656E-03, 5.15649E-03, 5.42096E-03,   F29850
     * 5.37193E-03, 5.23454E-03, 5.34557E-03, 5.50533E-03, 6.13216E-03,   F29860
     * 6.65129E-03, 7.09357E-03, 7.46042E-03, 7.68605E-03, 7.91866E-03,   F29870
     * 7.52953E-03, 7.48272E-03, 7.17800E-03, 6.80060E-03, 6.60427E-03,   F29880
     * 6.43049E-03, 6.45975E-03, 6.20534E-03, 5.93094E-03, 5.67360E-03,   F29890
     * 5.38584E-03, 5.19364E-03, 4.92599E-03, 4.60655E-03, 4.24669E-03,   F29900
     * 3.94253E-03, 3.55894E-03, 3.24256E-03, 2.92974E-03, 2.62760E-03,   F29910
     * 2.52238E-03, 2.24714E-03, 2.26350E-03, 2.44380E-03, 3.03798E-03/   F29920
      DATA C10321 /                                                       F29930
     * 3.50000E-03, 3.55416E-03, 3.43661E-03, 3.19814E-03, 3.02155E-03,   F29940
     * 2.73890E-03, 2.50078E-03, 2.34595E-03, 2.18282E-03, 2.19285E-03,   F29950
     * 2.49482E-03, 3.13434E-03, 4.18947E-03, 4.72069E-03, 5.29712E-03,   F29960
     * 5.39004E-03, 5.44846E-03, 5.37952E-03, 5.09935E-03, 5.08741E-03,   F29970
     * 5.05257E-03, 5.10339E-03, 5.17968E-03, 5.31841E-03, 5.58106E-03,   F29980
     * 5.65031E-03, 5.65680E-03, 5.76184E-03, 5.71213E-03, 5.48515E-03,   F29990
     * 5.32168E-03, 5.18505E-03, 4.99640E-03, 4.78746E-03, 4.57244E-03,   F30000
     * 4.32728E-03, 4.14464E-03, 3.97659E-03, 4.01874E-03, 4.10588E-03,   F30010
     * 3.99644E-03, 3.84584E-03, 3.64222E-03, 3.39590E-03, 3.00386E-03,   F30020
     * 2.73790E-03, 2.45095E-03, 2.29068E-03, 1.64530E-03, 1.68602E-03,   F30030
     * 2.32934E-03, 3.14851E-03, 3.65706E-03, 3.70878E-03, 3.75103E-03,   F30040
     * 3.79183E-03, 3.32032E-03, 2.42604E-03, 2.48775E-03, 2.34603E-03,   F30050
     * 2.36349E-03, 3.33744E-03, 3.44617E-03, 4.27280E-03, 4.61076E-03,   F30060
     * 5.20165E-03, 5.14851E-03, 5.22909E-03, 5.08278E-03, 5.16125E-03,   F30070
     * 5.01572E-03, 4.51685E-03, 4.67541E-03, 4.83421E-03, 4.57546E-03,   F30080
     * 4.55111E-03, 5.03093E-03, 4.67838E-03, 4.44282E-03, 4.40774E-03/   F30090
      DATA C10401 /                                                       F30100
     * 4.48123E-03, 4.24410E-03, 4.03559E-03, 3.73969E-03, 3.45458E-03,   F30110
     * 3.18217E-03, 3.16115E-03, 3.36877E-03, 3.62026E-03, 3.69898E-03,   F30120
     * 3.49845E-03, 3.13839E-03, 2.77731E-03, 2.40106E-03, 2.03935E-03,   F30130
     * 1.84377E-03, 2.07757E-03, 2.39550E-03, 2.86272E-03, 3.27900E-03,   F30140
     * 3.42304E-03, 3.50211E-03, 3.29197E-03, 3.24784E-03, 3.20864E-03,   F30150
     * 3.28063E-03, 3.01328E-03, 3.00379E-03, 3.19562E-03, 3.45113E-03,   F30160
     * 3.75149E-03, 3.98520E-03, 4.19181E-03, 4.15773E-03, 4.02490E-03,   F30170
     * 3.95936E-03, 3.79001E-03, 3.77647E-03, 3.48528E-03, 3.55768E-03,   F30180
     * 3.62812E-03, 3.48650E-03, 3.35434E-03, 3.20088E-03, 3.25316E-03,   F30190
     * 3.04467E-03, 3.12633E-03, 3.23602E-03, 3.07723E-03, 2.80070E-03,   F30200
     * 2.72498E-03, 2.74752E-03, 2.58943E-03, 2.32482E-03, 2.20218E-03,   F30210
     * 2.10846E-03, 2.05991E-03, 2.01844E-03, 2.16224E-03, 2.48456E-03,   F30220
     * 2.88022E-03, 2.93939E-03, 3.01176E-03, 2.98886E-03, 2.96947E-03,   F30230
     * 3.38082E-03, 3.61657E-03, 3.42654E-03, 3.41274E-03, 3.22475E-03,   F30240
     * 2.97658E-03, 3.21944E-03, 3.32032E-03, 3.33273E-03, 3.58854E-03,   F30250
     * 3.67023E-03, 3.64069E-03, 3.74557E-03, 3.77703E-03, 3.64042E-03/   F30260
      DATA C10481 /                                                       F30270
     * 3.39468E-03, 3.22657E-03, 3.16466E-03, 3.24224E-03, 3.24801E-03,   F30280
     * 3.19487E-03, 3.40155E-03, 3.16940E-03, 2.92293E-03, 3.00998E-03,   F30290
     * 2.82851E-03, 2.60381E-03, 2.59242E-03, 2.48530E-03, 2.76677E-03,   F30300
     * 2.45506E-03, 2.21845E-03, 2.30407E-03, 2.28136E-03, 2.37278E-03,   F30310
     * 2.25313E-03, 2.47836E-03, 2.77858E-03, 2.89803E-03, 2.86131E-03,   F30320
     * 3.14118E-03, 3.14119E-03, 2.88881E-03, 3.19502E-03, 2.99538E-03,   F30330
     * 2.91212E-03, 3.22739E-03, 3.05960E-03, 3.18901E-03, 3.05805E-03,   F30340
     * 3.12205E-03, 2.95636E-03, 3.24111E-03, 3.29433E-03, 3.09206E-03,   F30350
     * 3.06696E-03, 2.97735E-03, 2.90897E-03, 2.88979E-03, 2.75105E-03,   F30360
     * 2.92156E-03, 3.03445E-03, 2.91664E-03, 2.85559E-03, 2.98405E-03,   F30370
     * 2.95376E-03, 2.80234E-03, 2.78349E-03, 2.73421E-03, 2.70035E-03,   F30380
     * 2.60074E-03, 2.34840E-03, 2.37626E-03, 2.32927E-03, 2.20842E-03,   F30390
     * 2.31080E-03, 2.42771E-03, 2.43339E-03, 2.53280E-03, 2.37093E-03,   F30400
     * 2.37377E-03, 2.73453E-03, 2.60836E-03, 2.55568E-03, 2.44062E-03,   F30410
     * 2.71093E-03, 2.64421E-03, 2.66969E-03, 2.55560E-03, 2.71800E-03,   F30420
     * 2.79534E-03, 2.59070E-03, 2.55373E-03, 2.45272E-03, 2.55571E-03/   F30430
      DATA C10561 /                                                       F30440
     * 2.54606E-03, 2.57349E-03, 2.46807E-03, 2.35634E-03, 2.44470E-03,   F30450
     * 2.47050E-03, 2.57131E-03, 2.71649E-03, 2.58800E-03, 2.54524E-03,   F30460
     * 2.69505E-03, 2.89122E-03, 2.77399E-03, 2.63306E-03, 2.82269E-03,   F30470
     * 2.95684E-03, 3.07415E-03, 2.70594E-03, 2.65650E-03, 2.90613E-03,   F30480
     * 2.96666E-03, 2.94767E-03, 2.81765E-03, 2.64829E-03, 2.43062E-03,   F30490
     * 2.33816E-03, 2.38210E-03, 2.45701E-03, 2.38508E-03, 2.40746E-03,   F30500
     * 2.49779E-03, 2.28209E-03, 2.26185E-03, 2.26604E-03, 2.19232E-03,   F30510
     * 2.19160E-03, 2.32246E-03, 2.11108E-03, 2.26220E-03, 2.26849E-03,   F30520
     * 2.34787E-03, 2.49323E-03, 2.46872E-03, 2.52974E-03, 2.35858E-03,   F30530
     * 2.36865E-03, 2.33533E-03, 2.21338E-03, 2.24610E-03, 2.24776E-03,   F30540
     * 2.24423E-03, 2.29276E-03, 2.18487E-03, 2.27621E-03, 2.31141E-03,   F30550
     * 2.44095E-03, 2.45198E-03, 2.56919E-03, 2.56823E-03, 2.41982E-03,   F30560
     * 2.39968E-03, 2.62447E-03, 2.55339E-03, 2.51556E-03, 2.47477E-03,   F30570
     * 2.50276E-03, 2.48381E-03, 2.48484E-03, 2.48316E-03, 2.38541E-03,   F30580
     * 2.41183E-03, 2.55888E-03, 2.42810E-03, 2.43356E-03, 2.25996E-03,   F30590
     * 2.34736E-03, 2.10305E-03, 2.13870E-03, 2.17472E-03, 2.05354E-03/   F30600
      DATA C10641 /                                                       F30610
     * 2.11572E-03, 2.19557E-03, 2.09545E-03, 2.07831E-03, 1.94425E-03,   F30620
     * 1.89333E-03, 1.98025E-03, 1.98328E-03, 2.01702E-03, 1.98333E-03,   F30630
     * 2.01150E-03, 2.02484E-03, 2.10759E-03, 2.11892E-03, 2.10175E-03,   F30640
     * 2.05314E-03, 2.13338E-03, 2.25764E-03, 2.19055E-03, 2.10818E-03,   F30650
     * 2.05100E-03, 2.05685E-03, 2.10843E-03, 2.10228E-03, 2.10646E-03,   F30660
     * 2.22640E-03, 2.31253E-03, 2.31230E-03, 2.21885E-03, 2.19568E-03,   F30670
     * 2.23583E-03, 2.34754E-03, 2.28622E-03, 2.21876E-03, 2.26679E-03,   F30680
     * 2.30828E-03, 2.24944E-03, 2.13851E-03, 2.02938E-03, 1.96770E-03,   F30690
     * 2.05953E-03, 2.13814E-03, 2.03158E-03, 2.24655E-03, 1.95119E-03,   F30700
     * 2.12979E-03, 2.08581E-03, 2.02434E-03, 1.98926E-03, 1.98792E-03,   F30710
     * 1.97237E-03, 1.93397E-03, 1.92360E-03, 1.90805E-03, 1.89300E-03,   F30720
     * 1.83548E-03, 1.87215E-03, 1.85589E-03, 1.85718E-03, 1.79361E-03,   F30730
     * 1.77984E-03, 1.91506E-03, 2.04256E-03, 2.04095E-03, 1.94031E-03,   F30740
     * 1.90447E-03, 2.02049E-03, 1.98360E-03, 2.04364E-03, 2.02519E-03,   F30750
     * 2.20802E-03, 1.96964E-03, 1.94559E-03, 2.09922E-03, 2.11184E-03,   F30760
     * 2.05706E-03, 2.02257E-03, 2.01781E-03, 2.01055E-03, 1.86538E-03/   F30770
      DATA C10721 /                                                       F30780
     * 1.86899E-03, 1.76798E-03, 1.85871E-03, 1.95363E-03, 1.96404E-03,   F30790
     * 1.84169E-03, 1.82851E-03, 1.84582E-03, 1.81997E-03, 1.76461E-03,   F30800
     * 1.68384E-03, 1.65530E-03, 1.73550E-03, 1.62463E-03, 1.68793E-03,   F30810
     * 1.60472E-03, 1.67560E-03, 1.67431E-03, 1.61779E-03, 1.66446E-03,   F30820
     * 1.66403E-03, 1.55724E-03, 1.62351E-03, 1.71545E-03, 1.69645E-03,   F30830
     * 1.59540E-03, 1.62948E-03, 1.66784E-03, 1.66416E-03, 1.66131E-03,   F30840
     * 1.71502E-03, 1.76555E-03, 1.75182E-03, 1.72327E-03, 1.72338E-03,   F30850
     * 1.69993E-03, 1.78819E-03, 1.73517E-03, 1.74802E-03, 1.81751E-03,   F30860
     * 1.70973E-03, 1.65075E-03, 1.70784E-03, 1.73655E-03, 1.71670E-03,   F30870
     * 1.67367E-03, 1.69338E-03, 1.61772E-03, 1.54914E-03, 1.56009E-03,   F30880
     * 1.59467E-03, 1.60761E-03, 1.57117E-03, 1.54045E-03, 1.53102E-03,   F30890
     * 1.44516E-03, 1.49898E-03, 1.56048E-03, 1.60087E-03, 1.62636E-03,   F30900
     * 1.62472E-03, 1.53931E-03, 1.55536E-03, 1.61649E-03, 1.66493E-03,   F30910
     * 1.86915E-03, 1.59984E-03, 1.60483E-03, 1.66549E-03, 1.73449E-03,   F30920
     * 1.73673E-03, 1.68393E-03, 1.67434E-03, 1.77880E-03, 1.76154E-03,   F30930
     * 1.43028E-03, 1.69651E-03, 1.60934E-03, 1.69413E-03, 1.70514E-03/   F30940
      DATA C10801 /                                                       F30950
     * 1.62471E-03, 1.74854E-03, 1.76480E-03, 1.63495E-03, 1.59364E-03,   F30960
     * 1.39603E-03, 1.47897E-03, 1.49509E-03, 1.70002E-03, 1.63048E-03,   F30970
     * 1.44807E-03, 1.45071E-03, 1.53998E-03, 1.45276E-03, 1.29129E-03,   F30980
     * 1.52900E-03, 1.64444E-03, 1.37450E-03, 1.42574E-03, 1.47355E-03,   F30990
     * 1.51202E-03, 1.54376E-03, 1.51421E-03, 1.43989E-03, 1.45732E-03,   F31000
     * 1.42912E-03, 1.59906E-03, 1.56748E-03, 1.52383E-03, 1.47665E-03,   F31010
     * 1.51465E-03, 1.55582E-03, 1.54521E-03, 1.55189E-03, 1.56772E-03,   F31020
     * 1.45401E-03, 1.55775E-03, 1.43120E-03, 1.39659E-03, 1.41451E-03,   F31030
     * 1.45157E-03, 1.48303E-03, 1.42540E-03, 1.26387E-03, 1.37479E-03,   F31040
     * 1.46381E-03, 1.38134E-03, 1.32733E-03, 1.38030E-03, 1.44619E-03,   F31050
     * 1.41344E-03, 1.31982E-03, 1.24944E-03, 1.20096E-03, 1.21107E-03,   F31060
     * 1.27999E-03, 1.22523E-03, 1.22193E-03, 1.35957E-03, 1.41427E-03,   F31070
     * 1.35679E-03, 1.15438E-03, 1.41184E-03, 1.49093E-03, 1.32193E-03,   F31080
     * 1.25009E-03, 1.37625E-03, 1.49022E-03, 1.44180E-03, 1.27628E-03,   F31090
     * 1.29670E-03, 1.31636E-03, 1.28874E-03, 1.31177E-03, 1.35732E-03,   F31100
     * 1.33854E-03, 1.30253E-03, 1.31374E-03, 1.27379E-03, 1.18339E-03/   F31110
      DATA C10881 /                                                       F31120
     * 1.22016E-03, 1.26551E-03, 1.26371E-03, 1.28180E-03, 1.36024E-03,   F31130
     * 1.45759E-03, 1.29413E-03, 1.35858E-03, 1.26528E-03, 1.18623E-03,   F31140
     * 1.21812E-03, 1.28799E-03, 1.37028E-03, 1.29268E-03, 1.27639E-03,   F31150
     * 1.19487E-03, 1.23542E-03, 1.25010E-03, 1.17418E-03, 1.13914E-03,   F31160
     * 1.21951E-03, 1.13780E-03, 1.16443E-03, 1.17883E-03, 1.11982E-03,   F31170
     * 1.05708E-03, 1.04865E-03, 1.05884E-03, 1.06599E-03, 1.13828E-03,   F31180
     * 1.10373E-03, 1.07739E-03, 1.04632E-03, 1.06118E-03, 1.15445E-03,   F31190
     * 1.17300E-03, 1.00675E-03, 1.04235E-03, 1.08398E-03, 1.06587E-03,   F31200
     * 1.05536E-03, 1.08614E-03, 1.09026E-03, 1.09141E-03, 1.13051E-03,   F31210
     * 1.08667E-03, 1.04016E-03, 1.04897E-03, 1.08894E-03, 1.09682E-03,   F31220
     * 1.09638E-03, 9.79254E-04, 1.00668E-03, 1.02569E-03, 1.00581E-03,   F31230
     * 9.74433E-04, 9.66321E-04, 9.78440E-04, 9.01587E-04, 1.02149E-03,   F31240
     * 9.87464E-04, 9.41872E-04, 9.05021E-04, 8.59547E-04, 9.03963E-04,   F31250
     * 8.66415E-04, 8.84726E-04, 8.77087E-04, 8.70584E-04, 8.81338E-04,   F31260
     * 8.97658E-04, 8.97586E-04, 9.19028E-04, 8.82438E-04, 9.00710E-04,   F31270
     * 9.54329E-04, 9.54490E-04, 9.10940E-04, 9.95472E-04, 9.50134E-04/   F31280
      DATA C10961 /                                                       F31290
     * 9.17127E-04, 9.70916E-04, 9.87575E-04, 9.65026E-04, 9.71779E-04,   F31300
     * 1.00967E-03, 1.00053E-03, 9.26063E-04, 9.34721E-04, 9.76354E-04,   F31310
     * 9.78436E-04, 9.36012E-04, 9.64448E-04, 9.95903E-04, 9.89960E-04,   F31320
     * 9.41143E-04, 9.04393E-04, 8.84719E-04, 8.41396E-04, 8.67234E-04,   F31330
     * 8.55864E-04, 8.63314E-04, 8.72317E-04, 8.40899E-04, 7.79593E-04,   F31340
     * 7.88481E-04, 8.21075E-04, 7.38342E-04, 7.56537E-04, 7.57278E-04,   F31350
     * 7.35854E-04, 7.32765E-04, 6.67398E-04, 7.45338E-04, 7.33094E-04,   F31360
     * 7.01840E-04, 6.85595E-04, 6.95740E-04, 7.24015E-04, 7.00907E-04,   F31370
     * 7.28498E-04, 6.89410E-04, 6.91728E-04, 7.40601E-04, 7.62775E-04,   F31380
     * 7.40912E-04, 7.35021E-04, 7.07799E-04, 7.54113E-04, 8.44845E-04,   F31390
     * 8.53956E-04, 6.42186E-04, 7.40557E-04, 7.54340E-04, 7.55544E-04,   F31400
     * 7.88986E-04, 7.97902E-04, 6.98460E-04, 7.74873E-04, 6.81178E-04,   F31410
     * 7.15567E-04, 7.56723E-04, 7.98438E-04, 8.83150E-04, 8.45671E-04,   F31420
     * 7.40924E-04, 7.35498E-04, 7.77829E-04, 6.93566E-04, 5.10188E-04,   F31430
     * 7.52717E-04, 6.94185E-04, 6.71928E-04, 6.73286E-04, 6.89415E-04,   F31440
     * 7.22917E-04, 7.89448E-04, 8.53812E-04, 7.45132E-04, 7.68732E-04/   F31450
      DATA C11041 /                                                       F31460
     * 8.10104E-04, 7.55615E-04, 7.09145E-04, 6.80676E-04, 7.54594E-04,   F31470
     * 7.89416E-04, 7.88579E-04, 7.49805E-04, 6.13534E-04, 7.22491E-04,   F31480
     * 7.95410E-04, 7.80604E-04, 7.74283E-04, 7.93224E-04, 6.86522E-04,   F31490
     * 8.06038E-04, 8.30285E-04, 8.37763E-04, 8.03863E-04, 7.33526E-04,   F31500
     * 7.42588E-04, 6.31046E-04, 8.16153E-04, 8.95391E-04, 8.61330E-04,   F31510
     * 8.38726E-04, 8.16761E-04, 8.16118E-04, 6.37058E-04, 6.30868E-04,   F31520
     * 7.26410E-04, 7.03464E-04, 5.93454E-04, 6.01985E-04, 6.51157E-04,   F31530
     * 6.68569E-04, 6.56297E-04, 6.58732E-04, 5.99721E-04, 5.34301E-04,   F31540
     * 5.33271E-04, 5.57992E-04, 5.70096E-04, 5.59932E-04, 5.32110E-04,   F31550
     * 5.64713E-04, 6.25026E-04, 6.38973E-04, 6.05323E-04, 7.17460E-04,   F31560
     * 6.19407E-04, 5.90228E-04, 5.43682E-04, 5.38446E-04, 6.56146E-04,   F31570
     * 6.09081E-04, 6.04737E-04, 6.45526E-04, 6.46978E-04, 5.89738E-04,   F31580
     * 5.63852E-04, 6.18018E-04, 5.71768E-04, 5.75433E-04, 6.05766E-04,   F31590
     * 5.93065E-04, 5.31708E-04, 5.41187E-04, 5.76985E-04, 5.78176E-04,   F31600
     * 5.75339E-04, 6.85426E-04, 5.51038E-04, 6.02049E-04, 6.20406E-04,   F31610
     * 5.80169E-04, 5.36399E-04, 5.59608E-04, 5.46575E-04, 5.66979E-04/   F31620
      DATA C11121 /                                                       F31630
     * 5.94982E-04, 6.18469E-04, 6.56281E-04, 8.22124E-04, 7.81716E-04,   F31640
     * 7.29616E-04, 7.14460E-04, 7.08969E-04, 6.53794E-04, 7.33138E-04,   F31650
     * 8.29513E-04, 8.99395E-04, 9.05526E-04, 7.98257E-04, 7.86935E-04,   F31660
     * 6.10797E-04, 4.63912E-04, 4.05675E-04, 3.66230E-04, 4.86472E-04,   F31670
     * 5.31818E-04, 5.15865E-04, 4.87344E-04, 4.99857E-04, 5.35479E-04,   F31680
     * 5.27561E-04, 4.99000E-04, 4.77056E-04, 4.74242E-04, 4.66595E-04,   F31690
     * 4.66325E-04, 4.94704E-04, 5.12842E-04, 5.01795E-04, 4.80789E-04,   F31700
     * 5.73709E-04, 5.65214E-04, 5.11321E-04, 4.55242E-04, 4.29330E-04,   F31710
     * 5.09792E-04, 4.70489E-04, 4.82859E-04, 4.99195E-04, 4.07724E-04,   F31720
     * 4.99951E-04, 4.55755E-04, 4.42528E-04, 4.19433E-04, 3.31325E-04,   F31730
     * 3.70517E-04, 3.77708E-04, 2.97923E-04, 2.27470E-04, 2.47389E-04,   F31740
     * 2.38324E-04, 2.56706E-04, 2.45046E-04, 2.62539E-04, 3.37054E-04,   F31750
     * 3.33930E-04, 3.01390E-04, 3.08028E-04, 3.41464E-04, 3.70574E-04,   F31760
     * 3.47893E-04, 3.28433E-04, 3.46976E-04, 3.60351E-04, 3.50559E-04,   F31770
     * 3.56070E-04, 3.62782E-04, 3.37330E-04, 3.33763E-04, 3.57046E-04,   F31780
     * 3.08784E-04, 2.93898E-04, 2.80842E-04, 2.54114E-04, 2.38198E-04/   F31790
      DATA C11201 /                                                       F31800
     * 3.48753E-04, 2.97334E-04, 2.82929E-04, 2.94150E-04, 3.07875E-04,   F31810
     * 3.21129E-04, 3.38335E-04, 3.49826E-04, 3.47647E-04, 3.35438E-04,   F31820
     * 3.58145E-04, 3.72391E-04, 3.59372E-04, 3.64755E-04, 4.16867E-04,   F31830
     * 3.43614E-04, 3.34932E-04, 3.12782E-04, 3.28220E-04, 4.32595E-04,   F31840
     * 3.49513E-04, 3.51861E-04, 3.81166E-04, 3.91194E-04, 3.38944E-04,   F31850
     * 2.63445E-04, 2.49520E-04, 2.46184E-04, 2.33203E-04, 2.16315E-04,   F31860
     * 1.89536E-04, 1.95730E-04, 1.99664E-04, 1.77139E-04, 1.27969E-04,   F31870
     * 5.17216E-05, 7.60445E-05, 1.24418E-04, 1.30989E-04, 2.31539E-04,   F31880
     * 2.21334E-04, 2.08757E-04, 2.18351E-04, 2.46202E-04, 2.29824E-04,   F31890
     * 2.28909E-04, 2.88826E-04, 3.58039E-04, 2.60800E-04, 2.33025E-04,   F31900
     * 2.52667E-04, 2.61394E-04, 2.31384E-04, 2.29388E-04, 2.54701E-04,   F31910
     * 2.21158E-04, 1.61506E-04, 1.36752E-04, 1.69481E-04, 8.64539E-05,   F31920
     * 1.64407E-04, 3.65674E-04, 3.18233E-04, 4.00755E-04, 3.33375E-04,   F31930
     * 2.62930E-04, 2.87052E-04, 2.51395E-04, 2.85274E-04, 2.66915E-04,   F31940
     * 2.10866E-04, 1.89517E-04, 1.67378E-04, 2.79951E-04, 2.97224E-04,   F31950
     * 1.89222E-04, 3.33825E-04, 3.56386E-04, 3.89727E-04, 4.30407E-04/   F31960
      DATA C11281 /                                                       F31970
     * 4.45922E-04, 4.23446E-04, 4.41347E-04, 4.06723E-04, 3.00181E-04,   F31980
     * 1.85243E-04, 3.13176E-04, 4.08991E-04, 4.24776E-04, 3.56412E-04,   F31990
     * 3.84760E-04, 2.30602E-04, 1.77702E-04, 2.62329E-04, 2.49442E-04,   F32000
     * 3.76212E-04, 3.69176E-04, 2.97681E-04, 2.71662E-04, 2.05694E-04,   F32010
     * 2.11418E-04, 2.25439E-04, 2.27013E-04, 2.47845E-04, 3.14603E-04,   F32020
     * 2.68802E-04, 2.04334E-04, 2.77399E-04, 2.68273E-04, 2.04991E-04,   F32030
     * 2.24441E-04, 3.55074E-04, 2.90135E-04, 3.35680E-04, 3.59358E-04,   F32040
     * 3.44716E-04, 3.24496E-04, 3.48146E-04, 3.49042E-04, 3.54848E-04,   F32050
     * 3.86418E-04, 3.59198E-04, 3.47608E-04, 3.20522E-04, 2.78401E-04,   F32060
     * 2.64579E-04, 2.23694E-04, 2.34370E-04, 2.52559E-04, 1.88475E-04,   F32070
     * 2.01258E-04, 1.63979E-04, 1.45384E-04, 1.91215E-04, 1.76958E-04,   F32080
     * 1.69167E-04, 1.71767E-04, 1.86595E-04, 2.14969E-04, 2.48345E-04,   F32090
     * 2.46691E-04, 2.25234E-04, 2.26755E-04, 1.64112E-04, 1.87750E-04,   F32100
     * 2.22984E-04, 2.00443E-04, 2.38863E-04, 2.77590E-04, 2.91953E-04,   F32110
     * 2.80611E-04, 3.08215E-04, 1.79095E-04, 1.46920E-04, 2.29177E-04,   F32120
     * 2.54685E-04, 2.68866E-04, 2.13346E-04, 1.20122E-04, 5.55240E-05/   F32130
      DATA C11361 /                                                       F32140
     * 5.99017E-05, 1.07768E-04, 1.67810E-04, 2.06886E-04, 2.36232E-04,   F32150
     * 2.24598E-04, 2.30792E-04, 2.71274E-04, 1.29062E-04, 1.92624E-04,   F32160
     * 2.38438E-04, 1.98994E-04, 1.81687E-04, 2.55733E-04, 2.84379E-04,   F32170
     * 2.54459E-04, 2.30884E-04, 2.68873E-04, 3.07231E-04, 3.15063E-04,   F32180
     * 2.46725E-04, 2.60370E-04, 2.66391E-04, 2.50708E-04, 2.04296E-04,   F32190
     * 1.66011E-04, 1.19164E-04, 1.06700E-04, 1.77576E-04, 1.91741E-04,   F32200
     * 1.66618E-04, 1.49824E-04, 1.80699E-04, 2.20905E-04, 1.38754E-04,   F32210
     * 6.27971E-05, 7.52567E-05, 1.89995E-04, 1.72489E-04, 1.40424E-04,   F32220
     * 1.52384E-04, 1.63942E-04, 1.19901E-04, 1.49234E-04, 2.68313E-04,   F32230
     * 2.08815E-04, 1.17218E-04, 1.42235E-04, 2.71237E-04, 1.38192E-04,   F32240
     * 2.15643E-04, 2.84476E-04, 2.78117E-04, 2.19234E-04, 1.59128E-04,   F32250
     * 1.78819E-04, 2.67785E-04, 2.66786E-04, 2.58545E-04, 2.68476E-04,   F32260
     * 2.88542E-04, 2.59726E-04, 3.00936E-04, 3.11237E-04, 2.61275E-04,   F32270
     * 1.37136E-04, 2.76566E-04, 3.82888E-04, 3.97564E-04, 4.43655E-04,   F32280
     * 3.15415E-04, 2.60869E-04, 3.19171E-04, 3.34205E-04, 2.02914E-04,   F32290
     * 1.16223E-04, 1.14737E-04, 6.10978E-05,-8.03695E-06,-1.07062E-05/   F32300
      DATA C11441 /                                                       F32310
     * 6.50664E-05, 1.12586E-04, 1.56727E-04, 1.57927E-04, 1.05762E-04,   F32320
     * 1.03646E-04, 1.72520E-04, 2.23668E-04, 2.12775E-04, 2.33525E-04,   F32330
     * 2.75558E-04, 2.34256E-04, 5.10062E-05, 1.76007E-04, 1.70850E-04,   F32340
     * 1.43266E-04, 1.89626E-04, 2.97283E-04, 3.02773E-04, 2.74401E-04,   F32350
     * 3.00754E-04, 3.66813E-04, 3.54383E-04, 2.90580E-04, 2.32206E-04,   F32360
     * 1.58405E-04, 1.54663E-04, 1.84598E-04, 1.26408E-04, 2.14481E-04,   F32370
     * 2.00791E-04, 1.05796E-04, 2.39794E-04, 1.66105E-04, 7.88615E-05,   F32380
     * 4.30615E-05, 7.37518E-05, 1.24926E-04, 1.38295E-04, 8.54356E-05,   F32390
     * 6.12641E-05, 6.54466E-05, 6.17727E-05, 1.30688E-05, 6.00462E-05,   F32400
     * 1.52612E-04, 2.11656E-04, 9.67692E-05, 8.67858E-05, 1.34888E-04,   F32410
     * 1.90899E-04, 1.03234E-04, 1.03837E-04, 1.49767E-04, 2.19058E-04,   F32420
     * 2.26549E-04, 2.11506E-04, 1.85238E-04, 1.53774E-04, 1.32313E-04,   F32430
     * 6.10658E-05, 2.37782E-05, 1.24450E-04, 1.87610E-04, 1.44775E-04,   F32440
     * 5.60937E-05, 6.64032E-05, 1.28073E-04, 1.77512E-04, 1.84684E-04,   F32450
     * 5.73677E-05, 5.29679E-05, 9.95510E-05, 1.61423E-04, 3.19036E-04,   F32460
     * 3.17383E-04, 2.36505E-04, 1.80844E-04, 1.63722E-04, 1.21478E-04/   F32470
      DATA C11521 /                                                       F32480
     * 6.85823E-05, 7.42058E-05, 1.14838E-04, 1.21131E-04, 8.01009E-05,   F32490
     * 1.52058E-04, 2.18368E-04, 2.53416E-04, 2.27116E-04, 1.25336E-04,   F32500
     * 6.26421E-05, 5.32471E-05, 1.34705E-04, 2.07005E-05,-5.18630E-05,   F32510
     *-3.25696E-05,-8.06171E-05,-1.09430E-04,-1.05637E-04,-4.96066E-05,   F32520
     *-7.76138E-05,-4.85930E-05, 3.65111E-06,-2.86933E-05,-4.61366E-05,   F32530
     *-4.88820E-05,-3.08816E-05, 8.43778E-05, 1.40484E-04, 1.31125E-04,   F32540
     * 3.55198E-05, 8.47412E-05, 1.23408E-04, 1.36799E-04, 1.21147E-04,   F32550
     * 1.25585E-04, 1.32337E-04, 1.34092E-04, 1.26652E-04, 1.12131E-04,   F32560
     * 1.00927E-04, 1.13828E-04, 1.06053E-04, 9.43643E-05, 8.33628E-05,   F32570
     * 8.65842E-05, 7.59315E-05, 8.28623E-05, 1.39681E-04, 1.80492E-04,   F32580
     * 1.65779E-04, 1.03843E-04, 3.10284E-05, 1.94408E-05, 4.57525E-05,   F32590
     * 1.02436E-04, 1.39750E-04, 1.43342E-04, 1.11999E-04, 2.94197E-05,   F32600
     * 2.76980E-05, 5.51685E-05, 9.39909E-05, 1.16108E-04, 7.72703E-05,   F32610
     * 4.37409E-05, 1.13925E-04, 8.18872E-05, 2.87657E-05,-2.41413E-05,   F32620
     * 1.24699E-05, 2.19589E-05,-5.88247E-06,-9.66151E-05,-2.06255E-05,   F32630
     *-1.83148E-06,-5.63625E-05,-8.65590E-05,-8.26020E-05,-5.06239E-05/   F32640
      DATA C11601 /                                                       F32650
     * 1.28065E-05,-1.34669E-05, 1.59701E-05, 9.44755E-05, 1.63032E-05,   F32660
     * 2.51304E-05, 7.38226E-05, 1.28405E-04, 1.17413E-04, 9.92387E-05,   F32670
     * 9.51533E-05, 2.17008E-04, 2.25854E-04, 1.90448E-04, 1.77207E-04,   F32680
     * 1.80844E-04, 1.53501E-04, 9.80430E-05, 1.27404E-04, 1.16465E-04,   F32690
     * 9.98611E-05, 1.25556E-04, 1.73627E-04, 1.12347E-04,-7.73523E-05,   F32700
     * 5.66599E-05, 5.36347E-05, 1.20227E-06, 6.96325E-05, 4.79010E-05,   F32710
     *-1.09886E-05,-9.16457E-05,-7.09170E-05,-5.31410E-05,-2.68376E-05,   F32720
     * 6.32641E-05, 8.06052E-06,-4.99262E-05,-2.56644E-05,-8.76854E-05,   F32730
     *-8.21360E-05,-5.02403E-06, 4.66629E-05, 6.93127E-05, 5.53828E-05,   F32740
     *-2.32399E-05,-2.07514E-05,-7.33240E-05,-2.10483E-04,-1.53757E-04,   F32750
     *-7.13861E-05,-1.07356E-05,-1.26578E-04,-7.48854E-05, 3.25418E-06,   F32760
     * 2.97068E-05, 3.35685E-05, 3.15022E-05, 2.68904E-05, 3.87401E-05,   F32770
     * 5.12522E-05, 5.12172E-05, 1.05053E-05, 1.65321E-05, 3.47537E-05,   F32780
     * 5.62503E-05, 4.18666E-05, 3.13970E-05, 3.11750E-05, 7.21547E-05,   F32790
     * 2.55262E-05,-2.76061E-05, 5.43449E-06,-5.20575E-05,-1.08627E-04,   F32800
     *-1.40475E-04,-1.59926E-04,-1.32237E-04,-8.15458E-05,-1.31738E-04/   F32810
      DATA C11681 /                                                       F32820
     *-1.64036E-04,-1.69351E-04,-1.24797E-04,-1.61950E-04,-2.01904E-04,   F32830
     *-2.22995E-04,-1.87647E-04,-1.70817E-04,-1.64583E-04,-1.12811E-04,   F32840
     *-8.38306E-05,-8.62707E-05,-1.54362E-04,-1.98090E-04,-2.12920E-04,   F32850
     *-1.89358E-04,-2.02988E-04,-1.72791E-04,-1.02863E-04,-1.09877E-04,   F32860
     *-1.04257E-04,-8.20734E-05,-2.18346E-05,-2.94593E-05,-4.18226E-05,   F32870
     *-1.86891E-05,-6.14620E-05,-3.21912E-05, 1.00844E-04, 6.92419E-05,   F32880
     * 3.16713E-05, 5.62042E-07, 5.18900E-05, 7.48835E-05, 8.03381E-05,   F32890
     * 7.24685E-05, 9.55588E-05, 9.22801E-05, 2.87159E-05, 2.26234E-05,   F32900
     * 2.62790E-05, 3.58332E-05, 6.23297E-05, 5.01998E-05, 1.81446E-05,   F32910
     * 3.33564E-05, 3.97765E-06,-2.60624E-05, 7.01802E-06,-4.16797E-05,   F32920
     *-8.70108E-05,-8.22182E-05,-6.64886E-05,-7.88704E-05,-1.28305E-04,   F32930
     *-1.29990E-04,-1.12646E-04,-8.68394E-05,-1.29584E-04,-1.44352E-04,   F32940
     *-1.42082E-04,-1.33790E-04,-1.27963E-04,-1.21233E-04,-1.09965E-04,   F32950
     *-1.02233E-04,-1.03804E-04,-1.19503E-04,-7.74707E-05,-4.66805E-05,   F32960
     *-3.52201E-05,-4.07406E-05,-4.66887E-05,-5.05962E-05,-3.30333E-05,   F32970
     *-3.47981E-05,-3.60962E-05, 1.44242E-05, 4.10478E-05, 3.68984E-05/   F32980
      DATA C11761 /                                                       F32990
     *-2.81300E-05, 2.83171E-05, 7.48062E-05, 4.29333E-05, 8.50076E-06,   F33000
     * 4.98135E-06, 4.44854E-05, 2.51860E-05, 3.12189E-05, 6.39424E-05,   F33010
     * 7.20715E-05, 9.89688E-05, 1.33768E-04, 1.07781E-04, 9.76731E-05,   F33020
     * 9.21479E-05, 6.72624E-05, 5.41295E-05, 4.89022E-05, 5.28039E-05,   F33030
     *-4.48737E-06,-5.15409E-05,-3.57396E-05,-1.94752E-05,-2.09521E-05,   F33040
     *-5.13096E-05,-2.62781E-05,-2.75451E-05,-6.98423E-05,-1.25462E-04,   F33050
     *-1.68362E-04,-1.97456E-04,-1.90669E-04,-2.06890E-04,-2.36699E-04,   F33060
     *-1.97732E-04,-1.76504E-04,-1.67505E-04,-1.60694E-04,-1.85851E-04,   F33070
     *-2.01567E-04,-9.82507E-05,-1.33338E-04,-1.95199E-04,-1.40781E-04,   F33080
     *-8.90988E-05,-3.63239E-05, 2.16510E-05,-1.56807E-05,-4.21285E-05,   F33090
     * 5.50505E-06, 6.78937E-07, 3.12346E-06, 3.64202E-05, 3.50651E-05,   F33100
     * 6.20423E-05, 1.38667E-04, 7.74738E-05, 6.77036E-05, 1.38367E-04,   F33110
     * 1.17359E-04, 1.06637E-04, 1.12404E-04, 9.78586E-05, 1.03178E-04,   F33120
     * 1.28717E-04, 1.56642E-04, 1.62544E-04, 1.50109E-04, 1.43214E-04,   F33130
     * 1.33651E-04, 1.24352E-04, 1.41420E-04, 1.36340E-04, 1.18769E-04,   F33140
     * 1.31656E-04, 8.81533E-05, 1.55214E-05,-3.68736E-07,-1.76213E-05/   F33150
      DATA C11841 /                                                       F33160
     *-2.85341E-05, 4.65155E-06, 5.41350E-06,-7.01247E-06, 6.57918E-06,   F33170
     *-2.45784E-05,-6.89104E-05,-6.90953E-05,-6.23937E-05,-6.72978E-05,   F33180
     *-1.39547E-04,-1.44228E-04,-1.42543E-04,-2.31080E-04,-2.12756E-04,   F33190
     *-1.62089E-04,-1.66063E-04,-1.61872E-04,-1.59764E-04,-1.80217E-04,   F33200
     *-1.38355E-04,-8.45661E-05,-7.58308E-05,-4.65144E-05,-2.76855E-05,   F33210
     *-7.48714E-05,-8.28561E-05,-6.45277E-05,-7.08509E-06,-1.05566E-05,   F33220
     *-1.96352E-05, 3.55561E-05, 2.24676E-05,-1.25648E-05,-1.87661E-05,   F33230
     * 6.99061E-06, 2.33676E-05,-5.25111E-05,-3.86758E-05, 1.03585E-06,   F33240
     *-1.65901E-05,-1.04855E-05, 5.03694E-06, 1.25937E-05,-8.31340E-06,   F33250
     *-4.37906E-05,-7.91444E-05,-4.62167E-05, 5.14238E-06,-4.52863E-05,   F33260
     *-5.86455E-05,-4.98093E-05,-3.03495E-05,-5.09377E-05,-8.88116E-05,   F33270
     *-6.21360E-05,-7.38148E-05,-1.07502E-04,-7.55276E-05,-6.39257E-05,   F33280
     *-6.86921E-05,-8.05504E-05,-9.24178E-05,-1.03991E-04,-1.00468E-04,   F33290
     *-6.71447E-05,-3.84897E-06,-5.99067E-06,-2.21894E-05,-5.21766E-05,   F33300
     *-3.93796E-05,-4.06712E-05,-6.21649E-05,-1.13073E-04,-1.20560E-04,   F33310
     *-5.92397E-05, 5.24432E-05, 9.41628E-05,-3.47458E-07, 5.33267E-05/   F33320
      DATA C11921 /                                                       F33330
     * 8.92961E-05, 2.75694E-05,-7.48460E-06,-2.15504E-05, 1.05501E-06,   F33340
     * 6.30910E-06, 5.94620E-07,-2.45194E-05,-1.59657E-05, 7.93610E-07,   F33350
     *-1.05319E-05,-2.36584E-05,-3.95700E-05,-6.57225E-05,-5.23797E-05,   F33360
     *-1.82588E-05,-1.43240E-05,-3.29989E-05,-6.48909E-05,-2.41326E-05,   F33370
     *-1.89195E-05,-4.64607E-05,-1.00739E-05,-1.35033E-05,-6.49945E-05,   F33380
     *-5.19986E-05,-6.68505E-05,-1.31530E-04,-1.45464E-04,-1.46815E-04,   F33390
     *-1.39684E-04,-1.23205E-04,-1.26738E-04,-1.93822E-04,-2.37508E-04,   F33400
     *-2.52917E-04,-1.91110E-04,-1.36217E-04,-9.41093E-05,-1.20601E-04,   F33410
     *-1.17295E-04,-9.57420E-05,-1.57227E-04,-1.62795E-04,-1.12201E-04,   F33420
     *-1.20419E-04,-1.10597E-04,-7.61223E-05,-6.27167E-05,-5.54733E-05,   F33430
     *-5.50437E-05,-5.14148E-05,-3.59591E-05, 1.09906E-05, 5.94396E-06,   F33440
     *-1.38597E-05,-8.80857E-06,-3.13101E-05,-6.31715E-05,-4.04264E-05,   F33450
     *-1.66405E-05, 7.94396E-06,-3.41772E-06,-4.03175E-05,-1.06888E-04,   F33460
     *-9.50526E-05,-7.46111E-05,-5.09617E-05,-6.70981E-05,-7.93529E-05,   F33470
     *-5.58423E-05,-1.01523E-04,-1.62269E-04,-1.69958E-04,-1.37786E-04,   F33480
     *-8.79862E-05,-1.46838E-04,-1.66938E-04,-1.51380E-04,-1.62184E-04/   F33490
      DATA C12001 /                                                       F33500
     *-1.61105E-04,-1.42088E-04,-1.57033E-04,-1.65294E-04,-1.45079E-04,   F33510
     *-9.76982E-05,-6.09891E-05,-1.01719E-04,-1.03049E-04,-8.85546E-05,   F33520
     *-1.47754E-04,-1.44542E-04,-8.34620E-05,-8.99440E-05,-7.11901E-05,   F33530
     *-1.57480E-05,-8.81797E-05,-1.56314E-04,-1.65952E-04,-1.80986E-04,   F33540
     *-2.04610E-04,-2.58669E-04,-2.16016E-04,-1.21582E-04,-1.44929E-04,   F33550
     *-1.72886E-04,-2.05950E-04,-1.93829E-04,-1.67518E-04,-1.22969E-04,   F33560
     *-1.13060E-04,-1.14854E-04,-1.26198E-04,-1.24288E-04,-1.19519E-04,   F33570
     *-1.50456E-04,-1.53286E-04,-1.32231E-04,-7.42672E-05,-2.23129E-05,   F33580
     * 1.79115E-05, 1.42073E-05,-1.21676E-05,-7.56567E-05,-1.03423E-04,   F33590
     *-1.10373E-04,-8.77244E-05,-6.43485E-05,-4.05156E-05,-6.24405E-05,   F33600
     *-5.70375E-05,-2.36695E-06,-3.75929E-05,-7.97119E-05,-6.70419E-05,   F33610
     *-6.99475E-05,-8.19748E-05,-1.06895E-04,-1.31422E-04,-1.55438E-04,   F33620
     *-1.61937E-04,-1.62626E-04,-1.54977E-04,-1.77814E-04,-2.00386E-04,   F33630
     *-1.87407E-04,-2.07243E-04,-2.44672E-04,-2.19014E-04,-2.13695E-04,   F33640
     *-2.32440E-04,-1.85194E-04,-1.51172E-04,-1.69834E-04,-1.73780E-04,   F33650
     *-1.75232E-04,-2.00698E-04,-1.82826E-04,-1.27786E-04,-1.33633E-04/   F33660
      DATA C12081 /                                                       F33670
     *-1.21317E-04,-7.50390E-05,-1.06743E-04,-1.40805E-04,-1.06336E-04,   F33680
     *-9.46654E-05,-9.78182E-05,-1.19906E-04,-1.14160E-04,-7.28186E-05,   F33690
     *-1.07652E-04,-1.20978E-04,-3.79658E-05,-3.16113E-05,-6.02417E-05,   F33700
     *-7.51148E-05,-5.56145E-05,-6.77421E-06,-1.74321E-05,-4.67952E-05,   F33710
     *-1.05000E-04,-6.29932E-05,-4.74356E-06,-2.83397E-05,-4.65192E-05,   F33720
     *-6.04574E-05,-4.33970E-05,-3.18311E-05,-3.02321E-05,-4.49667E-05,   F33730
     *-6.85347E-05,-1.11375E-04,-1.16293E-04,-9.38757E-05,-1.38594E-04,   F33740
     *-1.60483E-04,-1.48344E-04,-1.33436E-04,-1.27387E-04,-1.59508E-04,   F33750
     *-1.74026E-04,-1.72170E-04,-1.49196E-04,-1.33233E-04,-1.22382E-04,   F33760
     *-1.78156E-04,-2.21349E-04,-2.41846E-04,-2.06549E-04,-1.68283E-04,   F33770
     *-1.89512E-04,-1.44523E-04,-4.67953E-05,-1.00334E-04,-1.23478E-04,   F33780
     *-8.14024E-05,-9.18016E-05,-1.17536E-04,-1.36160E-04,-1.38780E-04,   F33790
     *-1.27749E-04,-1.45598E-04,-1.55964E-04,-1.45120E-04,-1.25544E-04,   F33800
     *-1.05692E-04,-1.17639E-04,-1.24142E-04,-1.24749E-04,-1.63878E-04,   F33810
     *-1.97021E-04,-1.98617E-04,-2.69136E-04,-3.68357E-04,-2.33702E-04,   F33820
     *-1.61830E-04,-1.78578E-04,-2.01839E-04,-2.28731E-04,-2.63606E-04/   F33830
      DATA C12161 /                                                       F33840
     *-2.44698E-04,-1.86451E-04,-2.20546E-04,-2.22752E-04,-1.55169E-04,   F33850
     *-1.25100E-04,-1.09794E-04,-9.59016E-05,-1.03857E-04,-1.35573E-04,   F33860
     *-1.73780E-04,-1.82457E-04,-9.39821E-05,-1.18245E-04,-2.11563E-04,   F33870
     *-1.37392E-04,-9.28173E-05,-9.71073E-05,-9.72535E-05,-9.39557E-05,   F33880
     *-7.50117E-05,-6.70754E-05,-7.01186E-05,-5.76151E-05,-5.18785E-05,   F33890
     *-7.14209E-05,-7.01682E-05,-5.61614E-05,-8.92769E-05,-1.06238E-04,   F33900
     *-9.70294E-05,-6.70229E-05,-4.69214E-05,-1.53105E-04,-2.02326E-04,   F33910
     *-1.90395E-04,-2.04367E-04,-2.16787E-04,-2.08725E-04,-1.78119E-04,   F33920
     *-1.31043E-04,-1.32204E-04,-1.51522E-04,-2.05143E-04,-1.77144E-04,   F33930
     *-1.16130E-04,-1.44440E-04,-1.66010E-04,-1.78206E-04,-1.61163E-04,   F33940
     *-1.46351E-04,-1.96722E-04,-2.27027E-04,-2.37243E-04,-2.25235E-04,   F33950
     *-1.99552E-04,-1.40238E-04,-1.26311E-04,-1.42746E-04,-1.19028E-04,   F33960
     *-1.18750E-04,-1.72076E-04,-1.72120E-04,-1.48285E-04,-1.85116E-04,   F33970
     *-1.98602E-04,-1.74016E-04,-1.37913E-04,-1.01221E-04,-9.69581E-05,   F33980
     *-1.08794E-04,-1.39433E-04,-1.38575E-04,-1.32088E-04,-1.37431E-04,   F33990
     *-1.30033E-04,-1.10829E-04,-1.35604E-04,-1.66515E-04,-1.98167E-04/   F34000
      DATA C12241 /                                                       F34010
     *-1.97716E-04,-1.74019E-04,-1.64719E-04,-1.64779E-04,-1.85725E-04,   F34020
     *-2.28526E-04,-2.84329E-04,-1.82449E-04,-1.30747E-04,-1.93620E-04,   F34030
     *-2.28529E-04,-2.47361E-04,-1.90001E-04,-1.66278E-04,-2.02540E-04,   F34040
     *-2.31811E-04,-2.53772E-04,-2.08629E-04,-1.85021E-04,-1.93989E-04,   F34050
     *-2.16568E-04,-2.38288E-04,-1.94453E-04,-1.87154E-04,-2.30493E-04,   F34060
     *-2.34696E-04,-2.30351E-04,-2.60562E-04,-2.86427E-04,-3.06699E-04,   F34070
     *-2.79131E-04,-2.49392E-04,-3.03389E-04,-3.10346E-04,-2.61782E-04,   F34080
     *-2.30905E-04,-2.11669E-04,-2.37680E-04,-2.38194E-04,-2.10955E-04/   F34090
C                                                                         F34100
      END                                                                 F34110
C
C     --------------------------------------------------------------
C
      SUBROUTINE O3HHT2 (V1C,V2C,DVC,NPTC,C)                              F34120
C                                                                         F34130
      IMPLICIT REAL*8           (V)                                     ! F34140
C                                                                         F34150
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F34160
      COMMON /O3HH2/ V1S,V2S,DVS,NPTS,S(2687)                             F34170
      DIMENSION C(*)                                                      F34180
C                                                                         F34190
      DVC = DVS                                                           F34200
      V1C = V1ABS-DVC                                                     F34210
      V2C = V2ABS+DVC                                                     F34220
C                                                                         F34230
      I1 = (V1C-V1S)/DVS                                                  F34240
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F34310
         I = I1+J                                                         F34320
         C(J) = 0.                                                        F34330
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F34340
         C(J) = S(I)                                                      F34350
   10 CONTINUE                                                            F34360
C                                                                         F34370
      RETURN                                                              F34380
C                                                                         F34390
      END                                                                 F34400
C
C     --------------------------------------------------------------
C
      BLOCK DATA BO3HH2                                                   F34410
C                                                                         F34420
      IMPLICIT REAL*8           (V)                                     ! F34430
C                                                                         F34440
C     RATIO (C2/C0)                                                       F34450
C     DATA FROM BASS 1985                                                 F34460
C                                                                         F34470
C     NOW INCLUDES MOLINA & MOLINA AT 273K WITH THE TEMPERATURE           F34480
C     DEPENDENCE DETERMINED FROM THE 195K HARVARD MEASUREMENTS,           F34490
C     EMPLOYING THE BASS ALGORITHM
C
C              (CO(1+C1*(T-273.15)+C2*(T-273.15)**2);                     F34500
C
C     THIS IS ONLY FOR THE WAVELENGTH RANGE FROM .34 TO .35 MICRONS;      F34510
C     OTHERWISE, THE BASS DATA ALONE HAVE BEEN EMPLOYED BETWEEN           F34520
C     .34 AND .245 MICRONS.                                               F34530
C                                                                         F34540
C     NEW T-DEPENDENT X-SECTIONS BETWEEN .345 AND .36 MICRONS             F34550
C     HAVE NOW BEEN ADDED, BASED ON WORK BY CACCIANI, DISARRA             F34560
C     AND FIOCCO, UNIVERSITY OF ROME, 1987.  QUADRATIC TEMP               F34570
C     HAS BEEN DERIVED, AS ABOVE.                                         F34580
C                                                                         F34590
C     AGREEMENT AMONGST THE FOUR DATA SETS IS REASONABLE (<10%)           F34600
C     AND OFTEN EXCELLENT (0-3%)                                          F34610
C                                                                         F34620
C                                                                         F34630
      COMMON /O3HH2/ V1C,V2C,DVC,NC,                                      F34640
     *               O32001(85),C20086(80),C20166(80),C20246(65),         F34650
     *               C20311(16),C20327(80),C20407( 1),                    F34660
     *               C20001(80),C20081(80),C20161(80),C20241(80),         F34670
     *               C20321(80),C20401(80),C20481(80),C20561(80),         F34680
     *               C20641(80),C20721(80),C20801(80),C20881(80),         F34690
     *               C20961(80),C21041(80),C21121(80),C21201(80),         F34700
     *               C21281(80),C21361(80),C21441(80),C21521(80),         F34710
     *               C21601(80),C21681(80),C21761(80),C21841(80),         F34720
     *               C21921(80),C22001(80),C22081(80),C22161(80),         F34730
     *               C22241(40)                                           F34740
C                                                                         F34750
C     DATA V1C /29405./, V2C /40800./ ,DVC /5./, NC /2280/   BASS         F34760
C                                                                         F34770
      DATA V1C /27370./, V2C /40800./ ,DVC /5./, NC /2687/                F34780
C                                                                         F34790
      DATA O32001/85*1.0E-5/                                              F34800
C                                                                         F34810
      DATA C20086/                                                        F34820
     * 1.29359E-05, 1.55806E-05, 2.00719E-05, 2.64912E-05, 3.48207E-05,   F34830
     * 4.36986E-05, 5.31318E-05, 6.13173E-05, 6.89465E-05, 7.56793E-05,   F34840
     * 8.26345E-05, 8.90916E-05, 9.38759E-05, 9.22998E-05, 9.03184E-05,   F34850
     * 8.65369E-05, 8.58531E-05, 8.55635E-05, 8.40418E-05, 8.11983E-05,   F34860
     * 7.58246E-05, 7.29282E-05, 7.32629E-05, 7.04060E-05, 6.71451E-05,   F34870
     * 6.56515E-05, 6.68943E-05, 6.32785E-05, 5.88386E-05, 5.70860E-05,   F34880
     * 5.64435E-05, 5.49441E-05, 5.70845E-05, 5.89357E-05, 6.14433E-05,   F34890
     * 5.91790E-05, 5.31727E-05, 5.14007E-05, 4.74318E-05, 4.35356E-05,   F34900
     * 3.93903E-05, 3.70963E-05, 3.63867E-05, 4.05296E-05, 4.48891E-05,   F34910
     * 5.37190E-05, 5.70440E-05, 4.60408E-05, 5.25778E-05, 6.81728E-05,   F34920
     * 7.27275E-05, 6.81353E-05, 6.48386E-05, 5.46521E-05, 4.93098E-05,   F34930
     * 5.04438E-05, 5.30309E-05, 5.28788E-05, 5.47387E-05, 4.52523E-05,   F34940
     * 5.29451E-05, 7.42215E-05, 1.08971E-04, 1.40085E-04, 1.46553E-04,   F34950
     * 1.43526E-04, 1.39051E-04, 1.40983E-04, 1.45564E-04, 1.55589E-04,   F34960
     * 1.66142E-04, 1.82840E-04, 2.06486E-04, 2.24339E-04, 2.29268E-04,   F34970
     * 2.13109E-04, 2.00305E-04, 1.99955E-04, 2.18566E-04, 2.24182E-04/   F34980
      DATA C20166/                                                        F34990
     * 2.33505E-04, 2.31824E-04, 2.22666E-04, 2.23905E-04, 2.38131E-04,   F35000
     * 2.54322E-04, 2.69548E-04, 2.62953E-04, 2.67609E-04, 2.70567E-04,   F35010
     * 2.70689E-04, 2.68251E-04, 2.66029E-04, 2.60053E-04, 2.61689E-04,   F35020
     * 2.56582E-04, 2.43655E-04, 2.38792E-04, 2.45309E-04, 2.31061E-04,   F35030
     * 2.22837E-04, 2.16440E-04, 2.19032E-04, 1.85634E-04, 1.74638E-04,   F35040
     * 1.51767E-04, 1.38480E-04, 1.32506E-04, 1.28317E-04, 1.26855E-04,   F35050
     * 1.27123E-04, 1.24040E-04, 1.19202E-04, 1.28649E-04, 1.36271E-04,   F35060
     * 1.42080E-04, 1.47804E-04, 1.39534E-04, 1.27284E-04, 1.09554E-04,   F35070
     * 8.69470E-05, 6.72096E-05, 5.23407E-05, 5.12433E-05, 5.15794E-05,   F35080
     * 4.94683E-05, 4.95809E-05, 4.07499E-05, 3.14984E-05, 1.46457E-05,   F35090
     * 6.98660E-06, 1.85313E-05, 5.48879E-05, 1.09447E-04, 1.52536E-04,   F35100
     * 1.78778E-04, 1.91128E-04, 1.99161E-04, 2.02937E-04, 1.95527E-04,   F35110
     * 1.92214E-04, 1.83376E-04, 1.81710E-04, 1.82283E-04, 1.75182E-04,   F35120
     * 1.72406E-04, 1.68170E-04, 1.67400E-04, 1.69469E-04, 1.69092E-04,   F35130
     * 1.65985E-04, 1.66912E-04, 1.74226E-04, 1.85036E-04, 1.85517E-04,   F35140
     * 1.85805E-04, 1.73809E-04, 1.67628E-04, 1.57690E-04, 1.54952E-04/   F35150
      DATA C20246/                                                        F35160
     * 1.53707E-04, 1.57710E-04, 1.58175E-04, 1.67253E-04, 1.82079E-04,   F35170
     * 1.91285E-04, 1.96564E-04, 2.03822E-04, 1.93736E-04, 1.82924E-04,   F35180
     * 1.73610E-04, 1.69904E-04, 1.66725E-04, 1.63747E-04, 1.63129E-04,   F35190
     * 1.62435E-04, 1.67218E-04, 1.69507E-04, 1.70744E-04, 1.65839E-04,   F35200
     * 1.72077E-04, 1.67734E-04, 1.51487E-04, 1.43770E-04, 1.37435E-04,   F35210
     * 1.25172E-04, 1.12395E-04, 1.07991E-04, 1.00345E-04, 9.36611E-05,   F35220
     * 9.59763E-05, 9.26600E-05, 1.00120E-04, 1.04746E-04, 1.10222E-04,   F35230
     * 1.03308E-04, 8.97457E-05, 7.91634E-05, 7.50275E-05, 8.30832E-05,   F35240
     * 1.01191E-04, 1.21560E-04, 1.34840E-04, 1.38712E-04, 1.41746E-04,   F35250
     * 1.39578E-04, 1.37052E-04, 1.33850E-04, 1.26641E-04, 1.21342E-04,   F35260
     * 1.17669E-04, 1.25973E-04, 1.33623E-04, 1.33839E-04, 1.24427E-04,   F35270
     * 1.02462E-04, 8.76101E-05, 8.27912E-05, 8.29040E-05, 7.78590E-05,   F35280
     * 7.39042E-05, 6.45765E-05, 5.70151E-05, 5.11846E-05, 4.83163E-05/   F35290
      DATA C20311/                                                        F35300
     *                                                     5.4470E-05,    F35310
     * 5.3312E-05,  5.3135E-05,  5.3619E-05,  5.3686E-05,  5.2308E-05,    F35320
     * 5.0441E-05,  4.8402E-05,  4.7476E-05,  4.6215E-05,  4.4507E-05,    F35330
     * 4.3830E-05,  4.0508E-05,  3.8931E-05,  3.5525E-05,  3.4722E-05/    F35340
      DATA C20327/                                                        F35350
     * 3.2743E-05,  2.8456E-05,  2.8318E-05,  2.8132E-05,  2.6221E-05,    F35360
     * 2.5673E-05,  2.5521E-05,  2.4588E-05,  2.4093E-05,  2.2787E-05,    F35370
     * 2.1241E-05,  1.8553E-05,  1.5871E-05,  1.3462E-05,  1.2553E-05,    F35380
     * 1.6276E-05,  2.8296E-05,  3.8817E-05,  4.2733E-05,  4.2429E-05,    F35390
     * 4.0954E-05,  3.9868E-05,  3.7669E-05,  3.6312E-05,  3.5535E-05,    F35400
     * 3.5895E-05,  3.6349E-05,  3.9033E-05,  4.4512E-05,  5.0066E-05,    F35410
     * 5.4572E-05,  5.6710E-05,  5.6615E-05,  5.7520E-05,  5.8034E-05,    F35420
     * 5.7927E-05,  5.6027E-05,  5.5242E-05,  5.4974E-05,  5.2927E-05,    F35430
     * 5.1638E-05,  5.2027E-05,  5.1420E-05,  5.1618E-05,  5.0253E-05,    F35440
     * 5.0509E-05,  4.9376E-05,  5.0135E-05,  4.9191E-05,  4.9210E-05,    F35450
     * 4.8216E-05,  4.7487E-05,  4.5749E-05,  4.5884E-05,  4.3852E-05,    F35460
     * 4.3824E-05,  4.2612E-05,  4.0349E-05,  4.0177E-05,  3.7474E-05,    F35470
     * 3.8120E-05,  3.6915E-05,  3.5823E-05,  3.5186E-05,  3.3638E-05,    F35480
     * 3.3451E-05,  3.2428E-05,  3.2349E-05,  3.0183E-05,  2.8436E-05,    F35490
     * 2.6440E-05,  2.3597E-05,  2.1875E-05,  1.8164E-05,  1.6430E-05,    F35500
     * 1.3159E-05,  9.2907E-06,  7.4243E-06,  6.0469E-06,  5.4951E-06/    F35510
      DATA C20407/                                                        F35520
     * 8.7642E-06/                                                        F35530
      DATA C20001 /                                                       F35540
     * 2.16295E-05, 1.69111E-05, 5.39633E-05, 1.01866E-04, 8.28657E-05,   F35550
     * 9.16593E-05, 8.88666E-05, 1.37764E-04, 1.44322E-04, 1.20659E-04,   F35560
     * 1.10332E-04, 1.01317E-04, 9.09964E-05, 1.17148E-04, 1.18000E-04,   F35570
     * 7.21801E-05, 1.10550E-04, 1.32672E-04, 1.02474E-04, 1.10434E-04,   F35580
     * 1.38759E-04, 8.92135E-05, 9.18239E-05, 9.08256E-05, 7.02969E-05,   F35590
     * 1.12827E-04, 8.25561E-05, 1.39555E-04, 6.72239E-05, 7.82804E-05,   F35600
     * 8.56258E-05, 8.61068E-05, 7.16732E-05, 6.25720E-05, 5.23957E-05,   F35610
     * 3.78801E-05, 4.37281E-05, 4.99821E-05, 5.96976E-05, 7.19070E-05,   F35620
     * 3.89579E-05, 5.30171E-05, 3.92507E-05, 4.93901E-05, 4.53047E-05,   F35630
     * 4.89955E-05, 4.61649E-05, 3.75742E-05, 3.14124E-05, 2.37893E-05,   F35640
     * 3.34899E-06, 3.08080E-05, 5.35883E-05, 3.39838E-05, 7.02334E-05,   F35650
     * 7.24784E-05, 7.46533E-05, 6.22257E-05, 6.38945E-05, 6.73423E-05,   F35660
     * 4.51321E-05, 5.91854E-05, 5.51601E-05, 4.41923E-05, 3.59217E-05,   F35670
     * 4.08520E-05, 6.15981E-05, 6.66549E-05, 8.26031E-05, 1.13556E-04,   F35680
     * 8.72988E-05, 9.71052E-05, 9.31839E-05, 8.73745E-05, 8.61717E-05,   F35690
     * 6.05645E-05, 6.51131E-05, 6.93393E-05, 7.01096E-05, 6.43565E-05/   F35700
      DATA C20081 /                                                       F35710
     * 7.36929E-05, 7.66881E-05, 7.60815E-05, 7.13570E-05, 8.40487E-05,   F35720
     * 8.51489E-05, 7.54168E-05, 6.72694E-05, 4.75508E-05, 3.59379E-05,   F35730
     * 4.24698E-05, 4.17850E-05, 4.56047E-05, 4.12779E-05, 4.55933E-05,   F35740
     * 4.27941E-05, 4.42230E-05, 3.68525E-05, 3.83392E-05, 3.83722E-05,   F35750
     * 4.64904E-05, 3.33878E-05, 3.53027E-05, 3.54694E-05, 2.36233E-05,   F35760
     * 2.99641E-05, 2.56097E-05, 2.14134E-05, 2.74403E-05, 2.83896E-05,   F35770
     * 3.17082E-05, 1.75526E-05, 2.80382E-05, 3.18009E-05, 4.08715E-05,   F35780
     * 4.77807E-05, 5.00609E-05, 5.12459E-05, 4.44062E-05, 4.74942E-05,   F35790
     * 4.99882E-05, 5.18837E-05, 5.03246E-05, 5.55168E-05, 5.35853E-05,   F35800
     * 4.81834E-05, 6.66231E-05, 5.26670E-05, 6.84700E-05, 6.53412E-05,   F35810
     * 5.71740E-05, 4.61076E-05, 3.90239E-05, 4.72924E-05, 6.32194E-05,   F35820
     * 5.20868E-05, 4.81039E-05, 3.71748E-05, 4.37492E-05, 3.63959E-05,   F35830
     * 3.79823E-05, 3.72225E-05, 3.02360E-05, 3.22961E-05, 3.43398E-05,   F35840
     * 3.57176E-05, 2.65446E-05, 3.29388E-05, 1.65455E-05, 2.66173E-05,   F35850
     * 1.74277E-05, 1.74324E-05, 1.27879E-05, 1.46247E-05, 1.92378E-05,   F35860
     * 2.20049E-05, 1.44790E-05, 2.49244E-05, 2.29209E-05, 1.76192E-05/   F35870
      DATA C20161 /                                                       F35880
     * 1.84528E-05, 2.54350E-05, 3.33972E-05, 3.69190E-05, 2.92139E-05,   F35890
     * 2.47666E-05, 2.86764E-05, 1.48163E-05, 1.80461E-05, 2.84545E-05,   F35900
     * 2.41064E-05, 2.85721E-05, 3.31996E-05, 3.75973E-05, 3.73874E-05,   F35910
     * 4.69293E-05, 5.12665E-05, 5.35607E-05, 4.64577E-05, 3.59887E-05,   F35920
     * 3.39168E-05, 3.89746E-05, 3.12196E-05, 3.70907E-05, 3.95172E-05,   F35930
     * 4.61642E-05, 4.26029E-05, 4.17856E-05, 4.51437E-05, 4.04189E-05,   F35940
     * 4.19251E-05, 4.53977E-05, 3.69860E-05, 4.20904E-05, 3.69735E-05,   F35950
     * 3.57898E-05, 3.47729E-05, 3.14280E-05, 2.71197E-05, 3.34380E-05,   F35960
     * 2.69843E-05, 2.88036E-05, 2.51912E-05, 2.45699E-05, 2.23184E-05,   F35970
     * 2.50563E-05, 2.24493E-05, 1.77101E-05, 1.64763E-05, 1.34978E-05,   F35980
     * 1.57081E-05, 1.45966E-05, 1.02722E-05, 2.07177E-05, 1.47662E-05,   F35990
     * 1.50721E-05, 1.24431E-05, 1.51572E-05, 1.92210E-05, 2.06047E-05,   F36000
     * 2.02921E-05, 3.22062E-05, 2.37112E-05, 1.94803E-05, 2.40726E-05,   F36010
     * 2.11531E-05, 1.89158E-05, 2.46957E-05, 2.63175E-05, 2.57747E-05,   F36020
     * 2.22047E-05, 2.52755E-05, 2.80848E-05, 3.75157E-05, 4.09915E-05,   F36030
     * 4.04853E-05, 3.21661E-05, 3.15652E-05, 3.21576E-05, 3.67060E-05/   F36040
      DATA C20241 /                                                       F36050
     * 3.13071E-05, 2.84939E-05, 2.71169E-05, 2.99559E-05, 2.94631E-05,   F36060
     * 3.26716E-05, 2.99028E-05, 2.60045E-05, 3.15375E-05, 3.12895E-05,   F36070
     * 2.77767E-05, 2.43976E-05, 2.10764E-05, 2.22725E-05, 2.04581E-05,   F36080
     * 1.63509E-05, 1.60028E-05, 1.60294E-05, 1.62366E-05, 1.89293E-05,   F36090
     * 1.79675E-05, 1.89259E-05, 1.68300E-05, 1.99460E-05, 2.42370E-05,   F36100
     * 2.64738E-05, 1.93137E-05, 1.39460E-05, 1.32222E-05, 1.38752E-05,   F36110
     * 1.62071E-05, 1.79652E-05, 1.63772E-05, 1.56251E-05, 1.81918E-05,   F36120
     * 1.46111E-05, 2.92174E-05, 2.94263E-05, 2.46180E-05, 2.93333E-05,   F36130
     * 3.13657E-05, 2.97686E-05, 2.78387E-05, 2.40924E-05, 2.93369E-05,   F36140
     * 2.93747E-05, 2.77665E-05, 3.00814E-05, 3.01068E-05, 3.62275E-05,   F36150
     * 3.56613E-05, 3.66913E-05, 3.56280E-05, 3.52856E-05, 3.63928E-05,   F36160
     * 2.96738E-05, 2.90314E-05, 2.62972E-05, 2.15250E-05, 1.97910E-05,   F36170
     * 2.02314E-05, 2.20209E-05, 2.05131E-05, 2.12017E-05, 1.96689E-05,   F36180
     * 1.61907E-05, 1.57662E-05, 1.58239E-05, 1.54650E-05, 1.46376E-05,   F36190
     * 1.32891E-05, 1.30511E-05, 1.17635E-05, 1.28585E-05, 1.12887E-05,   F36200
     * 1.32627E-05, 1.31833E-05, 1.68679E-05, 1.98092E-05, 2.70744E-05/   F36210
      DATA C20321 /                                                       F36220
     * 2.22033E-05, 1.63430E-05, 1.61104E-05, 1.50865E-05, 1.54382E-05,   F36230
     * 1.55654E-05, 1.67924E-05, 1.89185E-05, 1.96791E-05, 2.14894E-05,   F36240
     * 2.76137E-05, 2.67339E-05, 2.79423E-05, 2.54664E-05, 3.10707E-05,   F36250
     * 2.72745E-05, 2.60940E-05, 2.47736E-05, 2.21105E-05, 2.20357E-05,   F36260
     * 2.26499E-05, 2.34137E-05, 2.29537E-05, 2.36157E-05, 2.48244E-05,   F36270
     * 2.26667E-05, 2.07781E-05, 2.11702E-05, 1.91214E-05, 1.62172E-05,   F36280
     * 1.61285E-05, 1.63952E-05, 1.68156E-05, 1.61236E-05, 1.56611E-05,   F36290
     * 1.47697E-05, 1.50856E-05, 1.44169E-05, 1.63816E-05, 1.74283E-05,   F36300
     * 1.49853E-05, 1.62444E-05, 1.70007E-05, 1.60371E-05, 1.22713E-05,   F36310
     * 1.45518E-05, 1.35051E-05, 1.40787E-05,-1.54925E-05,-2.15204E-05,   F36320
     *-4.04516E-06, 2.22439E-05, 3.21262E-05, 3.83792E-05, 4.44462E-05,   F36330
     * 4.44192E-05, 2.77328E-05, 4.10549E-06, 4.48758E-06,-1.27771E-05,   F36340
     *-2.17204E-05,-1.23979E-05,-1.04928E-05, 7.43085E-06, 1.55350E-05,   F36350
     * 3.15204E-05, 3.17601E-05, 2.93677E-05, 3.42485E-05, 3.87087E-05,   F36360
     * 3.61242E-05, 2.62406E-05, 3.31686E-05, 3.54314E-05, 2.50625E-05,   F36370
     * 2.60444E-05, 4.10729E-05, 3.47247E-05, 3.31716E-05, 3.34778E-05/   F36380
      DATA C20401 /                                                       F36390
     * 4.03029E-05, 4.09241E-05, 3.96717E-05, 3.53410E-05, 2.81048E-05,   F36400
     * 1.98891E-05, 1.92314E-05, 2.82525E-05, 3.76641E-05, 4.34135E-05,   F36410
     * 4.24570E-05, 3.98429E-05, 3.29417E-05, 2.16679E-05, 8.88085E-06,   F36420
     *-5.05319E-06,-8.14815E-06,-5.01930E-06, 7.13565E-06, 2.00949E-05,   F36430
     * 2.65988E-05, 2.77656E-05, 2.09299E-05, 1.98968E-05, 2.04835E-05,   F36440
     * 1.75254E-05, 6.48674E-06, 3.14323E-06, 1.93242E-06, 3.86745E-06,   F36450
     * 1.39727E-05, 2.10731E-05, 2.66432E-05, 2.69551E-05, 2.57453E-05,   F36460
     * 2.72834E-05, 2.58860E-05, 2.51266E-05, 1.76048E-05, 2.03072E-05,   F36470
     * 2.61960E-05, 2.36230E-05, 1.81172E-05, 1.33972E-05, 1.60959E-05,   F36480
     * 1.61081E-05, 2.34099E-05, 2.64979E-05, 2.36894E-05, 2.13665E-05,   F36490
     * 2.16774E-05, 2.52566E-05, 1.99785E-05, 1.40414E-05, 1.39948E-05,   F36500
     * 1.32637E-05, 7.24742E-06, 1.11395E-06,-1.27323E-06, 4.56637E-07,   F36510
     * 6.93250E-06, 5.07198E-06, 7.90632E-06, 9.08149E-06, 1.03602E-05,   F36520
     * 2.17425E-05, 2.71741E-05, 2.16875E-05, 1.95088E-05, 1.56568E-05,   F36530
     * 8.41152E-06, 1.26749E-05, 1.17673E-05, 9.96037E-06, 1.21982E-05,   F36540
     * 1.31854E-05, 1.50216E-05, 1.72214E-05, 2.02773E-05, 2.09625E-05/   F36550
      DATA C20481 /                                                       F36560
     * 1.66656E-05, 1.45666E-05, 1.66608E-05, 2.04989E-05, 2.21395E-05,   F36570
     * 2.35993E-05, 2.69390E-05, 2.13921E-05, 1.72643E-05, 1.70995E-05,   F36580
     * 1.78241E-05, 1.85308E-05, 1.80360E-05, 1.48619E-05, 1.90369E-05,   F36590
     * 1.51089E-05, 1.22705E-05, 1.62608E-05, 1.41637E-05, 1.23786E-05,   F36600
     * 7.02677E-06, 8.89811E-06, 1.07379E-05, 1.23677E-05, 1.48196E-05,   F36610
     * 2.05770E-05, 1.70994E-05, 1.00072E-05, 1.76119E-05, 1.41779E-05,   F36620
     * 1.34358E-05, 1.58674E-05, 1.65837E-05, 1.69569E-05, 1.40381E-05,   F36630
     * 1.46118E-05, 1.30556E-05, 1.97204E-05, 1.97488E-05, 1.64524E-05,   F36640
     * 1.73764E-05, 1.66355E-05, 1.64419E-05, 1.65486E-05, 1.21523E-05,   F36650
     * 1.51513E-05, 1.60354E-05, 1.38528E-05, 1.45538E-05, 1.71702E-05,   F36660
     * 1.56336E-05, 1.31279E-05, 1.47346E-05, 1.70719E-05, 1.75588E-05,   F36670
     * 1.55187E-05, 1.29598E-05, 1.38463E-05, 1.35382E-05, 1.16062E-05,   F36680
     * 1.37014E-05, 1.34487E-05, 1.15536E-05, 1.33597E-05, 9.24478E-06,   F36690
     * 7.28477E-06, 1.40321E-05, 1.31518E-05, 1.03118E-05, 8.59764E-06,   F36700
     * 1.57138E-05, 1.20792E-05, 1.49440E-05, 1.34375E-05, 1.54686E-05,   F36710
     * 1.65346E-05, 1.33823E-05, 1.37238E-05, 1.36128E-05, 1.46206E-05/   F36720
      DATA C20561 /                                                       F36730
     * 1.40777E-05, 1.59980E-05, 1.30180E-05, 1.01390E-05, 1.12366E-05,   F36740
     * 9.86099E-06, 1.10702E-05, 1.26783E-05, 9.51072E-06, 8.07299E-06,   F36750
     * 1.22955E-05, 1.53506E-05, 1.29711E-05, 9.78759E-06, 1.28800E-05,   F36760
     * 1.39702E-05, 1.64832E-05, 1.06473E-05, 1.15419E-05, 1.63795E-05,   F36770
     * 1.69837E-05, 1.72726E-05, 1.77231E-05, 1.62337E-05, 1.20881E-05,   F36780
     * 1.13210E-05, 1.20531E-05, 1.31374E-05, 1.22259E-05, 1.27802E-05,   F36790
     * 1.38962E-05, 8.87355E-06, 9.42264E-06, 1.02075E-05, 7.91816E-06,   F36800
     * 9.66835E-06, 1.24921E-05, 8.43227E-06, 1.10637E-05, 1.03958E-05,   F36810
     * 9.40996E-06, 1.22922E-05, 1.21088E-05, 1.30116E-05, 1.18776E-05,   F36820
     * 1.42245E-05, 1.34745E-05, 1.11165E-05, 1.29914E-05, 1.29801E-05,   F36830
     * 1.10895E-05, 1.12331E-05, 9.03490E-06, 9.33726E-06, 9.63923E-06,   F36840
     * 1.11299E-05, 9.53481E-06, 1.21708E-05, 1.11951E-05, 7.22558E-06,   F36850
     * 6.66928E-06, 1.08926E-05, 1.07870E-05, 9.23485E-06, 8.50452E-06,   F36860
     * 9.41914E-06, 8.74027E-06, 8.93322E-06, 9.79061E-06, 8.26490E-06,   F36870
     * 8.37630E-06, 1.17064E-05, 1.10176E-05, 1.11587E-05, 9.45563E-06,   F36880
     * 1.18352E-05, 7.79327E-06, 9.22766E-06, 1.01868E-05, 8.23925E-06/   F36890
      DATA C20641 /                                                       F36900
     * 9.23706E-06, 1.04428E-05, 8.80392E-06, 9.37098E-06, 7.43126E-06,   F36910
     * 7.01424E-06, 9.29360E-06, 8.97171E-06, 9.31718E-06, 9.87118E-06,   F36920
     * 8.11419E-06, 8.77416E-06, 9.96927E-06, 8.87533E-06, 9.33163E-06,   F36930
     * 7.41505E-06, 9.39988E-06, 1.17932E-05, 1.03287E-05, 9.17415E-06,   F36940
     * 8.43035E-06, 8.00040E-06, 8.33346E-06, 7.66787E-06, 7.18411E-06,   F36950
     * 1.06236E-05, 1.05559E-05, 8.49187E-06, 9.22472E-06, 8.16512E-06,   F36960
     * 8.35687E-06, 1.06325E-05, 9.80273E-06, 9.01599E-06, 9.20499E-06,   F36970
     * 9.98417E-06, 9.23191E-06, 6.98769E-06, 5.17748E-06, 4.57130E-06,   F36980
     * 8.18492E-06, 9.98095E-06, 7.52148E-06, 1.33038E-05, 8.17630E-06,   F36990
     * 1.02454E-05, 9.62706E-06, 9.44304E-06, 8.86704E-06, 8.88116E-06,   F37000
     * 8.79062E-06, 8.20042E-06, 8.55789E-06, 9.26249E-06, 1.00467E-05,   F37010
     * 7.96012E-06, 9.08773E-06, 1.01481E-05, 8.84360E-06, 7.94928E-06,   F37020
     * 6.68425E-06, 8.56576E-06, 1.05282E-05, 1.10647E-05, 9.91625E-06,   F37030
     * 7.95356E-06, 8.66443E-06, 9.13551E-06, 1.04870E-05, 9.79244E-06,   F37040
     * 1.26214E-05, 8.42148E-06, 8.13468E-06, 1.11338E-05, 1.06780E-05,   F37050
     * 8.54578E-06, 7.82119E-06, 8.33258E-06, 8.23644E-06, 5.95583E-06/   F37060
      DATA C20721 /                                                       F37070
     * 5.85592E-06, 4.05898E-06, 6.39260E-06, 8.43280E-06, 8.76251E-06,   F37080
     * 6.70423E-06, 6.81368E-06, 7.43506E-06, 7.14376E-06, 6.51065E-06,   F37090
     * 5.65633E-06, 5.42394E-06, 7.10817E-06, 4.78831E-06, 6.29380E-06,   F37100
     * 4.87344E-06, 6.81764E-06, 6.51611E-06, 5.70526E-06, 6.50590E-06,   F37110
     * 6.61568E-06, 5.39248E-06, 6.32002E-06, 7.98976E-06, 7.73795E-06,   F37120
     * 4.85788E-06, 5.83443E-06, 6.11694E-06, 5.40408E-06, 5.00946E-06,   F37130
     * 5.62153E-06, 6.30263E-06, 6.05764E-06, 5.53274E-06, 5.80664E-06,   F37140
     * 5.18684E-06, 6.85555E-06, 6.22889E-06, 6.06959E-06, 6.49228E-06,   F37150
     * 5.64064E-06, 4.92690E-06, 5.77661E-06, 7.18450E-06, 7.38658E-06,   F37160
     * 6.77379E-06, 5.74668E-06, 6.68355E-06, 6.13655E-06, 6.43266E-06,   F37170
     * 7.08896E-06, 7.71187E-06, 7.37273E-06, 6.75882E-06, 6.39307E-06,   F37180
     * 4.59520E-06, 5.10323E-06, 5.80178E-06, 6.88172E-06, 6.68825E-06,   F37190
     * 7.50416E-06, 6.14975E-06, 6.51422E-06, 7.74942E-06, 8.11492E-06,   F37200
     * 1.19607E-05, 7.92722E-06, 4.47848E-06, 6.02524E-06, 9.74067E-06,   F37210
     * 1.02429E-05, 8.60819E-06, 8.57044E-06, 1.09196E-05, 1.02048E-05,   F37220
     * 3.86222E-06, 9.26104E-06, 7.33341E-06, 9.08181E-06, 1.05569E-05/   F37230
      DATA C20801 /                                                       F37240
     * 1.06776E-05, 1.10247E-05, 1.04520E-05, 8.78328E-06, 7.60679E-06,   F37250
     * 7.27896E-06, 9.72776E-06, 5.16039E-06, 1.03134E-05, 1.09088E-05,   F37260
     * 8.12575E-06, 7.61685E-06, 8.16346E-06, 5.91269E-06, 3.61448E-06,   F37270
     * 8.74336E-06, 1.03990E-05, 6.25691E-06, 7.04541E-06, 7.94348E-06,   F37280
     * 8.39807E-06, 8.67342E-06, 8.32173E-06, 7.56015E-06, 8.31782E-06,   F37290
     * 6.36556E-06, 6.99328E-06, 6.24490E-06, 6.73080E-06, 6.95852E-06,   F37300
     * 7.55508E-06, 7.74168E-06, 7.90414E-06, 8.94934E-06, 7.99809E-06,   F37310
     * 6.12528E-06, 9.04115E-06, 7.14535E-06, 5.88625E-06, 6.43941E-06,   F37320
     * 7.11566E-06, 7.47425E-06, 8.23805E-06, 6.19919E-06, 7.31614E-06,   F37330
     * 8.24852E-06, 6.82172E-06, 5.45362E-06, 6.66115E-06, 8.44300E-06,   F37340
     * 8.07530E-06, 7.22735E-06, 5.85614E-06, 5.13900E-06, 6.03215E-06,   F37350
     * 6.59491E-06, 4.81592E-06, 4.48587E-06, 7.11525E-06, 8.36201E-06,   F37360
     * 7.11669E-06, 2.80033E-06, 6.50756E-06, 9.43974E-06, 5.22402E-06,   F37370
     * 3.82334E-06, 7.29963E-06, 8.62313E-06, 7.42018E-06, 4.56506E-06,   F37380
     * 5.29972E-06, 5.62787E-06, 4.63852E-06, 5.18329E-06, 7.01884E-06,   F37390
     * 7.24888E-06, 5.18157E-06, 5.40219E-06, 5.92412E-06, 4.97977E-06/   F37400
      DATA C20881 /                                                       F37410
     * 5.29040E-06, 5.33812E-06, 4.76620E-06, 4.65759E-06, 5.10546E-06,   F37420
     * 6.49525E-06, 4.43416E-06, 5.30223E-06, 3.27044E-06, 2.55324E-06,   F37430
     * 4.85017E-06, 7.46556E-06, 8.04448E-06, 5.14009E-06, 6.09755E-06,   F37440
     * 5.38381E-06, 6.41959E-06, 6.59233E-06, 4.83160E-06, 3.81289E-06,   F37450
     * 5.37013E-06, 5.69212E-06, 5.54983E-06, 5.73495E-06, 4.00639E-06,   F37460
     * 2.33817E-06, 2.55751E-06, 3.29627E-06, 3.59845E-06, 6.20623E-06,   F37470
     * 4.47088E-06, 3.49267E-06, 3.09273E-06, 3.32506E-06, 4.83353E-06,   F37480
     * 6.39001E-06, 3.78074E-06, 4.07848E-06, 4.01811E-06, 3.19767E-06,   F37490
     * 3.34053E-06, 4.34246E-06, 3.68003E-06, 3.01090E-06, 3.98545E-06,   F37500
     * 2.72338E-06, 1.90024E-06, 2.77553E-06, 3.73381E-06, 2.58685E-06,   F37510
     * 1.70987E-06,-5.48480E-07, 1.64591E-06, 2.43481E-06, 2.52116E-06,   F37520
     * 2.19316E-06, 1.32392E-06, 1.75370E-06, 2.65409E-07, 2.22278E-06,   F37530
     * 2.53079E-06, 2.87260E-06, 1.87600E-06,-3.84453E-07, 1.80836E-06,   F37540
     * 9.28123E-07, 1.94986E-06, 2.40483E-06, 2.79865E-06, 2.86361E-06,   F37550
     * 2.63868E-06, 3.34704E-06, 3.32132E-06, 2.58463E-06, 2.45684E-06,   F37560
     * 3.35043E-06, 3.19848E-06, 1.73037E-06, 2.98206E-06, 2.77491E-06/   F37570
      DATA C20961 /                                                       F37580
     * 6.51674E-07, 2.52219E-06, 2.97136E-06, 1.96700E-06, 2.29350E-06,   F37590
     * 3.01956E-06, 3.20811E-06, 1.30467E-06, 1.68172E-06, 2.56264E-06,   F37600
     * 2.46167E-06, 1.78221E-06, 2.31647E-06, 2.69480E-06, 2.63619E-06,   F37610
     * 1.81319E-06, 1.83448E-06, 2.23432E-06, 8.14045E-07, 8.75863E-07,   F37620
     * 1.61350E-06, 1.59796E-06, 2.08419E-06, 1.89665E-06, 6.93584E-07,   F37630
     * 1.09880E-06, 3.79031E-07,-3.36470E-07, 1.04326E-06, 1.06497E-06,   F37640
     * 2.15108E-07, 3.28774E-07,-5.17613E-07, 1.27762E-06, 8.22924E-07,   F37650
     * 4.92835E-07, 2.24698E-08,-1.99111E-07, 1.30262E-06,-3.81299E-07,   F37660
     * 9.55084E-07, 2.17641E-07,-6.03874E-08, 8.44121E-07, 1.72391E-06,   F37670
     * 1.66921E-06, 2.19855E-06, 1.17655E-06, 1.79637E-06, 3.31670E-06,   F37680
     * 3.40206E-06, 6.05670E-07, 2.08299E-06, 2.10121E-06, 1.68598E-06,   F37690
     * 2.21155E-06, 2.43221E-06, 5.81282E-08, 1.62613E-06,-5.49850E-07,   F37700
     * 2.14143E-07, 1.21751E-06, 2.30470E-06, 4.27911E-06, 2.96622E-06,   F37710
     * 8.67534E-07, 9.12041E-07, 2.48797E-06, 9.43519E-07,-3.60949E-06,   F37720
     * 2.01928E-06, 1.88873E-06, 8.06749E-07, 7.33519E-07, 1.17440E-06,   F37730
     * 1.69744E-06, 3.64492E-06, 3.11556E-06, 8.89471E-07, 1.93064E-06/   F37740
      DATA C21041 /                                                       F37750
     * 3.02787E-06, 1.92575E-06, 1.73720E-06,-1.32700E-07, 1.41743E-06,   F37760
     * 2.24632E-06, 2.47945E-06, 2.05151E-06,-9.56031E-07, 2.57317E-07,   F37770
     * 3.00980E-06, 3.07981E-06, 2.78202E-06, 3.02555E-06, 5.48784E-09,   F37780
     * 2.37693E-06, 2.90011E-06, 2.93608E-06, 2.14837E-06, 6.55832E-07,   F37790
     * 3.41155E-07,-2.13884E-06, 2.52553E-06, 4.27109E-06, 3.33766E-06,   F37800
     * 3.07708E-06, 2.66405E-06, 3.22850E-06,-5.78879E-07,-6.06194E-07,   F37810
     * 1.72864E-06, 1.57072E-06,-3.39701E-07, 7.21540E-08, 1.67012E-06,   F37820
     * 2.48568E-06, 2.70214E-06, 3.62383E-06, 2.20408E-06, 1.19395E-06,   F37830
     * 1.53825E-06, 2.37511E-06, 2.66754E-06, 1.77020E-06, 5.40420E-07,   F37840
     * 2.01156E-06, 3.27498E-06, 3.04720E-06, 1.96213E-06, 3.71633E-06,   F37850
     * 2.07886E-06, 1.60069E-06, 5.33370E-07, 1.33966E-07, 2.16073E-06,   F37860
     * 8.81457E-07, 1.12880E-06, 2.40509E-06, 2.94252E-06, 2.22899E-06,   F37870
     * 1.80941E-06, 2.68577E-06, 2.44584E-06, 2.51720E-06, 2.64857E-06,   F37880
     * 2.24182E-06, 1.62007E-06, 2.60421E-06, 3.09782E-06, 3.11099E-06,   F37890
     * 3.81513E-06, 6.91606E-06, 3.28767E-06, 3.44175E-06, 4.16771E-06,   F37900
     * 3.75452E-06, 2.21050E-06, 2.99939E-06, 2.86993E-06, 2.47080E-06/   F37910
      DATA C21121 /                                                       F37920
     * 2.33607E-06, 2.68568E-06, 3.39344E-06, 6.09518E-06, 5.10422E-06,   F37930
     * 4.04027E-06, 4.01363E-06, 4.53142E-06, 2.94424E-06, 4.76694E-06,   F37940
     * 6.44206E-06, 7.86435E-06, 8.55564E-06, 6.00857E-06, 5.48073E-06,   F37950
     * 1.56287E-06,-1.16619E-06,-1.85215E-06,-3.04762E-06,-3.45420E-07,   F37960
     * 2.48111E-07,-1.39302E-07,-6.27593E-07,-5.26792E-07, 4.81454E-08,   F37970
     *-3.08631E-08,-1.02976E-06,-1.54919E-06,-9.34044E-07,-1.02507E-06,   F37980
     *-1.39794E-06,-1.15709E-06,-1.04875E-06,-1.64379E-06,-2.97514E-06,   F37990
     *-3.22236E-07,-1.18978E-06,-2.85325E-06,-3.93143E-06,-4.15349E-06,   F38000
     *-2.33228E-06,-3.27125E-06,-2.44987E-06,-1.44460E-06,-3.59727E-06,   F38010
     *-7.18516E-07,-1.53237E-06,-1.53526E-06,-1.56450E-06,-2.91088E-06,   F38020
     *-8.52134E-07,-1.44575E-07,-1.50350E-06,-2.92806E-06,-2.47710E-06,   F38030
     *-9.71202E-07,-9.82754E-07,-1.09924E-06,-6.08199E-07, 3.62885E-07,   F38040
     *-6.67372E-07,-1.00033E-06,-1.12001E-06,-1.06624E-06,-9.23789E-07,   F38050
     *-9.83788E-07,-2.11656E-06,-2.45001E-06,-2.75874E-06,-3.36003E-06,   F38060
     *-3.38364E-06,-2.63747E-06,-3.11047E-06,-3.75258E-06,-3.83211E-06,   F38070
     *-3.52833E-06,-3.48464E-06,-3.77021E-06,-4.26887E-06,-4.23917E-06/   F38080
      DATA C21201 /                                                       F38090
     *-1.42438E-06,-2.48477E-06,-2.84719E-06,-2.70247E-06,-2.50588E-06,   F38100
     *-2.22900E-06,-1.78393E-06,-1.76826E-06,-2.16396E-06,-2.67543E-06,   F38110
     *-2.23706E-06,-2.31793E-06,-2.87590E-06,-3.07803E-06,-2.50493E-06,   F38120
     *-4.54223E-06,-5.15511E-06,-5.39690E-06,-4.89633E-06,-3.33710E-06,   F38130
     *-4.56583E-06,-4.78877E-06,-3.93508E-06,-3.29027E-06,-4.95668E-06,   F38140
     *-6.01801E-06,-5.76016E-06,-5.34657E-06,-5.29080E-06,-5.57133E-06,   F38150
     *-5.73135E-06,-5.39374E-06,-5.09808E-06,-5.12874E-06,-5.20269E-06,   F38160
     *-7.30702E-06,-7.04220E-06,-5.96514E-06,-5.74802E-06,-4.53961E-06,   F38170
     *-4.42127E-06,-4.63922E-06,-4.80622E-06,-4.69659E-06,-5.96786E-06,   F38180
     *-6.29800E-06,-4.75452E-06,-2.85907E-06,-5.33662E-06,-5.31681E-06,   F38190
     *-5.04646E-06,-5.21729E-06,-5.93409E-06,-5.73462E-06,-5.44926E-06,   F38200
     *-6.43325E-06,-7.74451E-06,-7.83147E-06,-5.51568E-06,-7.37048E-06,   F38210
     *-4.25726E-06, 2.32917E-06,-5.61131E-07, 2.05234E-06, 3.74631E-07,   F38220
     *-7.66493E-07, 1.42689E-06,-7.79683E-07, 9.06809E-07, 5.13642E-07,   F38230
     *-1.52504E-06,-2.12058E-06,-2.50316E-06, 1.03637E-08, 5.60002E-07,   F38240
     *-1.48075E-06, 1.94155E-06, 1.91846E-06, 2.78507E-06, 3.90146E-06/   F38250
      DATA C21281 /                                                       F38260
     * 3.61409E-06, 3.23677E-06, 4.00022E-06, 3.19157E-06, 4.03034E-07,   F38270
     *-2.03929E-06, 1.23366E-06, 3.28589E-06, 3.94168E-06, 3.94672E-06,   F38280
     * 3.84619E-06, 2.30400E-07,-2.07799E-06,-1.75115E-06,-5.71958E-07,   F38290
     * 2.33425E-06, 2.01664E-06, 6.05673E-07, 9.57363E-07,-8.89924E-07,   F38300
     *-4.71331E-07, 2.82826E-07, 5.10859E-07, 3.63512E-07, 9.86288E-07,   F38310
     *-4.86309E-07,-2.23163E-06,-1.23370E-06,-2.43131E-07,-2.11498E-06,   F38320
     *-1.56756E-06, 2.70905E-06, 1.87606E-08, 7.83721E-08, 1.58444E-06,   F38330
     * 2.88574E-06, 1.40306E-06, 2.40883E-06, 2.84063E-06, 3.13820E-06,   F38340
     * 3.71016E-06, 3.12975E-06, 3.21981E-06, 2.56191E-06, 1.04624E-06,   F38350
     * 1.87464E-07, 7.25329E-07, 1.03650E-06, 7.23663E-07,-4.18739E-07,   F38360
     * 9.95744E-07,-1.80878E-07,-1.04044E-06, 3.86965E-07,-9.36186E-07,   F38370
     *-4.02271E-07,-2.00231E-07,-5.94965E-07, 4.94038E-07, 3.34585E-07,   F38380
     * 4.82255E-07, 1.12599E-06, 2.11763E-06, 2.66807E-07, 2.29324E-07,   F38390
     * 7.07005E-07, 3.41907E-07,-1.17115E-07, 9.03089E-07, 1.76844E-06,   F38400
     * 1.87134E-06, 2.64057E-06, 4.00395E-07,-4.19679E-07, 6.30769E-07,   F38410
     * 1.02725E-06, 1.05876E-06,-4.08660E-07,-2.32668E-06,-2.73468E-06/   F38420
      DATA C21361 /                                                       F38430
     *-2.40600E-06,-1.81203E-06,-7.96431E-07, 7.40789E-07, 2.73188E-07,   F38440
     * 1.68367E-07,-1.27227E-07,-1.05041E-06,-3.51726E-06,-1.64956E-06,   F38450
     *-5.63840E-07,-1.61242E-06,-1.33264E-06, 1.56604E-06, 2.35083E-06,   F38460
     * 9.26708E-07, 5.41983E-07, 3.54277E-07, 8.53743E-07, 1.54196E-06,   F38470
     * 1.19902E-06, 1.10552E-06, 1.63179E-06, 1.96366E-06, 7.82848E-07,   F38480
     *-3.34741E-08,-7.90842E-07,-6.45131E-07, 1.36158E-06, 1.62453E-06,   F38490
     * 6.68965E-07,-4.86203E-08, 6.83561E-07, 1.89652E-06,-2.80988E-07,   F38500
     *-2.30536E-06,-1.90777E-06, 1.31617E-06, 1.27309E-06, 5.90825E-07,   F38510
     * 5.65686E-07, 1.23631E-07,-1.70279E-06,-1.60768E-06, 9.69543E-07,   F38520
     * 1.01108E-07,-2.02473E-06,-1.75146E-06, 6.33201E-07,-3.59110E-06,   F38530
     *-9.71706E-07, 9.16822E-07, 1.40681E-07,-7.16745E-07,-2.11376E-06,   F38540
     *-1.00951E-06, 2.12465E-06, 1.06982E-06, 1.44032E-06, 1.49692E-06,   F38550
     * 1.07277E-06, 1.37006E-06, 1.66932E-06, 1.75820E-06, 1.41859E-06,   F38560
     *-5.84947E-08, 2.17349E-06, 4.27053E-06, 5.27286E-06, 5.87085E-06,   F38570
     * 2.42692E-06, 2.39305E-06, 6.19573E-06, 5.12518E-06, 1.27171E-06,   F38580
     *-6.81963E-07, 4.16199E-08,-1.36608E-06,-2.53272E-06,-2.37700E-06/   F38590
      DATA C21441 /                                                       F38600
     *-7.96719E-07, 3.85367E-07,-1.08393E-07,-9.04587E-07,-1.54917E-06,   F38610
     *-3.11945E-06,-5.58484E-07, 1.61347E-06, 1.11736E-06, 2.11889E-06,   F38620
     * 2.43534E-06, 1.46709E-06,-1.05429E-06, 1.09978E-06, 7.22493E-07,   F38630
     * 8.53307E-08, 1.22733E-06, 2.99380E-06, 3.62416E-06, 3.81404E-06,   F38640
     * 4.46735E-06, 4.70753E-06, 4.54494E-06, 3.83002E-06, 2.28067E-06,   F38650
     * 2.03102E-06, 2.43844E-06, 2.93132E-06, 2.17555E-06, 3.92919E-06,   F38660
     * 3.53089E-06, 1.61388E-06, 5.09498E-06, 3.40067E-06, 1.58876E-06,   F38670
     * 1.17367E-06, 1.13344E-06, 1.17798E-06, 1.10976E-06, 7.90635E-07,   F38680
     *-4.15989E-07,-1.00581E-06,-9.60236E-07,-1.79111E-07,-5.70733E-07,   F38690
     * 1.49766E-06, 3.44374E-06, 6.45914E-07, 1.00532E-06, 2.01068E-06,   F38700
     * 2.59092E-06, 9.35770E-08, 6.00121E-07, 1.54409E-06, 2.03537E-06,   F38710
     * 8.10358E-07, 1.34126E-06, 1.88873E-06, 1.43283E-06,-2.05029E-07,   F38720
     *-1.09782E-06,-6.56149E-07, 2.01650E-06, 1.84770E-06, 4.39586E-08,   F38730
     *-2.03588E-06,-1.46366E-06,-3.45189E-07, 4.02577E-07, 3.10362E-07,   F38740
     *-2.16073E-06,-1.91861E-06,-2.90520E-07, 2.03692E-06, 3.47996E-06,   F38750
     * 4.21761E-06, 3.89000E-06, 1.86138E-06, 1.56143E-06, 4.88964E-07/   F38760
      DATA C21521 /                                                       F38770
     *-9.28184E-07,-4.34315E-07, 8.74954E-07, 1.58417E-06, 1.36880E-06,   F38780
     * 2.65016E-06, 4.62623E-06, 5.81990E-06, 4.72139E-06, 1.95905E-06,   F38790
     * 1.54151E-06, 2.95768E-06, 4.71536E-06, 2.62359E-06, 9.11513E-07,   F38800
     * 4.75677E-07,-1.53801E-06,-2.32382E-06,-2.25220E-06,-1.46641E-06,   F38810
     *-2.23014E-06,-2.12604E-06,-1.66259E-06,-2.48856E-06,-2.38895E-06,   F38820
     *-2.18158E-06,-1.95841E-06, 4.43899E-07, 1.08517E-06, 1.66370E-07,   F38830
     *-2.42342E-06,-7.19331E-07, 3.19532E-07, 3.58690E-07,-2.01979E-07,   F38840
     * 5.07242E-07, 1.10562E-06, 1.00419E-06, 1.22379E-06, 7.05180E-07,   F38850
     * 1.42283E-07, 8.61092E-07, 8.95236E-07, 1.18043E-07,-1.23589E-06,   F38860
     *-6.16316E-07,-1.18947E-06,-1.45838E-06,-1.47522E-09, 1.33867E-06,   F38870
     * 9.18310E-07,-8.98949E-07,-2.27314E-06,-1.71510E-06,-7.16704E-07,   F38880
     * 8.60666E-09, 5.68015E-07, 1.31219E-06, 1.75478E-06, 5.11790E-07,   F38890
     * 3.35270E-07, 5.39243E-07, 9.08467E-07, 1.39382E-06, 1.08806E-06,   F38900
     * 1.18589E-06, 3.58461E-06, 2.78668E-06, 1.25964E-06,-2.72255E-07,   F38910
     * 1.72305E-06, 1.82937E-06, 7.46252E-07,-1.10555E-06, 2.24967E-07,   F38920
     * 6.45674E-07,-1.87591E-07,-8.84068E-07,-1.75433E-06,-2.17670E-06/   F38930
      DATA C21601 /                                                       F38940
     *-1.37112E-06,-2.31722E-06,-2.23860E-06,-1.16796E-06,-2.23765E-06,   F38950
     *-1.86406E-06,-1.03517E-06,-5.90824E-07,-6.57710E-07,-7.00941E-07,   F38960
     *-4.46064E-07, 1.77205E-06, 2.45066E-06, 2.39371E-06, 2.30736E-06,   F38970
     * 2.35355E-06, 1.85070E-06, 9.62711E-07, 2.59644E-06, 2.05304E-06,   F38980
     * 9.70090E-07, 1.50942E-06, 3.79439E-06, 2.94597E-06,-1.91789E-06,   F38990
     * 6.44324E-08,-3.92094E-07,-1.55398E-06, 4.46701E-08,-4.78760E-07,   F39000
     *-1.70061E-06,-3.17252E-06,-2.93173E-06,-2.01455E-06,-7.76298E-07,   F39010
     *-2.74577E-07,-1.39907E-06,-2.16470E-06,-1.26010E-06,-2.76845E-06,   F39020
     *-2.38226E-06,-5.49068E-08, 9.65258E-07, 1.08650E-06, 5.64738E-07,   F39030
     *-5.78379E-07,-5.68918E-07,-1.90177E-06,-5.08874E-06,-3.03648E-06,   F39040
     *-1.30527E-06,-4.87669E-07,-2.83326E-06,-1.97823E-06,-5.94313E-07,   F39050
     *-1.50961E-07,-1.15908E-06,-1.43260E-06,-9.29331E-07,-1.39459E-06,   F39060
     *-1.27237E-06,-1.50189E-06,-3.79292E-06,-3.92038E-06,-3.58490E-06,   F39070
     *-3.26439E-06,-2.42138E-06,-2.70516E-06,-3.58080E-06,-1.71822E-06,   F39080
     *-2.41567E-06,-3.50193E-06,-2.62394E-06,-3.08424E-06,-3.89604E-06,   F39090
     *-4.84127E-06,-4.41385E-06,-3.22673E-06,-1.80987E-06,-2.93027E-06/   F39100
      DATA C21681 /                                                       F39110
     *-3.17366E-06,-2.79721E-06,-1.78848E-06,-2.80254E-06,-3.55572E-06,   F39120
     *-3.34632E-06,-2.83979E-06,-2.48022E-06,-2.15090E-06,-1.08311E-06,   F39130
     *-6.15216E-07,-7.13008E-07,-1.70841E-06,-2.96098E-06,-3.57134E-06,   F39140
     *-3.04405E-06,-3.35280E-06,-2.97780E-06,-1.97966E-06,-2.33197E-06,   F39150
     *-2.76708E-06,-2.70409E-06,-4.51094E-07,-1.43068E-06,-2.83719E-06,   F39160
     *-2.98921E-06,-4.14949E-06,-3.63780E-06,-8.10138E-07,-1.61597E-06,   F39170
     *-2.25394E-06,-2.58110E-06,-1.57781E-06,-1.71520E-06,-2.30016E-06,   F39180
     *-2.61268E-06,-1.96696E-06,-1.86744E-06,-3.15645E-06,-3.59354E-06,   F39190
     *-3.61015E-06,-3.21793E-06,-2.57436E-06,-2.74347E-06,-3.33319E-06,   F39200
     *-2.93400E-06,-3.25986E-06,-3.46384E-06,-2.22114E-06,-2.92650E-06,   F39210
     *-3.73666E-06,-3.70485E-06,-2.75963E-06,-2.40652E-06,-2.93107E-06,   F39220
     *-1.77517E-06,-1.57096E-06,-2.17533E-06,-2.80190E-06,-2.27942E-06,   F39230
     *-1.37371E-06,-1.65974E-06,-1.26079E-06,-8.08050E-07,-8.41278E-07,   F39240
     *-1.53860E-06,-1.66687E-06,-6.56592E-07,-3.05110E-08, 1.08623E-07,   F39250
     *-2.87222E-07,-2.63555E-07,-7.89575E-07,-1.56059E-06,-6.42174E-07,   F39260
     *-9.43333E-07,-1.38671E-06, 6.50443E-07, 1.35301E-06, 9.27981E-07/   F39270
      DATA C21761 /                                                       F39280
     *-1.21705E-06,-9.63848E-08, 8.73593E-07,-3.47278E-08,-1.79042E-06,   F39290
     *-2.15544E-06,-4.48668E-07,-1.17414E-06,-1.35437E-06,-8.90688E-07,   F39300
     *-4.54757E-07, 2.41484E-09, 3.88010E-07,-1.85349E-08, 1.58011E-07,   F39310
     * 3.70566E-07,-7.30268E-07,-8.42354E-07,-4.13738E-07, 3.96796E-07,   F39320
     *-5.55763E-07,-1.26877E-06,-2.89854E-07, 5.78676E-07, 9.51356E-07,   F39330
     * 5.56912E-07, 1.05014E-06, 9.75896E-07, 5.91573E-08,-6.15073E-07,   F39340
     *-1.48803E-06,-2.53397E-06,-1.77027E-06,-2.08546E-06,-3.10452E-06,   F39350
     *-1.65227E-06,-1.15981E-06,-1.25849E-06,-9.65711E-07,-1.90319E-06,   F39360
     *-2.71831E-06,-5.71559E-08,-1.20368E-06,-3.16820E-06,-2.22766E-06,   F39370
     *-1.19828E-06,-2.82573E-07, 2.53850E-07,-9.10547E-07,-1.65529E-06,   F39380
     *-6.00138E-07,-4.98898E-07,-3.45799E-07, 2.25160E-07, 1.14332E-07,   F39390
     * 3.16082E-07, 1.12681E-06,-6.04876E-07,-7.24616E-07, 1.48177E-06,   F39400
     * 1.05680E-06, 5.91076E-07, 2.07187E-07, 3.82385E-07, 5.91560E-07,   F39410
     * 8.26519E-07, 1.22139E-06, 1.63501E-06, 2.06423E-06, 2.50038E-06,   F39420
     * 2.38037E-06, 1.91688E-06, 2.46702E-06, 2.45066E-06, 2.16732E-06,   F39430
     * 3.13517E-06, 2.68221E-06, 1.39877E-06, 8.58945E-07, 6.83181E-07/   F39440
      DATA C21841 /                                                       F39450
     * 8.46816E-07, 1.73491E-06, 1.98732E-06, 1.94059E-06, 2.19284E-06,   F39460
     * 1.73215E-06, 1.06865E-06, 1.14117E-06, 1.43213E-06, 1.42275E-06,   F39470
     *-4.15449E-07,-2.39911E-07, 3.46498E-08,-2.75022E-06,-2.43736E-06,   F39480
     *-1.06489E-06,-7.81941E-07,-8.04801E-07,-1.04984E-06,-1.65734E-06,   F39490
     *-1.03167E-06,-3.18255E-08, 5.70283E-07, 6.19050E-07, 2.92257E-07,   F39500
     *-6.01436E-07,-7.04005E-07,-3.70875E-07, 4.12830E-07, 1.31319E-07,   F39510
     *-1.61570E-07, 9.76170E-07, 7.99907E-07, 1.41860E-07,-1.98022E-07,   F39520
     * 3.13766E-07, 7.43982E-07,-6.11287E-07,-5.21146E-07, 1.11156E-07,   F39530
     * 3.91719E-07, 5.45566E-07, 6.39059E-07, 7.29515E-07, 4.59167E-07,   F39540
     * 6.13179E-08,-3.48146E-08, 5.32046E-07, 1.19736E-06, 3.83982E-07,   F39550
     * 1.73267E-07, 3.54304E-07, 9.34657E-07, 5.53819E-07,-2.86678E-07,   F39560
     * 2.01853E-08,-1.56159E-07,-6.08130E-07,-2.14929E-07, 1.66317E-08,   F39570
     * 9.32462E-08,-4.83623E-07,-9.16323E-07,-1.22772E-06,-1.61586E-06,   F39580
     *-1.27409E-06,-1.98119E-07,-3.69182E-08,-1.41061E-07,-5.12562E-07,   F39590
     *-4.55495E-07,-8.12132E-07,-1.71772E-06,-2.70741E-06,-2.98751E-06,   F39600
     *-2.19520E-06, 3.01900E-07, 1.17806E-06,-1.23067E-06, 4.17086E-07/   F39610
      DATA C21921 /                                                       F39620
     * 1.68113E-06, 4.81677E-07,-1.55187E-07,-3.35287E-07, 2.94916E-07,   F39630
     * 4.57124E-07, 3.38692E-07,-2.49203E-07,-3.62585E-07,-2.39653E-07,   F39640
     * 3.72675E-08,-7.79964E-09,-2.83285E-07,-9.74713E-07,-6.91171E-07,   F39650
     * 1.21925E-07, 3.39940E-07, 3.68441E-08,-5.82188E-07, 2.12605E-07,   F39660
     * 4.65144E-07, 2.17190E-07, 7.50119E-07, 8.62008E-07, 4.63016E-07,   F39670
     * 1.25620E-06, 1.04567E-06,-8.17037E-07,-1.20023E-06,-1.06224E-06,   F39680
     *-3.77100E-07,-1.28057E-07,-2.76183E-07,-1.24304E-06,-2.56776E-06,   F39690
     *-3.36699E-06,-1.49408E-06,-1.01189E-07, 7.41870E-07,-6.45425E-07,   F39700
     *-7.47111E-07, 4.79055E-10,-1.32339E-06,-1.86135E-06,-1.61074E-06,   F39710
     *-1.82039E-06,-1.68040E-06,-1.08025E-06,-8.61965E-07,-7.00131E-07,   F39720
     *-5.63105E-07,-8.09843E-07,-8.09221E-07, 1.69474E-07,-1.33941E-07,   F39730
     *-7.49558E-07,-5.19013E-07,-8.53534E-07,-1.33703E-06,-3.11161E-07,   F39740
     * 8.99037E-07, 2.25330E-06, 1.44822E-06, 3.07437E-07,-1.22366E-06,   F39750
     *-7.64217E-07, 2.13156E-08, 1.07909E-06, 6.10755E-07, 1.81483E-07,   F39760
     * 8.12405E-07,-9.13283E-08,-1.35885E-06,-1.58366E-06,-7.88594E-07,   F39770
     * 4.48283E-07,-1.23754E-06,-1.65105E-06,-8.93014E-07,-1.48622E-06/   F39780
      DATA C22001 /                                                       F39790
     *-1.67948E-06,-1.24310E-06,-1.54411E-06,-1.65677E-06,-1.04998E-06,   F39800
     *-1.46985E-07, 4.61778E-07,-4.87832E-07,-4.89452E-07,-1.24840E-07,   F39810
     *-1.70101E-06,-1.66976E-06,-1.48528E-07,-1.12621E-07,-2.30607E-08,   F39820
     * 1.82301E-07,-8.58152E-07,-1.89794E-06,-2.46464E-06,-2.32745E-06,   F39830
     *-2.02112E-06,-2.07656E-06,-1.43824E-06,-5.16583E-07,-1.80702E-06,   F39840
     *-2.93490E-06,-3.89216E-06,-3.36211E-06,-2.41393E-06,-9.53406E-07,   F39850
     *-1.16269E-06,-1.66431E-06,-1.77150E-06,-1.82496E-06,-1.93095E-06,   F39860
     *-2.75759E-06,-2.83618E-06,-2.27908E-06,-6.33348E-07, 5.61257E-07,   F39870
     * 1.00142E-06, 7.73337E-07, 3.17721E-07,-3.69804E-07,-8.82058E-07,   F39880
     *-1.17364E-06,-4.53480E-07,-2.47824E-07,-4.79624E-07,-5.17032E-07,   F39890
     *-3.46498E-07, 1.42669E-07,-1.59168E-07,-5.06580E-07,-3.18573E-07,   F39900
     *-2.74092E-07,-2.68860E-07, 1.32811E-07,-2.35567E-09,-6.71971E-07,   F39910
     *-9.75302E-07,-8.70978E-07,-3.59071E-08,-3.01726E-07,-8.27641E-07,   F39920
     *-1.14899E-06,-1.50160E-06,-1.83660E-06,-1.26290E-06,-1.07659E-06,   F39930
     *-1.34878E-06,-5.24626E-07,-7.85094E-08,-8.79473E-07,-1.19291E-06,   F39940
     *-1.33298E-06,-1.59750E-06,-1.31836E-06,-5.73079E-07,-1.10349E-06/   F39950
      DATA C22081 /                                                       F39960
     *-1.11807E-06,-1.99530E-07,-8.10496E-07,-1.42679E-06,-5.34617E-07,   F39970
     *-2.05001E-07,-2.51690E-07,-1.01740E-06,-1.02841E-06,-7.48750E-08,   F39980
     *-1.01770E-06,-1.50413E-06, 1.80898E-07, 3.63788E-07,-1.97900E-07,   F39990
     *-1.16721E-06,-1.05497E-06,-2.07218E-08,-1.90590E-07,-8.25501E-07,   F40000
     *-2.21142E-06,-1.19905E-06, 2.16271E-07,-2.52574E-07,-4.35837E-07,   F40010
     *-3.95272E-07, 5.97065E-08, 2.76639E-07, 9.22569E-08, 1.20142E-07,   F40020
     *-2.95030E-09,-1.08216E-06,-1.32386E-06,-9.62248E-07,-1.99430E-06,   F40030
     *-2.13890E-06,-9.56082E-07,-6.94022E-07,-7.75721E-07,-1.31048E-06,   F40040
     *-1.50080E-06,-1.35873E-06,-7.48378E-07,-4.83436E-07,-4.69624E-07,   F40050
     *-1.51156E-06,-2.48221E-06,-3.30134E-06,-2.79114E-06,-2.08976E-06,   F40060
     *-2.24768E-06,-1.06947E-06, 1.17462E-06,-2.51423E-07,-7.85729E-07,   F40070
     * 5.37467E-07,-9.39876E-08,-1.11303E-06,-7.46860E-07,-9.36220E-07,   F40080
     *-1.59880E-06,-1.61420E-06,-1.54368E-06,-1.41036E-06,-7.20350E-07,   F40090
     * 1.35544E-07, 3.14481E-07, 6.29265E-07, 1.09161E-06,-1.36044E-07,   F40100
     *-1.22932E-06,-1.29847E-06,-3.26429E-06,-6.01062E-06,-2.09945E-06,   F40110
     * 1.26878E-07,-2.88050E-08,-6.82802E-07,-1.39340E-06,-1.82986E-06/   F40120
      DATA C22161 /                                                       F40130
     *-1.67208E-06,-1.07994E-06,-1.89195E-06,-2.10782E-06,-1.04519E-06,   F40140
     *-3.27672E-07, 1.95516E-07, 1.63838E-07,-2.29575E-07,-1.01609E-06,   F40150
     *-2.19286E-06,-2.71850E-06,-9.77485E-07,-1.48830E-06,-3.37826E-06,   F40160
     *-1.59130E-06,-5.74498E-07,-8.27962E-07,-9.92211E-07,-1.14422E-06,   F40170
     *-1.41420E-06,-1.11629E-06,-2.51575E-07, 1.60805E-07, 1.82934E-07,   F40180
     *-7.28868E-07,-2.57062E-07, 1.06520E-06, 4.16488E-07, 2.97049E-08,   F40190
     * 6.62797E-08, 8.29435E-07, 1.29657E-06,-2.27961E-06,-3.40386E-06,   F40200
     *-1.88594E-06,-2.29732E-06,-2.72594E-06,-2.09847E-06,-1.31771E-06,   F40210
     *-4.23693E-07,-4.96348E-07,-9.40209E-07,-2.08707E-06,-1.21368E-06,   F40220
     * 4.79409E-07,-1.12548E-08,-4.57316E-07,-8.40885E-07,-5.03210E-07,   F40230
     *-1.61036E-07,-1.05835E-06,-1.66417E-06,-1.97827E-06,-1.63737E-06,   F40240
     *-1.11711E-06,-3.16081E-07,-6.81746E-07,-1.82599E-06,-1.12895E-06,   F40250
     *-9.19712E-07,-1.91707E-06,-2.14767E-06,-2.03629E-06,-2.86441E-06,   F40260
     *-3.07735E-06,-2.28656E-06,-1.40256E-06,-5.50649E-07,-3.11627E-07,   F40270
     *-7.90261E-07,-2.10728E-06,-1.89739E-06,-1.53762E-06,-2.39947E-06,   F40280
     *-2.28765E-06,-1.27564E-06,-2.15154E-06,-3.17932E-06,-3.84234E-06/   F40290
      DATA C22241 /                                                       F40300
     *-3.65102E-06,-2.84055E-06,-2.48744E-06,-2.27683E-06,-2.33087E-06,   F40310
     *-3.44460E-06,-5.19613E-06,-2.85882E-06,-1.39921E-06,-2.00579E-06,   F40320
     *-2.80593E-06,-3.65940E-06,-2.39526E-06,-1.70389E-06,-2.03532E-06,   F40330
     *-2.71522E-06,-3.42227E-06,-2.23606E-06,-1.77845E-06,-2.42071E-06,   F40340
     *-2.61515E-06,-2.56413E-06,-1.49601E-06,-1.23245E-06,-2.08440E-06,   F40350
     *-2.11121E-06,-1.93424E-06,-2.27439E-06,-2.58183E-06,-2.84705E-06,   F40360
     *-2.32183E-06,-1.80966E-06,-3.04089E-06,-3.14334E-06,-1.91331E-06,   F40370
     *-1.51037E-06,-1.43610E-06,-2.11316E-06,-2.45184E-06,-2.42262E-06/   F40380
C                                                                         F40390
      END                                                                 F40400
C
C     --------------------------------------------------------------
C
      SUBROUTINE O3HHUV (V1C,V2C,DVC,NPTC,C)                              F40410
C                                                                         F40420
      IMPLICIT REAL*8           (V)                                     ! F40430
C                                                                         F40440
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F40450
      COMMON /O3HUV/ V1S,V2S,DVS,NPTS,S(133)                              F40460
      DIMENSION C(*)                                                      F40470
C                                                                         F40480
      DVC = DVS                                                           F40490
      V1C = V1ABS-DVC                                                     F40500
      V2C = V2ABS+DVC                                                     F40510
C                                                                         F40520
      I1 = (V1C-V1S)/DVS                                                  F40530
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      DO 10 J = 1, NPTC                                                   F40600
         I = I1+J                                                         F40610
         C(J) = 0.                                                        F40620
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F40630
         VJ = V1C+DVC* REAL(J-1)                                          F40640
         C(J) = S(I)/VJ                                                   F40650
C                                                                         F40660
C     RADIATION FLD REMOVED FROM U.V.    OZONE                            F40670
C                                                                         F40680
   10 CONTINUE                                                            F40690
C                                                                         F40700
      RETURN                                                              F40710
C                                                                         F40720
      END                                                                 F40730
C
C     --------------------------------------------------------------
C
      BLOCK DATA BO3HUV                                                   F40740
C                                                                         F40750
      IMPLICIT REAL*8           (V)                                     ! F40760
C                                                                         F40770
C     DATA DERIVED FROM MOLINA & MOLINA, JGR,91,14501-14508,1986.         F40780
C     VALUES BETWEEN 245 AND 185NM (40800 AND 54054CM-1) USED AS          F40790
C     DIRECT AVERAGE WITH NO TEMPERATURE DEPENDENCE.                      F40800
C                                                                         F40810
      COMMON /O3HUV/ V1C,V2C,DVC,NC,                                      F40820
     *               C02281(80),C02361(53)                                F40830
C                                                                         F40840
      DATA V1C /40800./, V2C /54000./ ,DVC /100./, NC /133/               F40850
C                                                                         F40860
      DATA C02281/                                                        F40870
     * 9.91204E-18, 9.76325E-18, 9.72050E-18, 9.51049E-18, 9.23530E-18,   F40880
     * 9.02306E-18, 8.90510E-18, 8.60115E-18, 8.39094E-18, 8.27926E-18,   F40890
     * 7.95525E-18, 7.73583E-18, 7.55018E-18, 7.31076E-18, 7.10415E-18,   F40900
     * 6.87747E-18, 6.66639E-18, 6.39484E-18, 6.27101E-18, 6.01019E-18,   F40910
     * 5.77594E-18, 5.60403E-18, 5.40837E-18, 5.21289E-18, 4.99329E-18,   F40920
     * 4.81742E-18, 4.61608E-18, 4.45707E-18, 4.28261E-18, 4.09672E-18,   F40930
     * 3.93701E-18, 3.77835E-18, 3.61440E-18, 3.45194E-18, 3.30219E-18,   F40940
     * 3.15347E-18, 3.01164E-18, 2.87788E-18, 2.74224E-18, 2.61339E-18,   F40950
     * 2.48868E-18, 2.36872E-18, 2.25747E-18, 2.14782E-18, 2.03997E-18,   F40960
     * 1.94281E-18, 1.84525E-18, 1.75275E-18, 1.67151E-18, 1.58813E-18,   F40970
     * 1.50725E-18, 1.43019E-18, 1.35825E-18, 1.28878E-18, 1.22084E-18,   F40980
     * 1.15515E-18, 1.09465E-18, 1.03841E-18, 9.83780E-19, 9.31932E-19,   F40990
     * 8.83466E-19, 8.38631E-19, 7.96631E-19, 7.54331E-19, 7.13805E-19,   F41000
     * 6.78474E-19, 6.44340E-19, 6.13104E-19, 5.81777E-19, 5.53766E-19,   F41010
     * 5.27036E-19, 5.03555E-19, 4.82633E-19, 4.61483E-19, 4.42014E-19,   F41020
     * 4.23517E-19, 4.07774E-19, 3.93060E-19, 3.80135E-19, 3.66348E-19/   F41030
      DATA C02361/                                                        F41040
     * 3.53665E-19, 3.47884E-19, 3.39690E-19, 3.34288E-19, 3.29135E-19,   F41050
     * 3.23104E-19, 3.18875E-19, 3.16800E-19, 3.15925E-19, 3.12932E-19,   F41060
     * 3.12956E-19, 3.15522E-19, 3.14950E-19, 3.15924E-19, 3.19059E-19,   F41070
     * 3.23109E-19, 3.27873E-19, 3.33788E-19, 3.39804E-19, 3.44925E-19,   F41080
     * 3.50502E-19, 3.55853E-19, 3.59416E-19, 3.68933E-19, 3.78284E-19,   F41090
     * 3.86413E-19, 3.98049E-19, 4.04700E-19, 4.12958E-19, 4.23482E-19,   F41100
     * 4.31203E-19, 4.41885E-19, 4.52651E-19, 4.61492E-19, 4.70493E-19,   F41110
     * 4.80497E-19, 4.90242E-19, 4.99652E-19, 5.10316E-19, 5.21510E-19,   F41120
     * 5.32130E-19, 5.43073E-19, 5.56207E-19, 5.61756E-19, 5.66799E-19,   F41130
     * 5.85545E-19, 5.92409E-19, 5.96168E-19, 6.12497E-19, 6.20231E-19,   F41140
     * 6.24621E-19, 6.34160E-19, 6.43622E-19/                             F41150
C                                                                         F41160
      END                                                                 F41170
C
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

      DATA To/ 296./, xlosmt/ 2.68675e+19/
c
      xktfac = (1./To)-(1./T)
c     
c     correct formulation for consistency with LBLRTM:
c
      factor = (1.e+20 /xlosmt) 
c
c     A factor of 0.21, the mixing ration of oxygen, in the Thibault et al.
c     formulation is not included here.  This factor is in the column amount.
C                           
      DVC = DVS             
      V1C = V1ABS-DVC       
      V2C = V2ABS+DVC       
C                           
      I1 = (V1C-V1S)/DVS        
      IF (V1C.LT.V1S) I1 = -1 
C                                    
      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       
c
      do 10 j=1,nptc
         i = i1+j
         C(J) = 0.                                                        F41620
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F41630
         VJ = V1C+DVC* REAL(J-1)                                          F41640
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
c        The units of these continua coefficients are  1 / (amagat_O2*amagat_air)
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
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = I1-1                                          
C                                                                        
      V1C = V1S+DVS* REAL(I1)                                            
      I2 = (V2C-V1S)/DVS                                                 
      NPTC = I2-I1+3                                                     
      V2C = V1C+DVS* REAL(NPTC-1)                                        
c
      DO 10 J = 1, NPTC                                                   F34310
         I = I1+J                                                         F34320
         C(J) = 0.                                                        F34330
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F34340
         vj = v1c + dvc* REAL(j-1)
         C(J) = xo2inf1(I)/vj
   10 CONTINUE                                                            F34360
C                                                                         F34370
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
      NPTC = (v2c-v1c)/dvc + 3.
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
C        O2 continuum formulated by Greenblatt et al. over the spectral region
C        8797-29870 cm-1:  "Absorption Coefficients of Oxygen Between 
c        330 and 1140 nm, G.D. Green blatt, J.J. Orlando, J.B. Burkholder,
c        and A.R. Ravishabkara,  J. Geophys. Res., 95, 18577-18582, 1990. 
c
c        The units conversion  is to (cm^2/molec)/atm(o2)
c
c      these are the conditions reported in the paper by Greenblatt et al. for 
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
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = I1-1                                          
C                                                                        
      V1C = V1S+DVS* REAL(I1)                                            
      I2 = (V2C-V1S)/DVS                                                 
      NPTC = I2-I1+3                                                     
      V2C = V1C+DVS* REAL(NPTC-1)                                        
c
      DO 10 J = 1, NPTC                                                   F34310
         I = I1+J                                                         F34320
         C(J) = 0.                                                        F34330
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                            F34340
         vj = v1c + dvc* REAL(j-1)
c
         C(J) = factor*S(I)/vj  

   10 CONTINUE                                                            F34360
C                                                                         F34370
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
      SUBROUTINE O2HERZ (V1C,V2C,DVC,NPTC,C,T,P)                          F42270
C                                                                         F42280
      IMPLICIT REAL*8           (V)                                     ! F42290
C                                                                         F42300
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)                F42310
      DIMENSION C(*)                                                      F42320
C                                                                         F42330
      V1S = 36000.                                                        F42340
      DVS = 10.                                                           F42350
      DVC = DVS                                                           F42360
C                                                                         F42370
      V1C = V1ABS-DVC                                                     F42380
      V2C = V2ABS+DVC                                                     F42390
C                                                                         F42400
      I1 = (V1C-V1S)/DVS                                                  F42410
      IF (V1C.LT.V1S) I1 = I1-1                                           F42420
C                                                                         F42430
      V1C = V1S+DVS* REAL(I1)                                             F42440
      I2 = (V2C-V1S)/DVS                                                  F42450
      NPTC = I2-I1+3                                                      F42460
      V2C = V1C+DVS* REAL(NPTC-1)                                         F42470
      DO 10 J = 1, NPTC                                                   F42480
         I = I1+J                                                         F42490
         C(J) = 0.                                                        F42500
         IF (I.LT.1) GO TO 10                                             F42510
         VJ = V1C+DVC* REAL(J-1)                                          F42520
         CALL HERTDA (HERZ,VJ)                                            F42530
         CALL HERPRS (HERZ,T,P)                                           F42540
         C(J) = HERZ/VJ                                                   F42550
   10 CONTINUE                                                            F42560
C                                                                         F42570
      RETURN                                                              F42580
C                                                                         F42590
      END                                                                 F42600
C
C     --------------------------------------------------------------
C
      SUBROUTINE HERTDA (HERZ,V)                                          F42610
C                                                                         F42620
      IMPLICIT REAL*8           (V)                                     ! F42630
C                                                                         F42640
C     HERZBERG O2 ABSORPTION                                              F42650
C     HALL,1987 PRIVATE COMMUNICATION, BASED ON:                          F42660
C                                                                         F42670
C     REF. JOHNSTON, ET AL., JGR,89,11661-11665,1984                      F42680
C          NICOLET, 1987 (RECENT STUDIES IN ATOMIC                        F42690
C                         & MOLECULAR PROCESSES,                          F42700
C                         PLENUM PUBLISHING CORP, NY 1987)                F42710
C                                                                         F42720
C     AND YOSHINO, ET AL., 1988 (PREPRINT OF "IMPROVED ABSORPTION         F42730
C          CROSS SECTIONS OF OXYGEN IN THE WAVELENGTH REGION 205-240NM    F42740
C          OF THE HERZBERG CONTINUUM")                                    F42750
C                                                                         F42760
C     **** NOTE:  CROSS SECTION AT 0 PRESSURE  ***                        F42770
C     THE PRESSURE DEPENDENT TERM IS IN SUBROUTINE HERPRS                 F42780
C                                                                         F42790
CC    COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                       F42800
C                                                                         F42810
      HERZ = 0.0                                                          F42820
      IF (V.LE.36000.00) RETURN                                           F42830
C                                                                         F42840
C     EXTRAPOLATE SMOOTHLY THROUGH THE HERZBERG BAND REGION               F42850
C     NOTE: HERZBERG BANDS ARE NOT CORRECTLY INCLUDED                     F42860
C                                                                         F42870
      CORR = 0.                                                           F42880
      IF (V.LE.40000.) CORR = ((40000.-V)/4000.)*7.917E-27                F42890
C                                                                         F42900
C     UNITS ARE (CM2)                                                     F42910
C                                                                         F42920
C     HALL'S NEW HERZBERG  (LEAST SQRS FIT, LN(P))                        F42930
C                                                                         F42940
C     YRATIO=2048.7/WL(I)  ****IN ANGSTOMS****                            F42950
C           =.20487/WN(I)     IN MICRONS                                  F42960
C           =WCM(I)/48811.0   IN CM-1                                     F42970
C                                                                         F42980
      YRATIO = V/48811.0                                                  F42990
      HERZ = 6.884E-24*(YRATIO)*EXP(-69.738*( LOG(YRATIO))**2)-CORR       F43000
C                                                                         F43010
      RETURN                                                              F43020
C                                                                         F43030
      END                                                                 F43040
C
C     --------------------------------------------------------------
C
      SUBROUTINE HERPRS (HERZ,T,P)                                        F43050
C                                                                         F43060
C     CORRECT THE HERZBERG CONTINUUM CROSS SECTION FOR PRESSURE           F43070
C     DEPENDENCE; BASED ON SHARDANAND, JQRST, 18, 525-530, 1977.          F43080
C                 FOR UN2| BROADENING                                     F43090
C                 AND YOSHINO ET AL 1988 FOR UO2| BROADENING              F43100
C                                                                         F43110
C     PO2= PARTIAL PRESSURE OF O2                                         F43120
C     PN2= PARTIAL PRESSURE OF N2; BN2=.45*BO2                            F43130
C                                                                         F43140
C     DATA BO2 / 1.72E-3 /                                                F43150
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
      DATA PO / 1013. /,TO / 273.16 /                                     F43160
C                                                                         F43170
C     NOTE:  THE HERZBERG CONTINUUM OBEYS BEER'S LAW                      F43180
C            OPTICAL DEPTH(TOTAL)=SUM OVER LAYER O.D.(I)                  F43190
C                                                                         F43200
C     BO2= RATIO OF SIGMA(O2-O2)/(SIGMA(O2)) * 760(TORR)*.2095            F43210
C     BN2=.45*BO2= RATIO OF SIGMA(O2-N2)/(SIGMA(O2)) * 760(TORR)*.78      F43220
C                                                                         F43230
C     BO2*760*(.2095+.45*.78) = .73 , AS BELOW                            F43240
C
C     Changed Herzberg continuum pressure (see above reference)
C
C     BO2*760*(.2095+.45*.78) = .83 , AS BELOW
C                                                                         F43250
C                                                                         F43260
      HERZ = HERZ*(1.+.83*(P/PO)*(TO/T))                                  F43270
C                                                                         F43280
      RETURN                                                              F43290
C                                                                         F43300
      END                                                                 F43310


