c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is CKDv2.4, obtained from LBLRTMv5.10
c http://metosrv2.umd.edu/~bobe/LBLRTM/
c modified to make it kCARTA compatible by Sergio De Souza-Machado

c modified so it can run with the Matlab LBL code
c also modified it so that it can be run in a way as to allow Jacobians done ie
c immediately spline interpolate onto finer kCARTA grid from SL296,SL260 grids

c************************************************************************
       include 'dspline.f'
c************************************************************************
       SUBROUTINE CALCON24_loc(raCON, IDGAS,NFREQ,raFreq,FSTEP,NLAY,
     $    raT, raP, raPartP, raAMNT, selfmult,formult, whichlayer)

      IMPLICIT NONE

      include '../FORTRANFILES/max.inc'

C COMMON BLOCKS CONTAINING CONTINUUM VALUES 


C Arguements
       INTEGER IDGAS, NFREQ, NLAY,whichlayer
       REAL*8 raFREQ(*), FSTEP, raT(*), raP(*),raPARTP(*),raAMNT(*), 
     $        raCON(*)
       REAL*8 selfmult, formult 

C Variables for new water continuum adjustment terms
       REAL*8 SFAC,fscal

c general terms used in self, foreign
      REAL*8 vs2,vf2,vf6,vf4
      REAL*8 v0f3,hwsqf3

C COMMON BLOCKS CONTAINING CONTINUUM VALUES 

      INTEGER NPTS1,NPTS2,NPTF
      REAL*8 V1S1,V2S1,DVS1
      REAL*8 V1S2,V2S2,DVS2
      REAL*8 V1F,V2F,DVF
      REAL*8 raSSFREQ(2003)
      REAL*8 H2OS96(2003),H2OS60(2003),H2OF(2003)

c do not really need these for water continuum
      INTEGER NPTFCO2,NPN2RT0,NPN2RT1
      INTEGER NPN2_F,NPO2_F
      INTEGER NPTO3CHAP,NPTO3HH0,NPTO3HH1,NPTO3HH2,NPTO3HUV
      REAL*8 V1FCO2,V2FCO2,DVFCO2
      REAL*8 V1N2RT0,V2N2RT0,DVFN2RT0
      REAL*8 V1N2RT1,V2N2RT1,DVFN2RT1
      REAL*8 V1N2_F,V2N2_F,DVFN2_F
      REAL*8 V1O2_F,V2O2_F,DVFO2_F
      REAL*8 V1O3CHAP,V2O3CHAP,DVO3CHAP
      REAL*8 V1O3HH0,V2O3HH0,DVO3HH0
      REAL*8 V1O3HH1,V2O3HH1,DVO3HH1
      REAL*8 V1O3HH2,V2O3HH2,DVO3HH2
      REAL*8 V1O3HUV,V2O3HUV,DVO3HUV
      REAL*8 CO2F(1003)
      REAL*8 N296(73),N220(73),xn(118),xnt(118)
      REAL*8 X(3150),Y(3150),Z(3150)
      REAL*8 o3HH0(2687),o3HH1(2687),o3HH2(2687),o3HUV(133),o2F(206)

C Local variables
       REAL*8 raXAMT(kMaxLayer),raC11(kMaxLayer),
     $    raTFAC(kMaxLayer), A1, A2, A3

C Variables for new water continuum adjustment terms
       INTEGER JFAC
       REAL*8 XFAC(0:50)

      INTEGER iL

c ---------------------------------------------------------------

       REAL*8 V0S1, V0S2, V0S3, HWSQ1, HWSQ2, HWSQ3, 
     $    BETAS1, BETAS3, 
     $    FACTRS1, FACTRS2, FACTRS3,FACTRF3

       REAL*8 VJ, rSH2O, VS4, V0F1, HWSQF1, BETAF1, 
     $    V0F1a, HWSQF1a, BETAF1a,  FACTRF1a, BETAF2, 
     $    FACTRF2, BETAF3, FACTRF, FACTRF1
      
      INTEGER npts_lo, npts_hi

      INTEGER nptabs
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)
      REAL*8 V1C,V2C,DVC
      REAL*8 raContSelf(MaxLen),raContFor(MaxLen)
      REAL*8 raSH2OT0(MaxLen),raSH2OT1(MaxLen),raFH2O(MaxLen)
      INTEGER NPTC
                                                   
      REAL*8 XSELF,XFRGN
      REAL*8 V0F2,HWSQF2 
 
      REAL*8 SH2OT0(5050),SH2OT1(5050),FH2O(5050)
cc     *          CN2T0(5050),FCO2(5050),CT1(5050),CT2(5050) 
cc      REAL*8 CCH0(5150),CCH1(5150),CCH2(5150)

cc      REAL*8 ABSBSV(5050), XLOSMT
                                                             
cc      EQUIVALENCE (C0,SH2OT0,CN2T0,FCO2) , (C1,SH2OT1,CT1),             
cc     *            (C2,FH2O,CT2)                                         
             
      REAL*8 P0,T0                                                

      REAL*8 SFACL,SFACH,FL

      REAL*8 V1
      INTEGER iI,iJ

      REAL*8 kPlanck2,TS,AVOG

      COMMON /SH2O/ V1S1,V2S1,DVS1,NPTS1,H2OS96
      COMMON /S260/ V1S2,V2S2,DVS2,NPTS2,H2OS60
      COMMON /FH2O/ V1F,V2F,DVF,NPTF,H2OF
      COMMON /FCO2/ V1FCO2,V2FCO2,DVFCO2,NPTFCO2,CO2F
      COMMON /N2RT0/V1N2RT0,V2N2RT0,DVFN2RT0,NPN2RT0,N296
      COMMON /N2RT1/V1N2RT1,V2N2RT1,DVFN2RT1,NPN2RT1,N220
      COMMON /n2_f/ V1N2_F,V2N2_F,DVFN2_F,NPN2_F,xn,xnt
      COMMON /O3CHAP/ V1O3CHAP,V2O3CHAP,DVO3CHAP,NPTO3CHAP,X,Y,Z
      COMMON /O3HH0/ V1O3HH0,V2O3HH0,DVO3HH0,NPTO3HH0,o3HH0
      COMMON /O3HH1/ V1O3HH1,V2O3HH1,DVO3HH1,NPTO3HH1,o3HH1
      COMMON /O3HH2/ V1O3HH2,V2O3HH2,DVO3HH2,NPTO3HH2,o3HH2
      COMMON /O3HUV/ V1O3HUV,V2O3HUV,DVO3HUV,NPTO3HUV,o3HUV
      COMMON /o2_f  /V1O2_F,V2O2_F,DVFO2_F,NPO2_F,o2F 

      COMMON /ABSORB/V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      DATA kPlanck2/1.4387863/
      DATA AVOG/6.022045E+26/
      DATA TS/296.0/

      DATA P0 / 1013. /,T0 / 296. /                                     
cc      DATA XLOSMT / 2.68675E+19 /                                       


c     These are self-continuum modification factors from 700-1200 cm-1
      DATA (XFAC(iI),iI=0,50)/
     1    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
     2    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
     3    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
     4    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
     5    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
     6    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
     7    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
     8    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
     9    1.00000,1.00000,1.00000/
                                                             
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


C=======================================================================

C               ********    WATER VAPOR   ********                      

C=======================================================================
C                             SELF

       DVABS = 1.          
       DVABS = FSTEP
       V1ABS = INT(raFreq(1))     
       V1=raFreq(1)
       IF (V1.LT.0.) V1ABS = V1ABS-1. 
       V1ABS = V1ABS-3.*DVABS         
       V2ABS = INT(raFreq(NFREQ)+3.*DVABS+0.5)   
       NPTABS = (V2ABS-V1ABS)/DVABS+1.5 
       XSELF=1.0                                   
       XFRGN=1.0                                   

      !this is to put it into Genln2/kCARTA units
      DO iL=whichlayer,whichlayer                  !loop over layer!!!!
        raTFAC(iL) = (raT(iL)-T0)/(260.-T0)              
        raXamt(iL)=raAmnt(iL)*avog
        raC11(iL)=0.5*kPlanck2/raT(iL)
        END DO

      IF ((raFreq(NFREQ).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN
         ! for self, foreign broadening the v1c,v2c,dvc,nptc parameters
         ! are all the same (sl296, sl260,frn296)
         CALL SL296(raSSFREQ,V1C,V2C,DVC,NPTC,SH2OT0)           
         CALL SL260(raSSFREQ,V1C,V2C,DVC,NPTC,SH2OT1)   
         !now spline them onto the finer grid
         CALL XSPL(raSSFreq,SH2OT0,nptc,raFreq,raSH2OT0,nfreq)
         CALL XSPL(raSSFreq,SH2OT1,nptc,raFreq,raSH2OT1,nfreq)
         END IF

C     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
      IF ((raFreq(NFREQ).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN

        DO iL=whichlayer,whichlayer                     !loop over layers

C--------------------------------------------------------------------
C     *****    Continuum Correction Patches    ********
C--------------------------------------------------------------------
c                   SELF
          V0S1 = 0.
          HWSQ1 = 100.**2
          BETAS1 = 1.E-04
c         FACTRS1 = 0.3                      ! CKD2.2 value
          FACTRS1 = 0.688
 
          V0S2 = 1050.
          HWSQ2 = 200.**2
          FACTRS2 = -0.2333
 
          V0S3 = 1310.
          HWSQ3 = 120.**2
          BETAS3 = 5.E-06
          FACTRS3 = -0.15
C--------------------------------------------------------------------
c         Loop calculating self continuum optical depth
          DO 20 iJ = 1, NFREQ
            VJ = raFreq(iJ)
            rSH2O = 0.                            
            IF (raSH2OT0(iJ).GT.0.) THEN            
              rSH2O=raSH2OT0(iJ)*(raSH2OT1(iJ)/raSH2OT0(iJ))**raTFAC(iL) 
      
              SFAC = 1.
              IF (VJ.GE.700. .AND.  VJ.LE.1200.) THEN 
                JFAC = (VJ-700.)/10. + 0.00001
                SFAC = XFAC(JFAC)
ccc this is new, copied from CKD2.1 to make interpolation smoother across 10cm-1
                   JFAC=INT( (VJ - 700.0)/10.0 + 0.0001)
                   SFACL=XFAC(JFAC)
                   SFACH=XFAC(JFAC+1)
                   FL=700.0 + FLOAT(JFAC*10)
                   SFAC=SFACL + (VJ - FL)*(SFACH - SFACL)/10.0

                ENDIF
      
C       ---------------------------------------------------------
C         Correction to self continuum (1 SEPT 85); factor of    
C         0.78 at 1000 and  .......
                                  
              VS2 = (VJ-V0S1)**2
              VS4 = VS2*VS2
              SFAC=SFAC*(1.+FACTRS1*(HWSQ1/(VJ**2+(BETAS1*VS4)+HWSQ1)))  
 
              VS2 = (VJ-V0S2)**2
              SFAC = SFAC*(1.+FACTRS2*(HWSQ2/(VS2+HWSQ2)))
    
              VS2 = (VJ-V0S3)**2
              VS4 = VS2*VS2
              SFAC = SFAC*(1.+FACTRS3*(HWSQ3/(VS2+(BETAS3*VS4)+HWSQ3))) 
                                            
              rSH2O = SFAC * rSH2O
C    ---------------------------------------------------------
              ENDIF                 !  IF (SH2OT0(iJ).GT.0.) THEN            
 
c            raContSelf(iJ) = (rSH2O*RH2O)*XSELF
            raContSelf(iJ) = (raPartP(iL)*rSH2O)*XSELF
                                  
 20         CONTINUE                                           

C=======================================================================
c                     FOREIGN

          V0F1 = 350.
          HWSQF1 = 200.**2
          BETAF1 = 5.e-09 
          FACTRF1 = -0.7
 
          V0F1a = 630.
          HWSQF1a = 65.**2
          BETAF1a = 2.e-08 
          FACTRF1a = +0.75
 
          V0F2 =1130.
          HWSQF2 = 330.**2
          BETAF2 = 8.E-11
          FACTRF2 = -0.97
 
          V0F3 = 1975.
          HWSQF3 = 250.**2
          BETAF3 = 5.E-06
          FACTRF3 = -0.65
 
c        ------------------------------------------------------------

          CALL FRN296(raSSFREQ,V1C,V2C,DVC,NPTC,FH2O)         
          CALL XSPL(raSSFreq,FH2O,nptc,raFreq,raFH2O,nfreq)

          DO 24 iJ = 1, NFREQ
            VJ = raFreq(iJ)
   
C           CORRECTION TO FOREIGN CONTINUUM
   
            VF2 = (VJ-V0F1)**2
            VF6 = VF2 * VF2 * VF2
            FSCAL = (1.+FACTRF1*(HWSQF1/(VF2+(BETAF1*VF6)+HWSQF1)))
 
            VF2 = (VJ-V0F1a)**2
            VF6 = VF2 * VF2 * VF2
            FSCAL = FSCAL* 
     *              (1.+FACTRF1a*(HWSQF1a/(VF2+(BETAF1a*VF6)+HWSQF1a)))
 
            VF2 = (VJ-V0F2)**2
            VF6 = VF2 * VF2 * VF2
            FSCAL = FSCAL* 
     *              (1.+FACTRF2*(HWSQF2/(VF2+(BETAF2*VF6)+HWSQF2)))
 
            VF2 = (VJ-V0F3)**2
            VF4 = VF2*VF2
            FSCAL = FSCAL* 
     *              (1.+FACTRF3*(HWSQF3/(VF2+BETAF3*VF4+HWSQF3)))
      
            raFH2O(iJ)=raFH2O(iJ)*FSCAL
      
c            raContFor(iJ) = (raFH2O(iJ)*RFRGN)*XFRGN
            raContFor(iJ) = (raFH2O(iJ)*(raP(iL)-raPartP(iL)))*XFRGN
 
 24         CONTINUE                                                  
 
C           ------------------------------------------------------------
c       Interpolate to total optical depth grid
c       If spectral range covers the CKD_2.3 1 cm-1 H2O continuum
c       as well, then stop interpolation at 2200 cm-1.

          npts_lo = 1
          if (v1abs .lt. 2200.) then
            npts_lo =  npts_hi + 1
            v1c = v1c+npts_hi
            endif

C         The factor of 1.e-20 is handled this way to avoid underflows
          DO iJ=1,NFREQ
            VJ = raFreq(iJ)
            a1=VJ*raXamt(iL)*tanh(raC11(iL)*VJ)
            a2=TS/raT(iL)
            a3=1.0e-20*(raContFor(iJ)*formult + raContSelf(iJ)*selfmult)
            raCon(iJ)=a1*a2*a3
            END DO
          END DO               !end loop over layers

        ENDIF
         
      RETURN
      END
c************************************************************************

ccccccccc      include 'calconCKD2.4.data.f'

C=======================================================================

      SUBROUTINE SL296(raF,V1C,V2C,DVC,NPTC,C)                            
                                      
      REAL*8 V1C,V2C,DVC
      REAL*8 C(*),raF(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(2003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /SH2O/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                        
      V1C = V1ABS-DVC                                                  
      V2C = V2ABS+DVC                                                  
                                          
      I1 = (V1C-V1S)/DVS                                               
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC 
         raF(J)=v1c + (j-1)*dvc                                               
         I = I1+J                                                      
         C(J) = 0.                                                     
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                         
         C(J) = S(I)                                                   
   10 CONTINUE                                                         
                                          
      RETURN                                                           
                                          
      END                                                              
 
c************************************************************************
      SUBROUTINE SL260 (raF,V1C,V2C,DVC,NPTC,C)                              

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*),raF(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(2003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /S260/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC 
         raF(J)=v1c + (j-1)*dvc
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         C(J) = S(I)                                                     
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE FRN296 (raF,V1C,V2C,DVC,NPTC,C)                             
                                            
      REAL*8 V1C,V2C,DVC
      REAL*8 C(*),raF(*)                                       
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(2003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /FH2O/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC 
         raF(J)=v1c + (j-1)*dvc
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.GE.1).AND.(I.LE.NPTS)) THEN                              
            C(J) = S(I)                                                  
         ENDIF                                                           
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************

      SUBROUTINE FRNCO2 (V1C,V2C,DVC,NPTC,C)                             

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(1003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /FCO2/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.GE.1).AND.(I.LE.NPTS)) THEN                              
            C(J) = S(I)                                                  
         ENDIF                                                           
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE N2R296 (V1C,V2C,DVC,NPTC,C)
 
C     Model used:
C      Borysow, A, and L. Frommhold, "Collision-induced
C         rototranslational absorption spectra of N2-N2
C         pairs for temperatures from 50 to 300 K", The
C         Astrophysical Journal, 311, 1043-1057, 1986.
 
      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(73) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /N2RT0/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS
      V1C = V1ABS-DVC
      V2C = V2ABS+DVC
 
      I1 = (V1C-V1S)/DVS
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
 
      DO 10 J = 1, NPTC
         I = I1+J
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
   10 CONTINUE
 
      RETURN
 
      END
 
c************************************************************************
      SUBROUTINE N2R220 (V1C,V2C,DVC,NPTC,C)
C
C     Model used:
C      Borysow, A, and L. Frommhold, "Collision-induced
C         rototranslational absorption spectra of N2-N2
C         pairs for temperatures from 50 to 300 K", The
C         Astrophysical Journal, 311, 1043-1057, 1986.
C

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(73) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /N2RT1/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS
      V1C = V1ABS-DVC
      V2C = V2ABS+DVC
 
      I1 = (V1C-V1S)/DVS
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
 
      DO 10 J = 1, NPTC
         I = I1+J
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
   10 CONTINUE
 
      RETURN
 
      END
 
c************************************************************************
      subroutine n2_ver_1 (v1c,v2c,dvc,nptc,c,T)

      REAL*8 V1C,V2C,DVC,T
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,xn2(118),xn2t(118) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /N2_F/ V1S,V2S,DVS,NPTS,xn2,xn2t                         
      
      INTEGER I1,I2,I,J
      REAL*8 xktfac,a1,a2,factor,vj

c     Nitrogen Collision Induced Fundamental

c     Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M. Hartman
c        Infrared collision-induced absorption by N2 near 4.3 microns for
c        atmospheric applications: measurements and emprirical modeling, 
c         Appl. Optics, 35, 5911-5917, (1996).
c
      REAL*8 To,xlosmt,vmr_n2
      DATA  To/ 296./, xlosmt/ 2.68675e+19/, vmr_n2/ 0.78 /
 
      xktfac = (1./To)-(1./T)
      
      a1 = 0.8387
      a2 = 0.0754
 
c     correct formulation for consistency with LBLRTM:
 
      factor = (1.e+20 /xlosmt) * (1./vmr_n2) * (a1-a2*(T/To))
 
c     Lafferty et al. reference  assumes that the
c     column amount is that for air 

                            
      DVC = DVS             
      V1C = V1ABS-DVC       
      V2C = V2ABS+DVC       
                            
      I1 = (V1C-V1S)/DVS    
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      do 10 j=1,nptc
         i = i1+j
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         VJ = V1C+DVC*FLOAT(J-1)                                         
c     the radiation field is removed with 1/vj
 
         c(j) = factor * xn2(i)* exp(xn2t(i)*xktfac) / vj
 
 10   end do
 920  format (f10.2,1p,e12.2,0p,f10.2,1p2e12.2)
      return

      end
c************************************************************************
      SUBROUTINE XO3CHP (V1C,V2C,DVC,NPTC,C0,C1,C2)                      

      REAL*8 V1C,V2C,DVC
      REAL*8 C0(*),C1(*),C2(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,x(3150),y(3150),z(3150)
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /O3CHAP/ V1S,V2S,DVS,NPTS,x,y,z                        
      
      INTEGER I1,I2,I,J
      REAL*8 VJ

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         IF ((I.LT.1).OR.(I.GT.NPTS)) THEN
             C0(J) = 0.
             C1(J)=0.
             C2(J)=0.
         ELSE
 
C            Remove radiation field from diffuse ozone
 
             VJ = V1C+DVC*FLOAT(J-1)
             C0(J)=X(I)/VJ
             C1(J)=Y(I)/VJ
             C2(J)=Z(I)/VJ
         ENDIF
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE O3HHT0 (V1C,V2C,DVC,NPTC,C)                             

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(2687)
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /O3HH0/ V1S,V2S,DVS,NPTS,S
      
      INTEGER I1,I2,I,J
      REAL*8 VJ
                                            
      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         VJ = V1C+DVC*FLOAT(J-1)                                         
         C(J) = S(I)/VJ                                                  
                                            
C     RADIATION FLD REMOVED FROM DIFFUSE OZONE                           
                                            
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
C
c************************************************************************
      SUBROUTINE O3HHT1 (V1C,V2C,DVC,NPTC,C)                             

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(2687)
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /O3HH1/ V1S,V2S,DVS,NPTS,S
      
      INTEGER I1,I2,I,J
                                            
      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         C(J) = S(I)                                                     
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE O3HHT2 (V1C,V2C,DVC,NPTC,C)                             
                                            
      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(2687)
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /O3HH2/ V1S,V2S,DVS,NPTS,S
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         C(J) = S(I)                                                     
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE O3HHUV (V1C,V2C,DVC,NPTC,C)                             

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*)
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(133)
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /O3HUV/ V1S,V2S,DVS,NPTS,S
      
      INTEGER I1,I2,I,J
      REAL*8 VJ
                                            
      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         VJ = V1C+DVC*FLOAT(J-1)                                         
         C(J) = S(I)/VJ                                                  
                                            
C     RADIATION FLD REMOVED FROM U.V.    OZONE                           
                                            
   10 CONTINUE                                                           
                                            
      RETURN
      END
 
c************************************************************************
      subroutine o2_ver_1 (v1c,v2c,dvc,nptc,c,T)
 
      REAL*8 V1C,V2C,DVC,T
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,xo2(103),xo2t(103) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /O2_F/ V1S,V2S,DVS,NPTS,xo2,xo2t                         
      
      INTEGER I1,I2,I,J
      REAL*8 xktfac,factor,vj

c     Oxygen Collision Induced Fundamental

c     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann,
c                                                         and Ch. Boulet
c        Infrared collision-induced absorption by O2 near 6.4 microns for
c        atmospheric applications: measurements and emprirical modeling, 
c         Appl. Optics, 35, 5911-5917, (1996).

      REAL*8 To,xlosmt
      DATA To/ 296./, xlosmt/ 2.68675e+19/
 
      xktfac = (1./To)-(1./T)
      
c     correct formulation for consistency with LBLRTM:
 
      factor = (1.e+20 /xlosmt) 
 
c     A factor of 0.21, the mixing ration of oxygen, in Thibault etal.
c     formulation is not included here.  This is in the column amt.
                            
      DVC = DVS             
      V1C = V1ABS-DVC       
      V2C = V2ABS+DVC       
                            
      I1 = (V1C-V1S)/DVS        
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      do 10 j=1,nptc
         i = i1+j
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         VJ = V1C+DVC*FLOAT(J-1)                                         
c     the radiation field is removed with 1/vj
 
         c(j) = factor * xo2(i)* exp(xo2t(i)*xktfac) / vj
 
 10   end do
 
 920  format (f10.2,1p,e12.2,0p,f10.2,1p2e12.2)
 
      return
      end

c************************************************************************
      SUBROUTINE O2INF1 (V1C,V2C,DVC,NPTC,C,T,P)                         

      REAL*8 V1C,V2C,DVC,T,P
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      INTEGER NPTABS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      
      INTEGER I1,I2,I,J
                                            
      REAL*8 V1S,DVS,VJ,O2INF

      V1S = 7600.                                                        
      DVS = 1.                                                          
      DVC = DVS                                                          
                                            
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = I1-1                                          
                                            
      V1C = V1S+DVS*FLOAT(I1)                                            
      I2 = (V2C-V1S)/DVS                                                 
      NPTC = I2-I1+3                                                     
      V2C = V1C+DVS*FLOAT(NPTC-1)                                        
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF (I.LT.1) GO TO 10                                            
         VJ = V1C+DVC*FLOAT(J-1)                                         
         CALL INFRD1 (O2INF,VJ)                                          
         C(J) = O2INF/VJ                                                 
   10 CONTINUE                                                           

      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE O2INF2 (V1C,V2C,DVC,NPTC,C,T,P)                         

      REAL*8 V1C,V2C,DVC,T,P
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      INTEGER NPTABS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      
      INTEGER I1,I2,I,J
      REAL*8 V1S,DVS,VJ,O2INF

      V1S = 9040.                                                        
      DVS = 1.                                                          
      DVC = DVS                                                          
                                            
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = I1-1                                          
                                            
      V1C = V1S+DVS*FLOAT(I1)                                            
      I2 = (V2C-V1S)/DVS                                                 
      NPTC = I2-I1+3                                                     
      V2C = V1C+DVS*FLOAT(NPTC-1)                                        
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF (I.LT.1) GO TO 10                                            
         VJ = V1C+DVC*FLOAT(J-1)                                         
         CALL INFRD2 (O2INF,VJ)                                          
         C(J) = O2INF/VJ                                                 
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE INFRD1 (O2INF,V)                                        

      REAL*8 O2INF,V

      REAL*8 DV
                                            
      O2INF = 0.0                                                        
      IF (V.LE.7500.00) RETURN                                           
      IF (V.GE.8300.00) RETURN
 
      DV = V - 7896.464      
      O2INF = 1.054*3.159e-26/(1.+(DV/54.93)**2+(DV/75.63)**4)

      RETURN

 9000 CONTINUE
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE INFRD2 (O2INF,V)                                        

      REAL*8 O2INF,V

      REAL*8 DAMP1,DAMP2,DV1,DV2

      REAL*8 V1,HW1,V2,HW2,S1,S2                                            
      DATA V1 /9375./, HW1 /58.96/, V2 /9439./, HW2 /45.04/
      DATA S1 /1.166E-24/, S2 /3.086E-25/

      O2INF = 0.0                                                        
      IF (V .LE. 9100.00) RETURN                                         
      IF (V .GE. 11000.00) RETURN                                        

      DV1 = V - V1
      DV2 = V - V2
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
     &     + ((S2 * DAMP2 / HW2)/(1. + (DV2/HW2)**2))) * 1.054
 9000 CONTINUE
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE O2HERZ (V1C,V2C,DVC,NPTC,C,T,P)                         
                                            
      REAL*8 V1C,V2C,DVC,T,P
      REAL*8 C(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS

      REAL*8 V1S,DVS,VJ,HERZ
      INTEGER I1,J,I2,I

      V1S = 36000.                                                       
      DVS = 10.                                                          
      DVC = DVS                                                          
                                            
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = I1-1                                          
                                            
      V1C = V1S+DVS*FLOAT(I1)                                            
      I2 = (V2C-V1S)/DVS                                                 
      NPTC = I2-I1+3                                                     
      V2C = V1C+DVS*FLOAT(NPTC-1)                                        
      DO 10 J = 1, NPTC                                                  
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF (I.LT.1) GO TO 10                                            
         VJ = V1C+DVC*FLOAT(J-1)                                         
         CALL HERTDA (HERZ,VJ)                                           
         CALL HERPRS (HERZ,T,P)                                          
         C(J) = HERZ/VJ                                                  
   10 CONTINUE                                                           
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE HERTDA (HERZ,V)                                         
                                            
C     HERZBERG O2 ABSORPTION                                             
C     HALL,1987 PRIVATE COMMUNICATION, BASED ON:                         
                                            
C     REF. JOHNSTON, ET AL., JGR,89,11661-11665,1984                     
C          NICOLET, 1987 (RECENT STUDIES IN ATOMIC                       
C                         & MOLECULAR PROCESSES,                         
C                         PLENUM PUBLISHING CORP, NY 1987)               
                                            
C     AND YOSHINO, ET AL., 1988 (PREPRINT OF "IMPROVED ABSORPTION        
C          CROSS SECTIONS OF OXYGEN IN THE WAVELENGTH REGION 205-240NM   
C          OF THE HERZBERG CONTINUUM")                                   
                                            
C     **** NOTE:  CROSS SECTION AT 0 PRESSURE  ***                       
C     THE PRESSURE DEPENDENT TERM IS IN SUBROUTINE HERPRS                
                                            
C     COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      

      REAL*8 HERZ,V
      REAL*8 CORR,YRATIO

      HERZ = 0.0                                                         
      IF (V.LE.36000.00) RETURN                                          
                                            
C     EXTRAPOLATE SMOOTHLY THROUGH THE HERZBERG BAND REGION              
C     NOTE: HERZBERG BANDS ARE NOT CORRECTLY INCLUDED                    
                                            
      CORR = 0.                                                          
      IF (V.LE.40000.) CORR = ((40000.-V)/4000.)*7.917E-27               
                                            
C     UNITS ARE (CM2)                                                    
                                            
C     HALL'S NEW HERZBERG  (LEAST SQRS FIT, LN(P))                       
                                            
C     YRATIO=2048.7/WL(I)  ****IN ANGSTOMS****                           
C           =.20487/WN(I)     IN MICRONS                                 
C           =WCM(I)/48811.0   IN CM-1                                    
                                            
      YRATIO = V/48811.0                                                 
      HERZ = 6.884E-24*(YRATIO)*EXP(-69.738*(DLOG(YRATIO))**2)-CORR      
                                            
      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE HERPRS (HERZ,T,P)                                       
                                            
C     CORRECT THE HERZBERG CONTINUUM CROSS SECTION FOR PRESSURE          
C     DEPENDENCE; BASED ON SHARDANAND, JQRST, 18, 525-530, 1977.         
C                 FOR UN2| BROADENING                                    
C                 AND YOSHINO ET AL 1988 FOR UO2| BROADENING             
                                            
C     PO2= PARTIAL PRESSURE OF O2                                        
C     PN2= PARTIAL PRESSURE OF N2; BN2=.45*BO2                           
                                            
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

      REAL*8 HERZ,T,P
      REAL*8 BO2,PO,TO


      DATA BO2 / 1.81E-3 /
      DATA PO / 1013. /,TO / 273.16 /                                    
                                            
C     NOTE:  THE HERZBERG CONTINUUM OBEYS BEER'S LAW                     
C            OPTICAL DEPTH(TOTAL)=SUM OVER LAYER O.D.(I)                 
                                            
C     BO2= RATIO OF SIGMA(O2-O2)/(SIGMA(O2)) * 760(TORR)*.2095           
C     BN2=.45*BO2= RATIO OF SIGMA(O2-N2)/(SIGMA(O2)) * 760(TORR)*.78     
                                            
C     BO2*760*(.2095+.45*.78) = .73 , AS BELOW                           
C
C     Changed Herzberg continuum pressure (see above reference)
C
C     BO2*760*(.2095+.45*.78) = .83 , AS BELOW
                                            
                                            
      HERZ = HERZ*(1.+.83*(P/PO)*(TO/T))                                 
                                            
      RETURN                                                             
                                            
      END                                                                

c************************************************************************
