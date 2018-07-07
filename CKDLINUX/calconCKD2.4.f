c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is CKDv2.4, obtained from LBLRTMv5.10
c http://metosrv2.umd.edu/~bobe/LBLRTM/
c modified to make it kCARTA compatible by Sergio De Souza-Machado

c modified so it can run with the Matlab LBL code
c also modified it so that it can be run in a way as to allow Jacobians done ie
c immediately spline interpolate onto finer kCARTA grid from SL296,SL260 grids

c THIS WAS USED UPTIL DEC 2002

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
     $    FACTRS1, FACTRS2, FACTRS3

       REAL*8 VJ, rSH2O, VS4, V0F1, HWSQF1, BETAF1, 
     $    V0F1a, HWSQF1a, BETAF1a,  FACTRF1a, BETAF2, 
     $    FACTRF2, BETAF3, FACTRF3, FACTRF1
      
      INTEGER npts_lo, npts_hi

      INTEGER NPTABS
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
             
      REAL*8 SFACL,SFACH,FL

      REAL*8 V1
      INTEGER iI,iJ

      REAL*8 P0,T0                                                
      REAL*8 kPlanck2,TS,AVOG

      COMMON /SH2O/ V1S1,V2S1,DVS1,NPTS1,H2OS96
      COMMON /S260/ V1S2,V2S2,DVS2,NPTS2,H2OS60
      COMMON /FH2O/ V1F,V2F,DVF,NPTF,H2OF

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

c      DO iJ = 1,10
c        print *,iJ,raSSFREQ(ij),sh2ot0(ij),sh2ot1(ij),fh2o(ij)
c        END DO


C         The factor of 1.e-20 is handled this way to avoid underflows
          DO iJ=1,NFREQ
            VJ = raFreq(iJ)
            a1=VJ*raXamt(iL)*tanh(raC11(iL)*VJ)
            a2=TS/raT(iL)
            a3=1.0e-20*(raContFor(iJ)*formult + raContSelf(iJ)*selfmult)
            raCon(iJ)=a1*a2*a3
c            IF (iJ .EQ. 1) THEN
c              print *,a1,a2,a3
c              print *,raContSelf(iJ),raContFor(iJ)       
c              print *,xself,xfrgn,VJ,raXamt(iL),raTFAC(iL)
c              print *,raP(iL),raPartP(iL),raT(iL)
c              print *,raSH2OT0(iJ),raSH2OT1(iJ),raFH2O(iJ)
c              END IF
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
