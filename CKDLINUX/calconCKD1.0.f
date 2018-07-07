c Copyright 2003 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is MTCKDv1.0 obtained from Dave Tobin
c modified to make it run7watercontinuum compatible
c by Sergio De Souza-Machado

c modified so it can run with the Matlab LBL code
c ALSO MODIFIED XFACREV SO IT STARTS AT 810 cm-1 instead of 820 cm-1
c ie IT IS 1.0 at 810, and 1.003 at 820 etc ....
C INSTEAD OF 1.003 at 820 etc
C so the arrays are (0:15) instead of (0:14) 

c     16 December 2002
c
c     This  version of the water vapor continuum, mt_ckd_1.00, is a completely new 
c     continuum model based on a collision induced  component and a sub-Lorentzian 
c     line wing component.  Continua realted to other species are the same as for 
c     ckd_2.4.2.
c
c     this is an updated version of the continuum program:
c     this version provides optical depths on file CNTNM.OPTDPT as before:
c     it also provides the continuum coefficients on file  WATER.COEF
c
c     the length of the header records may vary by version:
c         in this version the WATER.COEF header information is 47 records 
c         in this version the CNTNM.OPTDT header information is 34 records 
c
c     presumably the user will want to create an input file to address
c     individual requirements
C
C   THE FOLLOWING QUANTITIES MUST BE SPECIFIED: 
C
C          PRESSURE                   PAVE (MB)
C
C          TEMPERATURE                TAVE ( K)
C
C          COLUMN AMOUNT
C            NITROGEN                 WN2    (MOLEC/CM**2)
C            OXYGEN                   WK(7)  (MOLEC/CM**2)
C            CARBON DIOXIDE           WK(2)  (MOLEC/CM**2)
C            WATER VAPOR              WK(1)  (MOLEC/CM**2)
C
C          NUMBER OF MOLECULES        NMOL
C
C          BEGINNING WAVENUMBER       V1ABS (CM-1)
C
C          ENDING WAVENUMBER          V2ABS (CM-1)
C
C          SAMPLING INTERVAL          DVABS (CM-1)
C
C          NUMBER OF VALUES           NPTABS
C
c************************************************************************
       SUBROUTINE CALCON_MTCKD_01_loc(raCON, IDGAS,NFREQ,raFreq,FSTEP,
     $    NLAY, raT, raP, raPartP, raAMNT, selfmult,formult, whichlayer)

      IMPLICIT NONE

      include '../FORTRANFILES/max.inc'

C COMMON BLOCKS CONTAINING CONTINUUM VALUES 


C Arguments
       INTEGER IDGAS, NFREQ, NLAY,whichlayer
       REAL*8 raFREQ(*), FSTEP, raT(*), raP(*),raPARTP(*),raAMNT(*), 
     $        raCON(*)
       REAL*8 selfmult, formult 

C Variables for new water continuum adjustment terms
       REAL*8 SFAC,fscal 

c these are for common blocks
      INTEGER NPTS1,NPTS2,NPTF
      REAL*8 V1S1,V2S1,DVS1
      REAL*8 V1S2,V2S2,DVS2
      REAL*8 V1F,V2F,DVF

      REAL*8 raSSFREQ(2003)
      REAL*8 H2OS96(2003),H2OS60(2003),H2OF(2003)

C Local variables
       REAL*8 raXAMT(kMaxLayer),raC11(kMaxLayer),raTFAC(kMaxLayer)

C Variables for new water continuum adjustment terms
       INTEGER JFAC

      INTEGER iL

c ---------------------------------------------------------------

       REAL*8 V0S1, V0S2, V0S3, HWSQ1, HWSQ2, HWSQ3, 
     $    BETAS1, BETAS3, 
     $    FACTRS1, FACTRS2, FACTRS3

       REAL*8 VJ, rSH2O, VS4, V0F1, HWSQF1, BETAF1, 
     $    V0F1a, HWSQF1a, BETAF1a,  FACTRF1a, BETAF2, 
     $    FACTRF2, BETAF3, FACTRF3, FACTRF1

      REAL*8 XFACREV(1:16) 
      REAL*8 FREQFACREV(1:16) 
      REAL*8 raXFACREV(MaxLen)
      
      INTEGER npts_lo, npts_hi

      INTEGER NPTABS
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)
      REAL*8 V1C,V2C,DVC
      REAL*8 raContSelf(MaxLen),raContFor(MaxLen)
      REAL*8 raSH2OT0(MaxLen),raSH2OT1(MaxLen),raFH2O(MaxLen)
      INTEGER NPTC

      REAL*8 XSELF,XFRGN,PFRGN,RFRGN
      REAL*8 V0F2,HWSQF2
 
      REAL*8 SH2OT0(5050),SH2OT1(5050),FH2O(5050)

      REAL*8 V1,V2
      INTEGER iI,iJ

      REAL*8 P0,T0                                                
      REAL*8 kPlanck2,TS,AVOG

      REAL*8 a1,a2,a3

      REAL rDummy

      COMMON /ABSORB_MTCKD1/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB 
      COMMON /SH2O_MTCKD1/ V1S1,V2S1,DVS1,NPTS1,H2OS96
      COMMON /S260_MTCKD1/ V1S2,V2S2,DVS2,NPTS2,H2OS60
      COMMON /FH2O_MTCKD1/ V1F,V2F,DVF,NPTF,H2OF

      DATA P0 / 1013.25 /,T0 / 296. /
      DATA kPlanck2/1.4387863/
      DATA AVOG/6.022045E+26/
      DATA TS/296.0/

      DATA (XFACREV(iI),iI=1,16)/ 
     1     1.000,1.003,1.009,1.015,1.023,1.029,1.033, 
     2     1.037,1.039,1.040,1.046,1.036,1.027, 
     3     1.01,1.002,1.00/ 

      DO iI = 1,16
        FREQFACREV(iI) = 810.0 + (iI-1)*10.0
        END DO

      DVABS = 1.           
      DVABS = FSTEP 
      V1ABS = INT(raFreq(1))      
      V1 = raFreq(1) 
      V2 = raFreq(nfreq)
      IF (V1.LT.0.) V1ABS = V1ABS-1.  
      V1ABS = V1ABS-3.*DVABS          
      V2ABS = INT(raFreq(NFREQ)+3.*DVABS+0.5)    
      NPTABS = (V2ABS-V1ABS)/DVABS+1.5  
      XSELF=1.0                                    
      XFRGN=1.0                                    

      !this is to put it into Genln2/kCARTA units 
      iL=whichlayer                  !loop over layer!!!! 

      raTFAC(iL) = (raT(iL)-T0)/(260.-T0)               
      raXamt(iL)=raAmnt(iL)*avog 
      raC11(iL)=0.5*kPlanck2/raT(iL) 
 
C=======================================================================

C               ********    WATER VAPOR   ********
C                                                 
C=======================================================================
C                             SELF

C        Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
c
      IF ((raFreq(NFREQ).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN 
        CALL SL296_MTCKD1 (raSSFREQ,V1C,V2C,DVC,NPTC,SH2OT0)
        CALL SL260_MTCKD1 (raSSFREQ,V1C,V2C,DVC,NPTC,SH2OT1)
        !now spline them onto the finer grid 
        CALL XSPL(raSSFreq,SH2OT0,nptc,raFreq,raSH2OT0,nfreq) 
        CALL XSPL(raSSFreq,SH2OT1,nptc,raFreq,raSH2OT1,nfreq) 
c new!!!!!!!!
        CALL XSPL(FREQFACREV,XFACREV,16,raFreq,raXFACREV,nfreq) 

        DO iJ = 1, NFREQ
          VJ = raFreq(iJ)
          rSH2O = 0.    
          IF (raSH2OT0(iJ).GT.0.) THEN
            rSH2O = raSH2OT0(iJ)*
     $                (raSH2OT1(iJ)/raSH2OT0(iJ))**raTFAC(iL)
            SFAC = 1.
            IF (VJ .GE. 810. .AND. VJ .LE. 960.) THEN
c              JFAC = (VJ - 810.)/10. + 0.00001
c              SFAC = XFACREV(JFAC)
              SFAC = raXFACREV(iJ)
c              print *,vj,sfac
              ENDIF
            rSH2O = SFAC * rSH2O
            ENDIF
c *************** we really need csh2o(i) from cntnm_progr.f
          raContSelf(iJ) = (raPartP(iL)*rSH2O)*XSELF 
          END DO
        ENDIF

C=======================================================================
C                             FOREIGN
C
c      PFRGN = PATM-PH2O
c      RFRGN = PFRGN*(T0/TAVE)
       PFRGN = raP(iL)-raPartP(iL)
       RFRGN = PFRGN*(T0/raT(iL))

C     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
      IF ((raFreq(NFREQ).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN 
        CALL FRN296_MTCKD1 (raSSFREQ,V1C,V2C,DVC,NPTC,FH2O)         
        CALL XSPL(raSSFreq,FH2O,nptc,raFreq,raFH2O,nfreq) 
        DO iJ = 1, NFREQ
           VJ = raFreq(iJ)                                         
c          C(J) = W1*(FH2O(J)*RFRGN)*XFRG
c           raContFor(iJ) = (raFH2O(iJ)*(raP(iL)-raPartP(iL)))*XFRGN 
c           raContFor(iJ) = raFH2O(iJ)*RFRGN*XFRGN 
c *************** we really need cfh2o(i) from cntnm_progr.f
c but then we need to multiply by PFRGN = raP(iL)-raPartP(iL)
c           raContFor(iJ) = raFH2O(iJ)*XFRGN 
           raContFor(iJ) = raFH2O(iJ)*PFRGN*XFRGN 
           END DO
         ENDIF          

C        Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
      IF ((raFreq(NFREQ).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN 
C       The factor of 1.e-20 is handled this way to avoid underflows
        DO iJ=1,NFREQ
          VJ = raFreq(iJ)
          a1=VJ*raXamt(iL)*tanh(raC11(iL)*VJ)
          a2=TS/raT(iL)
          a3=1.0e-20*(raContFor(iJ)*formult + raContSelf(iJ)*selfmult)
          raCon(iJ)=a1*a2*a3
c          IF (iJ .LE. 35) THEN
c            print *,'    '
c            print *,iJ,raFreq(ij),a1*a2*a3
c            print *,a1,a2,a3
c            print *,raContSelf(iJ),raContFor(iJ)       
c            print *,xself,xfrgn,VJ,raXamt(iL),raTFAC(iL)
c            print *,raP(iL),raPartP(iL),raT(iL)
c            print *,raSH2OT0(iJ),raSH2OT1(iJ),raFH2O(iJ)
c            END IF
          END DO
        END IF               !end loop over layers
     
       RETURN
       END
C=======================================================================
C
      SUBROUTINE SL296_MTCKD1 (raF,V1C,V2C,DVC,NPTC,C)
C             
      IMPLICIT NONE

      REAL*8 V1C,V2C,DVC 
      REAL*8 C(*),raF(*)                                                    
      INTEGER NPTC 

      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)                   
      REAL*8 V1S,V2S,DVS,S(2003)  
      INTEGER NPTABS,NPTS 
 
      COMMON /ABSORB_MTCKD1/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB 
      COMMON /SH2O_MTCKD1/ V1S,V2S,DVS,NPTS,S

      INTEGER I1,I2,I,J 

      DVC = DVS
      V1C = V1ABS-DVC
      V2C = V2ABS+DVC

      I1 = (V1C-V1S)/DVS
      IF (V1C.LT.V1S) I1 = -1 

      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       

      DO 10 J = 1, NPTC
         raF(J)=v1c + (j-1)*dvc
 72      I = I1+J
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
   10 CONTINUE

      RETURN

      END

C     --------------------------------------------------------------

      SUBROUTINE SL260_MTCKD1 (raF, V1C,V2C,DVC,NPTC,C)

      IMPLICIT NONE

      REAL*8 V1C,V2C,DVC 
      REAL*8 C(*),raF(*)                                                    
      INTEGER NPTC 
     
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)                   
      REAL*8 V1S,V2S,DVS,S(2003)  
      INTEGER NPTABS,NPTS 
 
      COMMON /ABSORB_MTCKD1/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB 
      COMMON /S260_MTCKD1/ V1S,V2S,DVS,NPTS,S

      INTEGER I1,I2,I,J 

      DVC = DVS
      V1C = V1ABS-DVC
      V2C = V2ABS+DVC

      I1 = (V1C-V1S)/DVS
      IF (V1C.LT.V1S) I1 = -1 

      V1C = V1S+DVS* REAL(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       

      DO 10 J = 1, NPTC
         raF(J)=v1c + (j-1)*dvc
         I = I1+J
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
   10 CONTINUE

      RETURN
      END

C     --------------------------------------------------------------

      SUBROUTINE FRN296_MTCKD1 (raF,V1C,V2C,DVC,NPTC,C)

      IMPLICIT NONE

      REAL*8 V1C,V2C,DVC 
      REAL*8 C(*),raF(*)                                                    
      INTEGER NPTC 
     
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)                   
      REAL*8 V1S,V2S,DVS,S(2003)  
      INTEGER NPTABS,NPTS 
 
      COMMON /ABSORB_MTCKD1/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB 
      COMMON /FH2O_MTCKD1/ V1S,V2S,DVS,NPTS,S
 
      INTEGER I1,I2,I,J 

      DVC = DVS                          
      V1C = V1ABS-DVC                   
      V2C = V2ABS+DVC                  

      I1 = (V1C-V1S)/DVS             
      IF (V1C.LT.V1S) I1 = -1

      V1C = V1S+DVS* REAL(I1)
      I2 = (V2C-V1S)/DVS    
      NPTC = I2-I1+3       
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS* REAL(NPTC-1)       

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

C     --------------------------------------------------------------

