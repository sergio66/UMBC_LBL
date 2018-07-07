c SAME as calconCKD4.0.f except it applies the Dave Tobin PhD thesis 
c tuning multiplier both to SELF AND FOREIGN in the 1300-1600 cm-1 region,
c as opposed to CKD3 which applies the tuning multiplier
c only to the SELF

c so the SELF02/SELF04 tables should be the same while the FOR02/FOR04 
c should differ

c Copyright 2003 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is MTCKDv1.0 obtained from Dave Tobin
c modified to make it run7watercontinuum compatible
c MODIFIED BY SCOTT HANNON SO THAT IT AGREES WITH ECMWF  Aug 21, 2003
c v2a does this from /asl/packages/sartaV104/Src_tuning/tunmlt_con1_14jul03.txt
c v2  does this from /asl/packages/sartaV104/Src_tuning/tunmlt_31oct03.txt
c v3  does this from /asl/packages/sartaV104/Src_tuning/tunmlt_jan04deliv.txt
c v4  does this from /asl/packages/sartaV104/Src_tuning/tunmlt_jan04deliv.txt
c look at columns 2,5

c by Sergio De Souza-Machado
c************************************************************************
       SUBROUTINE CALCON_MTCKD_05_loc(raCON, IDGAS,NFREQ,raFreq,FSTEP,
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

c variables to read in scott's fixes
       INTEGER iLenX,iLenScottFile
       PARAMETER (iLenX = 10000)
       REAL*8 aa1,raAIRS_WN(iLenX),aa3,aa4,raAIRS_FIX(iLenX),
     $        aa6,aa7,aa8,aa9
       INTEGER iIOUN,iErr
       CHARACTER*80 FNAME
       REAL*8 raScottMultiplier(MaxLen)

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

c------------------------------------------------------------------------
 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80) 

      iLenScottFile = 2378
      FNAME = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_con1_14jul03.dat'
      FNAME = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_31oct03.dat'
      FNAME = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.dat'

      iIOUN = 11
      OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $     IOSTAT=IERR) 
      IF (IERR .NE. 0) THEN 
        write(*,1010) IERR, FNAME 
        STOP 
        ENDIF 
      DO iI = 1,iLenScottFile
        READ(iIOUN,*) aa1,raAIRS_WN(iI),aa3,aa4,raAIRS_FIX(iI),
     $                aa6,aa7,aa8,aa9
        END DO
      CLOSE(iIOUN)

      CALL xlinear(raAIRS_WN,raAIRS_FIX,iLenScottFile,raFreq,
     $             raScottMultiplier,NFREQ)
      DO iI=1,NFREQ
        IF (raFreq(iI) .LT. raAIRS_WN(1)) THEN
          raScottMultiplier(iI) = 1.0d0
        ELSEIF (raFreq(iI) .GT. raAIRS_WN(iLenScottFile)) THEN
          raScottMultiplier(iI) = 1.0d0
          END IF
        END DO
c------------------------------------------------------------------------
c this file comes from running use_tobin_continuum.m
c which I gave up doing and intead just made /home/sergio/SPECTRA/tobin5.m

      print *,'idiot : use /home/sergio/SPECTRA/tobin5.m'
      stop 

 1020 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
      FNAME = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_davetobin.dat'
      iIOUN = 11
      OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $     IOSTAT=IERR) 
      IF (IERR .NE. 0) THEN 
        write(*,1010) IERR, FNAME 
        STOP 
        ENDIF 
      DO iI = 1,iLenScottFile
        READ(iIOUN,*) aa1,raAIRS_WN(iI),aa3,aa4,raAIRS_FIX(iI),
     $                aa6,aa7,aa8,aa9
        END DO
      CLOSE(iIOUN)

      CALL xlinear(raAIRS_WN,raAIRS_FIX,iLenScottFile,raFreq,
     $             raScottMultiplier,NFREQ)
      DO iI=1,NFREQ
        IF (raFreq(iI) .LT. raAIRS_WN(1)) THEN
          raScottMultiplier(iI) = 1.0d0
        ELSEIF (raFreq(iI) .GT. raAIRS_WN(iLenScottFile)) THEN
          raScottMultiplier(iI) = 1.0d0
          END IF
        END DO
c------------------------------------------------------------------------

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
c ************** modified by Scott's multipliers
          raContSelf(iJ) = (raPartP(iL)*rSH2O)*XSELF * 
     $                     raScottMultiplier(iJ)
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
ccccc for CKD2           raContFor(iJ) = raFH2O(iJ)*PFRGN*XFRGN 
           raContFor(iJ) = raFH2O(iJ)*PFRGN*XFRGN*raScottMultiplier(iJ)
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
