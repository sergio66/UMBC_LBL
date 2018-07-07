c Copyright 2003 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is MTCKDv1.3 obtained from latest LBLRTM version (Dec 2006)
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
      h2o_fac  = WK(1)/Wtot
      Rself    =     h2o_fac  * RHOave * 1.e-20 * xself
      Rfrgn    = (1.-h2o_fac) * RHOave * 1.e-20 * xfrgn
      Rfrgn_aj =     h2o_fac  * RHOave * 1.e-20 * xfrgn
C
C=======================================================================
C                             SELF

C     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
c
      if ((V2.gt.-20.0).and.(V1.lt.20000.)) then
c
            CALL SL296 (V1C,V2C,DVC,NPTC,SH2OT0)                          F00410
            CALL SL260 (V1C,V2C,DVC,NPTC,SH2OT1)                          F00420
C                                                                         F00440
c           Loop calculating self continuum optical depth
C
            TFAC = (TAVE-T0)/(260.-T0)                                    F00540
C
            DO 20 J = 1, NPTC                                             F00560
               VJ = V1C+DVC* REAL(J-1)                                    F00570
               SH2O = 0.                                                  F00580
               IF (SH2OT0(J).GT.0.) THEN                                  F00590
                  SH2O = SH2OT0(J)*(SH2OT1(J)/SH2OT0(J))**TFAC            F00600
                  SFAC = 1.
c
                  IF (VJ .GE. 820. .AND. VJ .LE. 960.) THEN
                     JFAC = (VJ - 820.)/10. + 0.00001
                     SFAC = XFACREV(JFAC)
                  ENDIF
C                                                                         F00630
                  SH2O = SFAC * SH2O
c
               ENDIF
C              ---------------------------------------------------------
c
               cself(j) = WK(1)*(SH2O*Rself)
C                                                                         F00720
c********************************************
               v1h=V1C
               dvh=DVC
               npth=NPTC
c
               csh2o(j)=1.e-20 * sh2o * xself  
c********************************************
C                                                                         F00720
C              ---------------------------------------------------------
C              Radiation field                                            F00730
C                                                                         F00740
               IF (JRAD.EQ.1) cself(j) = cself(j)*RADFN(VJ,XKT)                   F00750
C              ---------------------------------------------------------

 20         CONTINUE                                                      F00760
C
c           Interpolate to total optical depth grid

            CALL XINT (V1C,V2C,DVC,cself,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)      F00770

         endif
C                                                                         F00780
C=======================================================================
C                             FOREIGN
C
C--------------------------------------------------------------------
C
C        Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
c
         if ((V2.gt.-20.0).and.(V1.lt.20000.)) then

C--------------------------------------------------------------------
C     *****    Continuum Correction Patches    ********
C--------------------------------------------------------------------

            V0F1 = 370.
            HWSQF1 = 190.**2
            BETAF1 = 1.e-08 
            FACTRF1 = -0.25
C
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
               FH2O(J)=FH2O(J)*FSCAL
C     
               C(J)        = WK(1) * (FH2O(J)*RFRGN)

               cfrgn_aj(j) = wk(1) * fh2o(j) * rfrgn_aj
C                                          
c********************************************
               cfh2o(j)=1.e-20 * fh2o(j) * xfrgn
c********************************************
C              ---------------------------------------------------------
C              Radiation field                                                  
C                                                                    
               IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)               
C              ---------------------------------------------------------
C
 24         CONTINUE                                                  
C
            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)
C
C           ------------------------------------------------------------
c
            if  (icflg.eq.1) then

               do j=1,nptc

                  if (jrad.eq.1) then
                     vj = v1c + dvc*real(j-1)
                     c(j) = (cself(j)-cfrgn_aj(j)) * radfn(vj,xkt)
                  else
                     c(j) =  cself(j)-cfrgn_aj(j)
                  endif
               enddo

               Call XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)

            endif

C           ------------------------------------------------------------
C                                                                         F00780
         endif

c compute H2O continuum derivatives
         if (icflg.eq.1) then

            if ((V2.gt.-20.0).and.(V1.lt.20000.)) then

               do j=1,nptc
c w.r.t. ln(q)
c                dqh2o must be returned with the radiation field included
               enddo
            else
               write(ipr,*) 'WARNING:  ANALYTIC DERIVATIVE / CONTNM'
               write(ipr,*) ' v1 - v2 out of range for H2O continuum'
               write(ipr,*) '  (error trap - 1)'
            endif               ! v1,v2 range
         endif                  ! icflg


      RETURN
      END IF
c========================================================================
