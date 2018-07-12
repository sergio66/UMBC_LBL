c basically copied from 
c /home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/src/cntnm_progr_sergio.f

c see /home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/src/Readme_adapt_to_run8watercontinuum_code
c see /home/sergio/SPECTRA/CKDLINUX/MT_CKD3.2/cntnm/build/make_cntnm_sergio_run8

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcon_mtckd_32_loc(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY,
     $    T,P, PARTP,AMT,rXSelf,rXForn,whichlayer)

      IMPLICIT NONE

      INTEGER NPTABS

      INTEGER IDGAS, NFREQ, NLAY, whichlayer
      REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*)
      real*8 rXSelf,rXForn

      INTEGER iVers

      integer iGasID,kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts
      PARAMETER (kMaxPtsDVS = 25000)    !!! at coarser spacing
      PARAMETER (kMaxPts    = 2500010)   !!! at user spacing
      real*8 num_kmolesIN,ppaveIN,taveIN,paveIN
      real*8 v1absIN,v2absIN,dvabsIN
      real*8 raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),
     $                           raFornDVS(kMaxPtsDVS)
      real*8 raFreq(kMaxPts),raAbs(kMaxPts)
      integer ii

c      print *,'Units as per run8 : gasID atm atm K kmoles/cm2'
c      print *,'Enter [iGasID paveIN ppaveIN taveIN  num_kmolesIN] : '
c      read *,iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
      iGasID = IDGAS
      paveIN  = P(whichlayer)     * 1013.25  !! change from atm to mb
      ppaveIN = PARTP(whichlayer) * 1013.25  !! change from atm to mb
      taveIN  = T(whichlayer)
      num_kmolesIN = AMT(whichlayer)

c      print *,'Enter [v1 v2 dv] : '
c      read *,v1absIN,v2absIN,dvabsIN
      v1absIN = freq(1)
      v2absIn = freq(NFREQ)
      dvabsIN = FSTEP !! note this is "low" resolution eg 1 cm-1 wide

c      print *,'here1 ',iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
c      print *,'here2 ',v1absIN,v2absIn,dvabsIN,NFREQ

      !! [Self Forn   o3 o2 n2]      
      IF (iGasID .EQ. 1) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                 v1absIN,v2absIN,dvabsIN,NFREQ,
     $                 raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                 raFreq,raAbs,iNumPts,
     $                 iGasID, rXSelf, rXForn, 0.0d0, 0.0d0, 0.0d0)
      ELSE
        print *,'Error : need iGasID = 1 for WV not',iGasID
        STOP
      END IF

      !! now interp raAbs onto CON
c      print *,'now need to interp iNumPts onto NFREQ',iNumPts,NFREQ
c      print *,raFreq(1),raAbs(1),raFreq(iNumPts),raAbs(iNumPts)
c      print *,FREQ(1),FREQ(2),FREQ(3),FREQ(NFREQ)
c      print *,raFREQ(1),raFREQ(2),raFREQ(3),FREQ(iNumPts)
      !CALL xspl(raFreq,raAbs,iNumPts,FREQ,CON,NFREQ)
      DO ii = 1,NFREQ
        CON(ii) = raAbs(ii)
      END DO

      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calconOXYNIT_12_1(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY,T,P, 
     $    PARTP,AMT,whichlayer)

      INTEGER NPTABS

      INTEGER IDGAS, NFREQ, NLAY, whichlayer
      REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*)

      INTEGER iVers

      integer iGasID,kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts
      PARAMETER (kMaxPtsDVS = 25000)    !!! at coarser spacing
      PARAMETER (kMaxPts    = 2500010)   !!! at user spacing
      real*8 num_kmolesIN,ppaveIN,taveIN,paveIN
      real*8 v1absIN,v2absIN,dvabsIN
      real*8 raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),
     $                           raFornDVS(kMaxPtsDVS)
      real*8 raFreq(kMaxPts),raAbs(kMaxPts)
      integer ii

c      print *,'Units as per run8 : gasID atm atm K kmoles/cm2'
c      print *,'Enter [iGasID paveIN ppaveIN taveIN  num_kmolesIN] : '
c      read *,iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
      iGasID = IDGAS
      paveIN  = P(whichlayer)     * 1013.25  !! change from atm to mb
      ppaveIN = PARTP(whichlayer) * 1013.25  !! change from atm to mb
      taveIN  = T(whichlayer)
      num_kmolesIN = AMT(whichlayer)

c      print *,'Enter [v1 v2 dv] : '
c      read *,v1absIN,v2absIN,dvabsIN
      v1absIN = freq(1)
      v2absIn = freq(NFREQ)
      dvabsIN = FSTEP !! note this is "low" resolution eg 1 cm-1 wide
c      print *,'here ',v1absIN,v2absIn,dvabsIN

      !! [Self Forn   o3 o2 n2]      
c      IF (iGasID .EQ. 1) THEN
c        call CALCON_MTCKD_12_1_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
c     $                      v1absIN,v2absIN,dvabsIN,NFREQ,
c     $                      raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
c     $                      raFreq,raAbs,iNumPts,
c     $                      iGasID, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0)
c      ELSEIF (iGasID .EQ. 3) THEN
      IF (iGasID .EQ. 3) THEN
        call CALCON_MTCKD_12_1_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                      v1absIN,v2absIN,dvabsIN,NFREQ,
     $                      raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                      raFreq,raAbs,iNumPts,
     $                      iGasID, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0)
      ELSEIF (iGasID .EQ. 7) THEN
        call CALCON_MTCKD_12_1_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                      v1absIN,v2absIN,dvabsIN,NFREQ,
     $                      raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                      raFreq,raAbs,iNumPts,
     $                      iGasID, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0)
      ELSEIF (iGasID .EQ. 22) THEN
c        print *,'calling calcon_mtckd_12_1_loc'
        call CALCON_MTCKD_12_1_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                      v1absIN,v2absIN,dvabsIN,NFREQ,
     $                      raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                      raFreq,raAbs,iNumPts,
     $                      iGasID, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0)
      ELSE
        print *,'Error : need iGasID = 3,7,22 for WV,O3,O3,N2, not',iGasID
        STOP
      END IF

      print *,'blah 1'
      !! now interp raAbs onto CON
      CALL xspl(raFreq,raAbs,iNumPts,FREQ,CON,NFREQ)
      print *,'blah 2'

      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE CALCON_MTCKD_12_1_loc(
     $                paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                v1absIN,v2absIN,dvabsIN,NFREQ,
     $                raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                raFreq,raAbs,iNumPts,
     $                iGasID,rXSelf,rXForn,rXozone,rXoxygen,rXnitrogen)

c input pave,ppave in mb, 
c       tave       in K 
c       num_kmoles in kilomoles/cm2
c       v1absIN,v2absIN,dvabsIN in cm-1
c
c output coarser output wavenumber (at spacing 10 cm-1) : raFreqDVS
c        self continuum coeffs at coarse spacing        : raSelfDVS
c        forn continuum coeffs at coarse spacing        : raFornDVS
c        number of pts         at coarse spacing        : iNumPtsDVS
c
c        user specified wavenumber : raFreq
c        optical depth             : raAbs
C        num of pts                : iNumPts
C
c   The mt_ckd water vapor continuum is a completely new continuum
c   formulation based on a collision induced component and a sub-Lorentzian 
c   line wing component. Both the water vapor continuum coefficients and 
c   those for other molecules are constrained to agree with accurate
c   measurements of continuum absorption in the spectral regions where such
c   measurements exist.
c
c     This is an updated version of the continuum program:
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
C
      IMPLICIT REAL*8           (V)                                  ! F00030

      REAL*8    rXSelf,rXForn,rXozone,rXoxygen,rXnitrogen
      integer kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts
      PARAMETER (kMaxPtsDVS = 25000)    !!! at coarser spacing
      PARAMETER (kMaxPts    = 2500010)   !!! at user spacing
      real*8 num_kmolesIN,ppaveIN,taveIN,paveIN
      real*8 v1absIN,v2absIN,dvabsIN
      real*8 raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),
     $                           raFornDVS(kMaxPtsDVS)
      real*8 raFreq(kMaxPts),raAbs(kMaxPts)
      INTEGER iGasID,NFREQ

      INTEGER n_absrb,nc
      parameter (n_absrb=2500010)
      parameter (nc=2500010)

c
      real dvabs
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)              

C
      COMMON /CVRCNT/ HNAMCNT,HVRCNT
c
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       F00130
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   F00140
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    F00150
c
      Common /share/ HOLN2
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY  
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           A03020
      COMMON /XCONT/  V1C,V2C,DVC,NPTC,C(nc) 
c
c********************************************
      COMMON /cnth2o/ V1h,V2h,DVh,NPTh,Ch(n_absrb),csh2o(n_absrb),
     $                cfh2o(n_absrb)
c********************************************
c
      COMMON /IFIL/ IRD,IPRcnt,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        F00180
     *              NLTEFL,LNFIL4,LNGTH4                                  F00190

      common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

c------------------------------------
c for analytic derivative calculation
c note: ipts  = same dimension as ABSRB
c       ipts2 = same dimension as C
      parameter (ipts=n_absrb,ipts2=nc)
      common /CDERIV/ icflg,iuf,v1absc,v2absc,dvabsc,nptabsc,
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o,
     &    dqco2C(ipts),dTco2C(ipts),
     &    dqo3C(ipts),dTo3C(ipts),
     &    dqo2C(ipts),dTo2C(ipts),
     &    dqn2C(ipts),dTn2C(ipts)

      real ctmp(ipts2),cself(ipts),cforeign(ipts),ch2o(ipts2)
      real ctmp2(ipts2),ctmp3(ipts2),ctmp4(ipts),ctmp5(ipts)

c------------------------------------
c
      dimension xcnt(7)
c
      equivalence (xself,xcnt(1))
c
      CHARACTER*18 HNAMCNT,HVRCNT
C                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
c
      character*8 holn2
C                                                                         F00120
      REAL XLOSMT
      DATA XLOSMT/2.68675E+19/
C
      RADCN1=2.*PLANCK*CLIGHT*CLIGHT*1.E-07                                 3070
      RADCN2=PLANCK*CLIGHT/BOLTZ                                            3080
c
      icflg = -999
C
      !!! initialize the multipliers (weights or nm_weight) for the 7 imp gases 
      !!! common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
      do 1, i=1,7
c         xcnt(i)=1.     !!use all gases
         xcnt(i)=0.      !!use  no gases, except as given below
 1    continue

      IF (iGASID .EQ. 1) THEN
        xcnt(1) = 1.0      
        xcnt(2) = 1.0      
      ELSEIF (iGASID .EQ. 3) THEN
        xcnt(4) = 1.0      
      ELSEIF (iGASID .EQ. 7)  THEN
        xcnt(5) = 1.0      
        xcnt(6) = 1.0      
      ELSEIF (iGASID .EQ. 22) THEN
        xcnt(5) = 1.0      
        xcnt(6) = 1.0      
      END IF

c      !! zero out the optical depth for all wavenumbers
c      !! this is done again below, after figuring out nptabs
      do 2, i=1,n_absrb
         absrb(i)=0.0
 2    continue

      !!! initializ the mizing ratios for all gases
      do 3, i=1,60
         wk(i)=0.0
 3    continue
c
      ird = 55
      ipr = 66
      ipu = 7
c

c      OPEN (ipr,FILE='CNTNM.OPTDPT')
       OPEN (ipu,FILE='WATER.COEF')

c
cc      print *
cc      print *, '  This version is limited to ', n_absrb-50 ,' values  '
C
C   THIS PROGRAM CALCULATES THE CONTINUUM OPTICAL DEPTH
C         FOR AN HOMOGENEOUS LAYER
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
C
C   THE RESULTS ARE IN ARRAY ABSORB
C
C   NOTE THAT FOR AN ATMOSPHERIC LAYER: 
C
C            WTOT   = XLOSMT * (PAVE/1013) * (273/TAVE) * (PATH LENGTH)
C
C            WBROAD = the column amount for all species not explicitly provided
C
C            WK(M)  = (VOLUME MIXING RATIO) * (COLUMN OF DRY AIR)
C
C
      iprcnt = ipr
                   CALL PRCNTM 
      iprcnt = ipu
                   CALL PRCNTM
C
C   THE FOLLOWING IS AN EXAMPLE FOR A ONE CM PATH (SEE CNTNM.OPTDPT FOR RESULTS)
C

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c      PAVE = 1013.
c      TAVE =  296.
cc
c      VMRH2O = 0.01
cC
c      xlength = 1.
c

c      print *
c      print *,' *** For this program, vmr_h2o is taken ',
c     * 'with respect to the total column ***'

c      print *
c      print *,' read: pressure (mb)  if negative use default values'
c      read *, press_rd
      
c      if (press_rd .gt. 0.) then
c         pave = press_rd
c         print *,' read:   temperature (K)'
c         read *, tave
c         print *,' read:   path length (cm)'
c         read *, xlength
c         print *,' read:   vmr h2o '
c         read *, vmrh2o
c      endif
c      print *,
c     * 'Pressure (mb), Temperature (K), Path Length (cm),    VMR H2O'
c      print 911, pave,tave,xlength,vmrh2o
c 911  format(1x,f13.6,f17.4,f18.4,f12.8)

      pave = paveIN
      tave = taveIN
      VMRH2O = ppaveIN/paveIN         !! assume we are doing WV

      xlength = 1000 * num_kmolesIN   !! change from kmoles/cm2 to moles/cm2
      xlength = xlength * 10000       !! change to moles/m2
      xlength = xlength * 8.31 * tave/(ppaveIN*100)  !! change from mb to N/m2
      xlength = xlength * 100         !! change from m to cm

c      print *,'P (mb), PP(mb), Temp (K), kmols, Path Length (cm) = ',
c     $          pave,ppaveIN,tave,num_kmolesIN, xlength
c      print *,' xlength = ',xlength/100,' meters'

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C
c     It may be preferable to specifiy the water column directly!
c
      WTOT = XLOSMT*(PAVE/1013.)*(273./TAVE)* xlength

      IF (iGasID .EQ. 1) THEN
        W_dry = WTOT * (1.-VMRH2O)
        !ww is column of dry air;  vol. mix. ratios are based on dry air
        WA     = 0.009     * W_dry ! argon
        WN2    = 0.78      * W_dry ! N2
        WK(7)  = 0.21      * W_dry ! O2
        WK(2)  = 345.E-06  * W_dry ! CO2
        !WK(2) = 0.

        !water vapor:
        if (abs(vmrh2o-1.) .lt. 1.e-05) then
           wk(1) = wtot
        else
           WK(1) = VMRH2O * W_dry
           !for our purposes we want  : W_dry + WK(1) = WTOT
           WK(1) = VMRH2O * WTOT
        endif

      ELSE   !! doing continuum for O3, N2 or O2
        W_dry = WTOT
        !ww is column of dry air;  vol. mix. ratios are based on dry air
        WA     = 0.009     * W_dry ! argon
        WN2    = 0.78      * W_dry ! N2
        WK(7)  = 0.21      * W_dry ! O2
        WK(2)  = 345.E-06  * W_dry ! CO2
        !WK(2) = 0.
        wk(1)  = 0.0               !WV
      END IF

      print *,'A ',VMRH2O
      print *,'B ',WTOT,W_DRY,WK(1),W_DRY+WK(1)
      print *,'C ',WK(1),WA,WK(7),WN2
      print *,'D ',WK(1)/(wk(1) + WN2+ WA + WK(7)), PPAVEIN/PAVEIN
      print *,' ' 

      IF (iGasID .EQ. 1) THEN
        ! remember we set WK(i) = 0 for i >= 2, but we need O2 amt!!
        WBROAD = WN2+ WA + WK(7)
      ELSEIF (iGasID .EQ. 3) THEN
        ! remember we set WK(i) = 0 for i=1, i >= 3, but we need O2 amt!!
        WBROAD = WN2+ WA + WK(7)
      ELSE
        !! WTOT in contnm_12_1.F, line 245, will now be set correctly
        WBROAD = WN2+ WA
      END IF

      print *,'before ',iGasID,WK(1),WK(2),WK(3),WK(7),WK(22)

      IF (iGasID .EQ. 1) THEN
        WK(1) = WK(1)
        WK(1) = wtot * ppaveIN/paveIN
        WK(2) = 0.0
        WK(3) = 0.0
        WK(7) = 0.0
        WK(22) = 0.0
        WN2   = 0.0
      ELSEIF (iGasID .EQ. 3) THEN
        WK(3) = WK(1)
        WK(3) = wtot * ppaveIN/paveIN
        WK(1) = 0.0
        WK(2) = 0.0
        WK(7) = 0.0
        WK(22) = 0.0
        WN2   = 0.0
      ELSEIF (iGasID .EQ. 7) THEN
        WK(7) = WK(1)
        WK(7) = wtot * ppaveIN/paveIN
        WK(1) = 0.0
        WK(2) = 0.0
        WK(3) = 0.0
        WK(22) = 0.0
        WK(22) = WK(7)*0.79/0.21   !! assume O2 : N2 is 21 : 779
        WN2   = 0.0
      ELSEIF (iGasID .EQ. 22) THEN
        WN2   = WK(1)
        WN2   = wtot * ppaveIN/paveIN
        WK(22) = WN2
        WK(1) = 0.0
        WK(2) = 0.0
        WK(3) = 0.0
        WK(7) = 0.0
        WK(7) = WN2*0.21/0.79      !! assume O2 : N2 is 21 : 779
      END IF

      print *,'after ',iGasID,WK(1),WK(2),WK(3),WK(7),WK(22),WTOT
c
      NMOL = 7
c

c      V1ABS =    0.
c      V2ABS = 10000.  
c      DVABS =    2.

       V1ABS = v1absIN
       V2ABS = v2absIN
       DVABS = real(dvabsIN)

c ..........................................................
c      write (*,*) '  v1abs,  v2abs,  dvabs  '
c      read  (*,*)    v1abs,  v2abs,  dvabs
c ..........................................................

      NPTABS =  1. + (v2abs-v1abs)/dvabs
      NPTABS = NFREQ

c      print *,V1abs,v2abs,dvabs,NPTABS

      do 85 i=1,nptabs
         absrb(i) =0.
 85   continue

cc
c      WRITE (IPR,970) PAVE,TAVE
c      WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2                             A23490
c      WRITE (IPR,980) (WK(M),M=1,7),WBROAD
c
c      WRITE (IPu,970) PAVE,TAVE
c      WRITE (IPu,975) (HMOLID(I),I=1,7),HOLN2                             A23490
c      WRITE (IPu,980) (WK(M),M=1,7),WBROAD
C
  970 FORMAT (/,29x, 'P(MB)',7X,'T(K)', //,23x,0P,F12.3,F9.2)
  975 FORMAT (/, 9X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',//,
     *         8(1X,A6,3X))
  980 FORMAT (/,1P,8E10.3,//)
C
      jrad=1
c
      v1 = v1abs
      v2 = v2abs

c************************************************************************
      !! rXSelf,rXForn are UMBC_LBL mults for self and forn
      print *,'CALL CONTNM'
      CALL CONTNM(JRAD,rXSelf,rXForn)  
      print *,'END CALL CONTNM'
c************************************************************************

      iNumPts = NPTABS
      print *,'printing out wavnumber/OD at high res',NPTABS
      DO 100 I = 1,NPTABS
        VI=V1ABS+FLOAT(I-1)*DVABS
        raFreq(I) = VI
        raAbs(I)  = ABSRB(I)
c        IF ((I .GT. 24995) .AND. (I .LT. 25015)) THEN
c           print *, 'zdzd',I, nptabs,VI, ABSRB(I)
c           END IF
c        write(*,*) I, nptabs,VI, ABSRB(I) !! can comment this out
c        WRITE (ipr, 910) VI, ABSRB(I) !! can comment this out
 100    CONTINUE        
 910  FORMAT(F10.3,1P,E12.3)

C
c      WRITE (7,920) tave
c  920 FORMAT(//,' self and foreign water vapor continuum coefficients ',/,
c     x       'for  ',f8.2,'K - ',   //,
c     x       ' the self-continuum scales as ( Rself/Ro ) ',/,
c     x       ' the foreign continuum scales as ( (Rtot-Rself)/Ro ) ',/,
c     x       ' where R is the density rho [ R = (P/Po)*(To/T) ]. ',//,  
c     x   10x,'     without radiation field:  ',
c     x   10x,'       with radiation field:   ',      /,
c     x   10x,'      self         foreign     ',
c     x   10x,'      self         foreign     ',      /,
c     x       '    cm-1  ',                    
c     x       '       1/(cm-1 molec/cm**2)    ',
c     x   10x,'       1/(molec/cm**2) '        ,//)
c
      xkt=tave/radcn2

cc this is going to output at DVS = 10 cm-1,  closest points multiples of 10
cc which span v1absIN,v2absIN
cc eg if v1absIN = 605, v2absIN = 655, then continuum at 610,620,630,640,650 
cc are output
cc so may as well send in [1 10000 100] for [v1absIN v2absIN dvabsIN] and
cc get out 1000 points, from 10 cm-1 to 10000 cm-1, then interp them!!!
      iNumPtsDVS = NPTH
      print *,'printing out CKD coeffs at (coarse) 10 cm-1 res',NPTH
      do 200 i=1,npth
        vi=v1h+float(i-1)*dvh
        raFreqDVS(i) = vi
        if (vi.ge.v1abs .and. vi.le.v2abs) then
          radfld=radfn(vi,xkt)
          csh2or=csh2o(i) * radfld
          cfh2or=cfh2o(i) * radfld
c          print *,i,npth,vi,csh2o(i),cfh2o(i)
          write (ipu,930) vi, csh2o(i), cfh2o(i), csh2or, cfh2or  !! can comment this out
          raSelfDVS(i) = csh2o(i)
          raFornDVS(i) = cfh2o(i)
        endif
  200 continue
  930 format(f10.2, 1p, 2e15.4,10x, 1p, 2e15.4)
      CLOSE(ipu)
c
      RETURN
      END 

c**********************************************************************
      Block Data phys_consts
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY 
c
      DATA PI /3.1415926535898 /
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     Pi was obtained from   PI = 2.*ASIN(1.)                             A03980
c
c     units are generally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      A03990
c                            RADCN2 = PLANCK*CLIGHT/BOLTZ                 A04000
      end
c
      BLOCK DATA                                                          A07600
C
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           A03020
C
      DATA ONEPL/1.001/, ONEMI/0.999/, ARGMIN/34./
C                                                                         A07710

      END                                                                 A07720
      BLOCK DATA cntnm
c
      IMPLICIT REAL*8           (V)                                     ! F00030
c
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       F00130
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   F00140
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    F00150
c
      Common /share/ HOLN2
C                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
c
      character*8 holn2
c
      DATA HMOlid/ '  H2O   ' , '  CO2   ' , '   O3   ' , '  N2O   ' ,    A12260
     *             '   CO   ' , '  CH4   ' , '   O2   ' , 53*'        '/
c
      DATA HOLN2 / ' OTHER'/                                              A19810
c
      end
      SUBROUTINE XINT (V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3)         B17520
C                                                                         B17530
C                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! F00030
C                                                                         B17550
C     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED                     B17560
C     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE         B17570
C     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN       B17580
C     INCREMENTS OF DVR3                                                  B17590
C                                                                         B17600
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B17610
      DIMENSION A(*),R3(*)                                                B17620
C                                                                         B17630
      RECDVA = 1./DVA                                                     B17640
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI                                   B17650
      ILO = MAX(ILO,N1R3)                                                 B17660
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI                                      B17670
      IHI = MIN(IHI,N2R3)                                                 B17680
C                                                                         B17690
      DO 10 I = ILO, IHI                                                  B17700
         VI = VFT+DVR3*FLOAT(I-1)                                         B17710
         J = (VI-V1A)*RECDVA+ONEPL                                        B17720
         VJ = V1A+DVA*FLOAT(J-1)                                          B17730
         P = RECDVA*(VI-VJ)                                               B17740
         C = (3.-2.*P)*P*P                                                B17750
         B = 0.5*P*(1.-P)                                                 B17760
         B1 = B*(1.-P)                                                    B17770
         B2 = B*P                                                         B17780
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2        B17790
         R3(I) = R3(I)+CONTI*AFACT                                        B17800
   10 CONTINUE                                                            B17810
C                                                                         B17820
      RETURN                                                              B17830
C                                                                         B17840
      END                                                                 B17850

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION RADFN (VI,XKT)                                             B17860
C                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! F00030
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B17920
C                                                                         B17930
C               LAST MODIFICATION:    12 AUGUST 1991                      B17940
C                                                                         B17950
C                  IMPLEMENTATION:    R.D. WORSHAM                        B17960
C                                                                         B17970
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
C                                     R.D. WORSHAM                        B17990
C                                     J.L. MONCET                         B18000
C                                                                         B18010
C                                                                         B18020
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
C                                                                         B18050
C----------------------------------------------------------------------   B18060
C                                                                         B18070
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
C                                     OFFICE OF ENERGY RESEARCH           B18090
C                                     DEPARTMENT OF ENERGY                B18100
C                                                                         B18110
C                                                                         B18120
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
C                                                                         B18140
C                                             FASCOD3                     B18150
C                                                                         B18160
C                                                                         B18170
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18180
C                                                                         B18190
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B18200
C                                                                         B18210
C      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                         B18220
C                                                                         B18230
      XVI = VI                                                            B18240
C                                                                         B18250
      IF (XKT.GT.0.0) THEN                                                B18260
C                                                                         B18270
         XVIOKT = XVI/XKT                                                 B18280
C                                                                         B18290
         IF (XVIOKT.LE.0.01) THEN                                         B18300
            RADFN = 0.5*XVIOKT*XVI                                        B18310
C                                                                         B18320
         ELSEIF (XVIOKT.LE.10.0) THEN                                     B18330
            EXPVKT = EXP(-XVIOKT)                                         B18340
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT)                           B18350
C                                                                         B18360
         ELSE                                                             B18370
            RADFN = XVI                                                   B18380
         ENDIF                                                            B18390
C                                                                         B18400
      ELSE                                                                B18410
         RADFN = XVI                                                      B18420
      ENDIF                                                               B18430
C                                                                         B18440
      RETURN                                                              B18450
C                                                                         B18460
      END                                                                 B18470

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION RADFN_KCARTA (VI,XKT)                                      B17860

C     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE     B17900
c     SERGIO MACHADO

      IMPLICIT REAL*8           (V)                                     ! F00030

      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY  

C                                                                         B17890
C     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE     B17900
C                                                                         B17910
      XVI = VI                                                            B18240
C                                                                         B18250
      IF (XKT.GT.0.0) THEN                                                B18260
C                                                                         B18270
         XVIOKT = XVI/XKT                                                 B18280
C                                                                         B18290
         IF (XVIOKT.LE.0.01) THEN                                         B18300
            RADFNx = 0.5*XVIOKT*XVI                                        B18310
C                                                                         B18320
         ELSEIF (XVIOKT.LE.10.0) THEN                                     B18330
            EXPVKT = EXP(-XVIOKT)                                         B18340
            RADFNx = XVI*(1.-EXPVKT)/(1.+EXPVKT)                           B18350
C                                                                         B18360
         ELSE                                                             B18370
            RADFNx = XVI                                                   B18380
         ENDIF                                                            B18390
C                                                                         B18400
      ELSE                                                                B18410
         RADFNx = XVI                                                      B18420
      ENDIF                                                               B18430

c kcarta uses AVOG * q * v * tanh(c2 v/2/T) * (296/T) * press
c tanhx = (e^+x - e^-x)/(e^+x + e^-x) = (1-e^(-2x))/(1+e^(2x))
c so we have the [v*tanh(gx)] covered above
c contm_12_1 covers AVOG * q  as well as press = ppress or (total-p)press
c all we now need is 296/T
c      RADFN = RADFNx
c      print *,XKT,RADCN2,XKT*RADCN2
      RADFN = RADFNx * (296.0/(XKT*RADCN2))

      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Include 'contnm_12_1.F'


