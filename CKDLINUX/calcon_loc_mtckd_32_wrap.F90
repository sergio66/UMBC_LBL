! Copyright 2018
! University of Maryland Baltimore County 
! All Rights Reserved

! these water subroutines calculate the continuum for WATER(ID=1)
! and are based entirely in the KCARTA routines, except NO jacobians

! similarly the O2 and N2 subroutines are based ENTIRELY on GENLN2 stuff

! note we are using        REAL CON(*),WHICHLAYER instead of
!                          REAL CON(KMAXLAYER,KMAXPTS)

! also real ------> real*8

! because of memeory problems 
! and so we do        DO 10 ILAY=whichlayer,whichlayer
! instead of          DO 10 ILAY=1,NLAY

! NOTE THIS IS SAME BASE CODE AS CALCONWATER.F --- THE ONLY DIFFERENCE
! IS THAT THE INCLUDE FILES ARE DIFFERENT!!!!!

!MODULE calconwater_loc_ckd3p2
!IMPLICIT NONE
!CONTAINS

!************************************************************************
       include '/home/sergio/SPECTRA/CKDLINUX/dspline.F90'
!************************************************************************
! this is MT_CKD3.2
! directly taken from AER MT_CKD3.2

      SUBROUTINE calcon_loc_mtckd_32_wrap(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY, &
             T,P, PARTP,AMT,rXSelf,rXForn,whichlayer)

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
      real*8 raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),raFornDVS(kMaxPtsDVS)
      real*8 raFreq(kMaxPts),raAbs(kMaxPts)
      integer ii

!      print *,'Units as per run8 : gasID atm atm K kmoles/cm2'
!      print *,'Enter [iGasID paveIN ppaveIN taveIN  num_kmolesIN] : '
!      read *,iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
      iGasID = IDGAS
      paveIN  = P(whichlayer)     * 1013.25  !! change from atm to mb
      ppaveIN = PARTP(whichlayer) * 1013.25  !! change from atm to mb
      taveIN  = T(whichlayer)
      num_kmolesIN = AMT(whichlayer)

!      print *,'Enter [v1 v2 dv] : '
!      read *,v1absIN,v2absIN,dvabsIN
      v1absIN = freq(1)
      v2absIn = freq(NFREQ)
      dvabsIN = FSTEP !! note this is "low" resolution eg 1 cm-1 wide

!      print *,'here1 ',iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
!      print *,'here2 ',v1absIN,v2absIn,dvabsIN,NFREQ

      !! [Self Forn   o3 o2 n2]      
      IF (iGasID .EQ. 1) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN, &
                      v1absIN,v2absIN,dvabsIN,                       &
                      raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,      &
                      raFreq,raAbs,iNumPts,                          &
                      iGasID, rXSelf, rXForn, 0.0d0, 0.0d0, 0.0d0)
      ELSE
        print *,'Error : need iGasID = 1 for WV not',iGasID
        STOP
      END IF
      
      !! now interp raAbs onto CON
!      print *,'now need to interp iNumPts onto NFREQ',iNumPts,NFREQ
!      print *,raFreq(1),raAbs(1),raFreq(iNumPts),raAbs(iNumPts)
!      print *,FREQ(1),FREQ(2),FREQ(3),FREQ(NFREQ)
!      print *,raFREQ(1),raFREQ(2),raFREQ(3),FREQ(iNumPts)
      !CALL xspl(raFreq,raAbs,iNumPts,FREQ,CON,NFREQ)
      DO ii = 1,NFREQ
        FREQ(ii) = FREQ(ii)
        CON(ii)  = raAbs(ii)
	print *,ii,FREQ(ii),raFreq(ii)
      END DO

      RETURN
      END

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      include 'MT_CKD3.2/cntnm/src/cntnm_progr_sergio.f90'

!END MODULE calconwater_loc_ckd3p2

