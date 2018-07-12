      integer iGasID,kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts
      PARAMETER (kMaxPtsDVS = 1000)     !!! at coarser spacing
      PARAMETER (kMaxPts    = 100000)   !!! at user spacing
      real num_kmolesIN,ppaveIN,taveIN,paveIN
      real v1absIN,v2absIN,dvabsIN
      real raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),raFornDVS(kMaxPtsDVS)
      real raFreq(kMaxPts),raAbs(kMaxPts)

      print *,'Units as per run8 : gasID atm atm K kmoles/cm2'
      print *,'Enter [iGasID paveIN ppaveIN taveIN  num_kmolesIN] : '
      read *,iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
      paveIN  = paveIN * 1013.25
      ppaveIN = ppaveIN * 1013.25

      print *,'Enter [v1 v2 dv] : '
      read *,v1absIN,v2absIN,dvabsIN

      !! [Self Forn   o3 o2 n2]      
      IF (iGasID .EQ. 1) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,     &
                               v1absIN,v2absIN,dvabsIN,                  &
                               raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS, &
                               raFreq,raAbs,iNumPts,                     &
                               iGasID, 1.0, 1.0, 0.0, 0.0, 0.0)
      ELSEIF (iGasID .EQ. 3) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,     &
                               v1absIN,v2absIN,dvabsIN,                  &
                               raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS, &
                               raFreq,raAbs,iNumPts,                     &
                               iGasID, 0.0, 0.0, 1.0, 0.0, 0.0)
      ELSEIF (iGasID .EQ. 7) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,     &
                               v1absIN,v2absIN,dvabsIN,                  &
                               raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS, &
                               raFreq,raAbs,iNumPts,                     &
                               iGasID, 0.0, 0.0, 0.0, 1.0, 0.0)          &
      ELSEIF (iGasID .EQ. 22) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,     &
                               v1absIN,v2absIN,dvabsIN,                  &
                               raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS, &
                               raFreq,raAbs,iNumPts,                     &
                               iGasID, 0.0, 0.0, 0.0, 0.0, 1.0)
      ELSE
        print *,'Error : need iGasID = 1,3,7,22 for WV,O3,O3,N2, not',iGasID
        STOP
      END IF

      END

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE XYZCALCON_MTCKD_32_loc(                                &
                      paveIN,ppaveIN,taveIN,num_kmolesIN,            &
                      v1absIN,v2absIN,dvabsIN,                       &
                      raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,      &
                      raFreq,raAbs,iNumPts,                          &
                      iGasID,rXSelf,rXForn,rXozone,rXoxygen,rXnitrogen)

! input pave,ppave in mb, 
!       tave       in K 
!       num_kmoles in kilomoles/cm2
!       v1absIN,v2absIN,dvabsIN in cm-1
!
! output coarser output wavenumber (at spacing 10 cm-1) : raFreqDVS
!        self continuum coeffs at coarse spacing        : raSelfDVS
!        forn continuum coeffs at coarse spacing        : raFornDVS
!        number of pts         at coarse spacing        : iNumPtsDVS
!
!        user specified wavenumber : raFreq
!        optical depth             : raAbs
!        num of pts                : iNumPts
!
!   The mt_ckd water vapor continuum is a completely new continuum
!   formulation based on a collision induced component and a sub-Lorentzian 
!   line wing component. Both the water vapor continuum coefficients and 
!   those for other molecules are constrained to agree with accurate
!   measurements of continuum absorption in the spectral regions where such
!   measurements exist.
!
!     This is an updated version of the continuum program:
!     this version provides optical depths on file CNTNM.OPTDPT as before:
!     it also provides the continuum coefficients on file  WATER.COEF
c
!     the length of the header records may vary by version:
!         in this version the WATER.COEF header information is 47 records 
!         in this version the CNTNM.OPTDT header information is 34 records 
!
!     presumably the user will want to create an input file to address
!     individual requirements
!
!
      IMPLICIT REAL*8           (V)                                     ! F00030

      REAL    rXSelf,rXForn,rXozone,rXoxygen,rXnitrogen
      integer kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts
      PARAMETER (kMaxPtsDVS = 1000)     !!! at coarser spacing
      PARAMETER (kMaxPts    = 100000)   !!! at user spacing
      real num_kmolesIN,ppaveIN,taveIN,paveIN
      real v1absIN,v2absIN,dvabsIN
      real raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),raFornDVS(kMaxPtsDVS)
      real raFreq(kMaxPts),raAbs(kMaxPts)
      INTEGER iGasID

! this is part of lblparams.f90
!      INTEGER n_absrb,nc
!      parameter (n_absrb=150050)
!      parameter (nc=160000)

