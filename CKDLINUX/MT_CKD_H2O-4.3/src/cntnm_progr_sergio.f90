! similar to MT_CKD3.2/cntnm/src/cntnm_progr_sergio.f90
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! input pave,ppave in mb,
!       tave       in K
!       num_kmoles in kilomoles/cm2
!       v1absIN,v2absIN,dvabsIN in cm-1
!       irand is an integer random placeholder
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
!     The mt_ckd water vapor continuum is a completely new continuum  
!     formulation based on a collision induced  component and a sub-Lorentzian 
!     line wing component.  Both the water vapor continuum coefficients and 
!     those for other molecules are constrained to agree with accurate
!     measurements of continuum absorption in the spectral regions where such
!     measurements exist.
!
!     This is an updated version of the continuum program:
!     this version provides optical depths on file CNTNM.OPTDPT as before:
!     it also provides the continuum coefficients on file  WATER.COEF
!
!     the length of the header records may vary by version:
!         in this version the WATER.COEF header information is 47 records 
!         in this version the CNTNM.OPTDT header information is 34 records 
!
!     presumably the user will want to create an input file to address
!     individual requirements
!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! this is from MT_CKD_H2O-4.3/src/run8_drive_mt_ckd_h2o.f90

  SUBROUTINE CALCON_MTCKD_43_loc(                           &
                paveIN,ppaveIN,taveIN,num_kmolesIN,         &
                v1absIN,v2absIN,dvabsIN,                    &
                raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,   &
                raFreq,raAbs,iNumPts,                       &
                iGasID,rXSelf,rXForn,rXozone,rXoxygen,rXnitrogen, &
		irand)

   use mt_ckd_h2o
   USE phys_consts
!   use lblparams

   implicit none
!
!      IMPLICIT REAL*8           (V)                                     ! F00030

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin junk1.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     
! this is from MT_CKD4.3/cntnm/src/cntnm_progr_sergio.f90

      include '/home/sergio/SPECTRA/FORTRANLINUX/max.incF90'

! see ../build/addlibs.inc
!   include "netcdf.inc"
   include "/usr/ebuild/installs/software/netCDF-Fortran/4.6.1-iimpi-2023a/include/netcdf.inc"
 
      integer kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts

!      PARAMETER (kMaxPtsDVS = 1000)     !!! at coarser spacing
!      PARAMETER (kMaxPts    = 100000)   !!! at user spacing
      PARAMETER (kMaxPtsDVS = 250000)    !!! at coarser spacing
      PARAMETER (kMaxPts    = 2500010)   !!! at user spacing
      real*8 num_kmolesIN,ppaveIN,taveIN,paveIN
      real*8 v1absIN,v2absIN,dvabsIN
      real*8 raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),raFornDVS(kMaxPtsDVS)
      real*8 raFreq(kMaxPts),raAbs(kMaxPts)
      REAL*8    rXSelf,rXForn,rXozone,rXoxygen,rXnitrogen
      INTEGER iGasID,irand
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin junk1.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     

   double precision :: p_atm,t_atm,h2o_frac
   double precision :: dwv
   double precision :: wv1,wv2
   real*8, dimension(:),allocatable :: sh2o,fh2o
   double precision,dimension(:),allocatable :: wvn
   logical(kind=8) :: radflag=.TRUE.
   character(len=81) :: mt_version

   integer :: nwv,i
   integer :: istat,ncid,id_wv,idvar(3)
   character :: FRGNX
   integer :: iprint
   real*8 gamma_T,gamma_Ts

      !INTEGER n_absrb,nc
      ! this is set in MT_CKD3.2/cntnm/src/lblparams.f90 
      !      parameter (n_absrb=150050)
       
       INTEGER nc,n_absrb
!       parameter (nc=160000,N_ABSRB=2500010)
       parameter (N_ABSRB=2500010)

       integer nptabs,jrad
       real*8 pave,tave
       real*8 ABSRB(n_absrb),VMRH2O
       real*8 V1,V2,xlength,VI
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin junk1.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     

  pave = paveIN
  tave = taveIN
 
  IF (iGasID .EQ. 1) THEN
    VMRH2O = ppaveIN/paveIN         !! assume we are doing WV
  ELSE
    VMRH2O = 0.0
  END IF
  
  xlength = 1000 * num_kmolesIN   !! change from kmoles/cm2 to moles/cm2
  xlength = xlength * 10000       !! change to moles/m2
  xlength = xlength * 8.3144598 * tave/(ppaveIN*100)  !! pressure changed from mb to N/m2
  xlength = xlength * 100         !! change from m to cm
  print *,'  '
  print *,' inside cntnm_progr_sergio.f90'
  print *,' xlength = ',xlength/100,' meters'
  
   mt_version = '4.3'
   wv2 = v2absIN
   wv1 = v1absIN
   dwv = DVABSIN
   nwv = (wv2-wv1)/dwv+1
   allocate(wvn(nwv))
   allocate (sh2o(nwv))
   allocate (fh2o(nwv))

    write(*,'(A,I3,3(F12.5),ES12.5)') 'cntnm_progr_sergio.f90  here1 ',iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
    write(*,'(A,3(F12.5),I10)') 'cntnm_progr_sergio.f90  here2 ',v1absIN,v2absIn,dvabsIN,NWV

!   print *,'wv1, wv2, dwv, nwv = ',wv1, wv2, dwv, nwv

      nptabs = nwv
      do 85 i=1,nptabs
         absrb(i) =0.
 85   continue


  970 FORMAT (/,29x, 'P(MB)',7X,'T(K)', //,23x,0P,F12.3,F9.2)
  975 FORMAT (/, 9X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',//, &
              8(1X,A6,3X))
  980 FORMAT (/,1P,8E10.3,//)
!
      jrad=1
!
      v1 = v1absIN
      v2 = v2absIN
!
!      CALL CONTNM(JRAD,rXSelf,rXforn)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! this is /home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/src/drive_mt_ckd_h2o.f90

!   namelist /mt_ckd_input/ p_atm,t_atm,h2o_frac,wv1,wv2,dwv
!   read (*,mt_ckd_input)


   p_atm = pave
   t_atm = tave
   h2o_frac = VMRH2O
   
   !!! I have to use dble since mt_ckd_h2o was compiled using   make makefile.bkah intelDBL
   !!! I have to use dble since mt_ckd_h2o was compiled using   make makefile.bkah intelDBL
   call mt_ckd_h2o_absco (dble(p_atm),dble(t_atm),dble(h2o_frac),wv1,wv2, &
                          dble(dwv),sh2o,fh2o, &
                          FRGNX,radflag=radflag,mt_version=mt_version)
   wvn = wv1+[(i,i=0,nwv-1)]*dwv
   
!! see /home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/src/drive_mt_ckd_h2o.f90
      iNumPts = NPTABS
      gamma_T = 1.0   !!! correction from 296 K to tave
      gamma_Ts = 1.0  !!! part of pself

      print *,'mt_version in cntnm_progr_sergio.f90 = ',mt_version
      print *,' nwv and NPTABS = ',nwv,NPTABS
      print *,' gamma_T, gamma_Ts = ',gamma_T,gamma_Ts
      print *,' h2o_frac num_kmolesIN = ',h2o_frac, num_kmolesIN
      print *,' p_atm, tatm  = ',p_atm, t_atm 

      DO I=1,NPTABS
        VI = V1ABSIN+FLOAT(I-1)*DVABSIN
        raFreq(I) = VI
        raAbs(I)  = sh2o(I)
        raAbs(I)  = sh2o(I) + fh2o(I)
        raAbs(I)  = gamma_T*num_kmolesIN*(rXSelf*sh2o(I)*num_kmolesIN*gamma_Ts + rXforn*fh2o(I)*(paveIN-ppaveIN))
!        write(*,'(A,I15,F12.5,ES12.5)') 'in MT_CKD_H2O-4.3/src/cntnm_progr_sergio.f90 ',I,raFreq(I),raAbs(I)
      END DO
  910 FORMAT(F20.8,1P,E20.12)
   
end
  
