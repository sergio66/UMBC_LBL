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
!!! see mt_ckd_h2o_module.f90
!    The MT_CKD reference continuum coefficients are read from the netCDF file absco-ref_mt_ckd_h2o.nc.
!    The coefficients have units cm^2/molecule/cm^-1 and are normalized to a reference density
!    (P=1013 mb, T=296K) of the relevant collisional partner (self - water vapor, foreign - all
!    gases except water vapor). The radiation field (units of cm-1), must be applied to these
!    coefficients if optical depths are computed -- this operation is optional in this routine
!    (depends on radflag).  Once the radiation term is applied, optical depths can be obtained by
!    multiplying the continuum coefficients by the water vapor column amount (molecules/cm^2). (This
!    final step is not performed in this routine.)
!!! see mt_ckd_h2o_module.f90
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
      INTEGER iGasID,irand,iTest
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin junk1.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     

   double precision :: p_atm,t_atm,h2o_frac
   double precision :: dwv
   double precision :: wv1,wv2
   double precision,dimension(:),allocatable :: sh2o,fh2o
   double precision,dimension(:),allocatable :: wvn
   logical(kind=8) :: radflag=.TRUE.    !!! TURN THIS OFF SINCE I DO THIS MESELF
!   logical(kind=8) :: radflag=.FALSE.
   character(len=81) :: mt_version

   integer :: nwv,i
   integer :: istat,ncid,id_wv,idvar(3)
   character :: FRGNX='0'
   integer :: iprint
   real*8 gamma_T,gamma_Ts

!!!!!!!!! see read_module.f90
!    if (FRGNX.EQ.'1') then
!      call readVarNC(ncid,"for_closure_absco_ref",    dat%for_absco_ref)
!    else
!      call readVarNC(ncid,"for_absco_ref",    dat%for_absco_ref)
!    endif
!!!!!!!!! see read_module.f90

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
       real c2,AVOG,rho_rat
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
   !!! note radflag = .TRUE. above YESYESYES
   !!! note radflag = .FALSE.  above  NONON
   !!!   radflag - (optional) if true, multiply by radiation term (default); if false, do not
   call mt_ckd_h2o_absco (dble(p_atm),dble(t_atm),dble(h2o_frac),wv1,wv2, &
                          dble(dwv),sh2o,fh2o, &
                          FRGNX,radflag=radflag,mt_version=mt_version)
   !!! note this routine says rho_rat = (p_atm/dat%ref_press)*(dat%ref_temp/t_atm)
   !!!    sh2o_coeff --> sh2o_coeff * h2o_vmr * rho_rat
   !!!    fh2o_coeff --> fh2o_coeff * (1-h2o_vmr) * rho_rat
   !!! rho_rat = (p_atm/dat%ref_press)*(dat%ref_temp/t_atm)
   rho_rat = (p_atm/1013.00) * (296/t_atm)
   sh2o = sh2o/(rho_rat * h2o_frac)
   fh2o = fh2o/(rho_rat * (1-h2o_frac))

!    self_absco - computed water vapor self continuum absorption coefficients  (cm2/molec)
!    for_absco - computed water vapor foreign continuum absorption coefficients (cm2/molec)
   wvn = wv1+[(i,i=0,nwv-1)]*dwv
   
!! see /home/sergio/SPECTRA/CKDLINUX/MT_CKD_H2O-4.3/src/drive_mt_ckd_h2o.f90
      iNumPts = NPTABS

!! see KCARTA/DOC/kcarta1.22.pdf, Sect 22.5 Cross section and water continuum jacobians
!! OD = q v tanh(c2v/2T)(296/T)(ps CS + pf CF) = gamma(T) q (CS(T) ps + CF (p-ps))
!! but recall   pV = nRT ==> n/V = p/RT 
!!                       ==> q = n/V L = ps/RT L ==> ps = q RT/L
!! Thus  OD = q v tanh(c2v/2T)(296/T)(ps CS + pf CF)

      gamma_T = 1.0   !!! correction from 296 K to tave
      gamma_Ts = 1.0  !!! part of pself

      iTest = +1  !!! output SH2O
      iTest = +2  !!! output FH2O

      iTest = -1  !!! default, ouput OD

      if ((rXSelf .ge. 0.5) .and. (rXforn .ge. 0.5)) then
        iTest = -1  !!! default, ouput OD
      elseif ((rXSelf .lt. 0.5) .and. (rXforn .ge. 0.5)) then
        iTest = +2  !!! default, ouput FH2O
      elseif ((rXSelf .ge. 0.5) .and. (rXforn .lt. 0.5)) then
        iTest = +1  !!! default, ouput SH2O
      else
        print *,'cntnm_progr_sergio.f90 cannot figure out iTest'
        STOP
      end if

      iTest = -1
      c2 = 1.4387863        !! radiation constant
      AVOG = 6.022045E+26   !! molecules/kmol
      !! our units are in kmoles/cm2 ==> AVOG * Q = (molecules/kmol) (kmol/cm2) = molecules/cm2

!! see run8watercontinuum
!! v tanh(c2 v/2T) (296/T) is done by do_local_lineshape_CKD.m
      if (iTest .LT. 0) THEN
        !!! radiation term v tanh(c2v/2T) is taken care of by AER code
!        gamma_T  = 296/t_atm
!        gamma_Ts = 296/t_atm
        DO I=1,NPTABS
          VI = V1ABSIN + FLOAT(I-1)*DVABSIN
          raFreq(I) = VI
          raAbs(I)  = gamma_T*(rXSelf*sh2o(I)*ppaveIN + rXforn*fh2o(I)*(paveIN-ppaveIN))
          IF ((I .EQ. 1) .OR. (I .EQ. NPTABS)) THEN
            write(*,'(A,I15,F12.5,ES12.5)') 'in MT_CKD_H2O-4.3/src/cntnm_progr_sergio.f90 OD final ',I,raFreq(I),raAbs(I)
          END IF
        END DO
        raAbs = raAbs * num_kmolesIN * AVOG/1000     !!! remember self_absco, for_absco  are in cm2/molecule while we send in kmoles/cm2
        !raAbs = raAbs / (raFreq * tanh(c2 * raFreq)/2/t_atm) * (296/T) !! this is redone in do_local_lineshape_CKD.m
      ELSEIF (iTest .EQ. +1) then
        !! basic test
        DO I=1,NPTABS
          VI = V1ABSIN+FLOAT(I-1)*DVABSIN
          raFreq(I) = VI
          raAbs(I) = AVOG * sh2o(I)
          IF ((I .EQ. 1) .OR. (I .EQ. NPTABS)) THEN
            write(*,'(A,I15,F12.5,3(A,ES12.5))') 'in MT_CKD_H2O-4.3/src/cntnm_progr_sergio.f90  SH2O', &
                 I,raFreq(I),' --> ',sh2o(I),' * ',AVOG,' = ',raAbs(I)
          END IF
        END DO
      ELSEIF (iTest .EQ. +2) then
        !! basic test
        DO I=1,NPTABS
          VI = V1ABSIN+FLOAT(I-1)*DVABSIN
          raFreq(I) = VI
          raAbs(I) = AVOG * fh2o(I)
          IF ((I .EQ. 1) .OR. (I .EQ. NPTABS)) THEN
            write(*,'(A,I15,F12.5,3(A,ES12.5))') 'in MT_CKD_H2O-4.3/src/cntnm_progr_sergio.f90  FH2O', & 
                 I,raFreq(I),' --> ',fh2o(I),' * ',AVOG,' = ',raAbs(I)
          END IF
        END DO
      END IF

      write(*,'(A,A)')       'mt_version in cntnm_progr_sergio.f90 = ',mt_version
      write(*,'(A,2I12)')    ' nwv and NPTABS            = ',nwv,NPTABS
      write(*,'(A,2F12.5)')  ' gamma_T, gamma_Ts         = ',gamma_T,gamma_Ts
      write(*,'(A,2ES12.5)') ' h2o_frac num_kmolesIN/cm2 = ',h2o_frac, num_kmolesIN
      write(*,'(A,2F12.5)')  ' p_atm, tatm               = ',p_atm, t_atm 
      write(*,'(A)') ' '

! from DOC/lbl1.tex
! 
! the Ideal Gas Law (with pressure in atmospheres changed to $Nm^{-2}$;
! dividing by $10^{6}$ changes the units to moles per cubic cm, and dividing
! by a further $10^{3}$ changes it to kiloMoles per cubic cm. After
! multiplying this density by the path cell length  (or layer thickness) $L$
! in cm, the final units are in the GENLN2 units of $kiloMoles/cm^{2}$
! 
! Note that the $stren$ parameter is used in units of
! $cm^{-1}/(moleculescm^{-2})$ while the gas amounts $U$ are in units of
! $kiloMoles cm-2$. To get absorption coefficients and optical depths that
! are in the correct units, the program eventually multiplies $stren$ by
! $6.023 \times 10^{26}$, which is the number of molecules per kilomole.

   
end
  
