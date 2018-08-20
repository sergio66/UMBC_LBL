c ifort -extend-source 132 -O2 -names lowercase -check bounds N2_routines.f -o N2WV_routines.x

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NVAL=901,NSIGVECT=10000)
      COMMON/N2CIA/SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),
     ,   B0H2O(NVAL),BETA0H2O(NVAL)    

      INTEGER iNpts,iI
      DOUBLE PRECISION SIGVECT(NSIGVECT),X(NSIGVECT),Y(NSIGVECT),Z(NSIGVECT),CTN2
      DOUBLE PRECISION A(NSIGVECT)
      
      PTOT = 1.0        !! atm
      PN2 = PTOT*0.79   !! atm
      
      PH2O = 1.0D-1     !! atm
      PH2O = 1.0D-2     !! atm
      PH2O = 5.0D-1     !! atm
      PH2O = 1.0D-7     !! atm
      PH2O = 5000/1D6   !! tropical WV mix ratio at surface ~ 5000 ppm
      
      T    = 296.0      !! K
      
      SIGMA0 = 1930.0D0
      SIGMAF = 2830.0D0
      DSIGMA = 1.0D-1

      iNpts = (SIGMAF-SIGMA0)/DSIGMA + 1
!      print *,iNpts,NVAL

      CALL LECN2
      
      DO iI = 1,iNpts
        SIGVECT(iI) = SIGMA0 + (iI-1)*DSIGMA
        X(iI) = CTN2(SIGVECT(iI),PN2,PH2O,PTOT,T,-1)         ! WV/N2
        Y(iI) = CTN2(SIGVECT(iI),PN2,PH2O,PTOT,T,0)          ! N2/N2
        Z(iI) = CTN2(SIGVECT(iI),PN2,PH2O,PTOT,T,+1)         ! N2/N2 + WV/N2
!        write(6,10) iI,SIGVECT(iI),X(iI),Y(iI),Z(iI)
      END DO

      CALL CTN2_vect(SIGVECT,PN2,PH2O,PTOT,T,-1,A,iNpts)
      DO iI = 1,iNpts
        write(6,20) iI,SIGVECT(iI),X(iI),Y(iI),Z(iI),A(iI),A(iI)/X(iI)
      END DO
      
  10  FORMAT(I5,' ',F10.3,' ',3(ES12.5,' '))
  20  FORMAT(I5,' ',F10.3,' ',4(ES12.5,' '),F10.3)  
  
      END
c************************************************************************
      SUBROUTINE LECN2
C
C This routine reads the data that enable calculations
C of the N2-N2 and N2-H2O collision-induced absorptions
C in the fundamental band of N2
C
C These data have been generated as explained in the paper 
C "Indirect influence of of humidity on atmospheric emission 
C  and transmission spectra near 4 microns"
C
C Creted by J-M Hartmann, March 2018
C jmhartmann@lmd.polytechnique.fr
C
C
C Number of tabulated values
      PARAMETER (NVAL=901)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/N2CIA/SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),
     ,   B0H2O(NVAL),BETA0H2O(NVAL)    
C      
C OOpen files and read data
C
      OPEN(UNIT=3,FILE='CT-N2.N2',
     ,          STATUS='OLD',FORM='FORMATTED')
C Read header then read data
          DO 1 I=1,11
          READ(3,*)
  1       CONTINUE
      DO 2 I=1,NVAL
      READ(3,*)SIGRF(I),B0air(I),BETA0air(I)
  2   CONTINUE
      CLOSE(3)
C
      OPEN(UNIT=3,FILE='CT-N2.H2O',
     ,          STATUS='OLD',FORM='FORMATTED')

C Read header then read data
          DO 3 I=1,11
          READ(3,*)
  3       CONTINUE
      DO 4 I=1,NVAL
      READ(3,*)SIGRF(I),B0H2O(I),BETA0H2O(I)
  4   CONTINUE
      CLOSE(3)

      RETURN
      END

c************************************************************************

	DOUBLE PRECISION FUNCTION 
     *         CTN2(SIGMA,PN2,PH2O,PTOT,T,iWhich)
C
C
C This routine computes the absorption coefficient in the collision
C nnduced fundamental absorption band of N2 for air in the presence
C of some humidity.using the data that have been read by Subroutine "LECN2"
C
C     The arguments and their units are the following
C  input        
C    Sigma: wavenumber in units of "1/cm" (reciprocal centimeter)
C    PN2  : partial pressure of N2 in units of "atm" (atmosphere)
C    PH2O : partial pressure of H2O in units of "atm" (atmosphere)
C    PTOT : total pressure in units of "atm" (atmosphere)
C     T   : temperature in units of "K" (Kelvin)
C  iWhich : +1 for N2/N2 and N2/WV and 0 for N2/N2 and -1 for N2/WV
C  output        
C    CTN2:  absorption coefficient for the considered conditions
C           in units of "1/cm" (reciprocal centimeter). Hence, for
C           an optical path of legth L, the transmission is
C           trans = exp(-CTN2*L) where L has to be in centimeer units 
C
C Important: if the nominal N2 vmr of 0.781 is used to compute
C PN2, then PN2=0.781*(PTOT-PH2O) and NOT PN2=0.781*PTOT
C
C This routine uses a model that is described in the paper 
C "Indirect influence of of humidity on atmospheric emission 
C  and transmission spectra near 4 microns"
C
C Creted by J-M Hartmann, March 2018
C jmhartmann@lmd.polytechnique.fr
C
C
C
C Number of tabulated values and wavenumber step
      PARAMETER (NVAL=901 , StpSig=1.D0)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
C Tabulated values of the data for N2-N2 and N2-H2O
C These have been read by routine "LECN2"
C
      COMMON/N2CIA/SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),
     ,   B0H2O(NVAL),BETA0H2O(NVAL)    
C
      DATA T0/273.16D0/
      DATA TREF/296.D0/
C
C
      CTN2=0.D0
      IF ( T.GT.350.D0 ) RETURN
      IF((SIGMA.LT.SIGRF(1)).OR.(SIGMA.GT.SIGRF(NVAL)))RETURN
C 
C Compute the N2-N2 and N2-H2O CIA absorption coefficients
C (Bair and BH2o, respectively) by using the exponential Temperature
C dependence from the tabulated values and a liner interpolation versus
C wavenumber using the two sorrounding points (INF and SUP)
c  The CIA at T is computed from B0*exp[BETA0*(1/Tref-1/T)]
C
      IINF=INT( (SIGMA-SIGRF(1)+0.1D-4)/StpSig ) + 1
      IINF=MIN0(IINF,NVAL-1)
      ISUP=IINF+1
      D1ST=(1.D0/TREF)-(1.D0/T)
      BINFair=B0air(IINF)*DEXP(BETA0air(IINF)*D1ST)
      BSUPair=B0air(ISUP)*DEXP(BETA0air(ISUP)*D1ST)
      Bair=BINFair+(BSUPair-BINFair)*(SIGMA-SIGRF(IINF))/StpSig
      BINFH2O=B0H2O(IINF)*DEXP(BETA0H2O(IINF)*D1ST)
      BSUPH2O=B0H2O(ISUP)*DEXP(BETA0H2O(ISUP)*D1ST)
      BH2O=BINFH2O+(BSUPH2O-BINFH2O)*(SIGMA-SIGRF(IINF))/StpSig
C
C Then correct Bair by introducing the contribution of the N2-O2 CIA
C
      Bair=Bair*(0.79 + 0.21*(1.294D0-0.4545D0*T/TREF))
C
C Switch from pressures (in atm) to densities (in amagat)
C and compute CIA by combining dry air (N2+O2) and H2O
C contributions
C
      DTOT=PTOT*(T0/T)
      DN2=PN2*(T0/T)
      DH2O=PH2O*(T0/T)

      IF (iWhich .EQ. 1) THEN
        !! both N2/N2 and N2/WV
        CTN2=DN2*(Bair*(DTOT-DH2O)+BH2O*DH2O)
      ELSEIF (iWhich .EQ. 0) THEN
        !! only N2/N2         
        CTN2=DN2*(Bair*(DTOT-DH2O)+0.0D0*BH2O*DH2O)
      ELSEIF (iWhich .EQ. -1) THEN
        !! only N2/WV                  
         CTN2=DN2*(Bair*(DTOT-DH2O)*0.0D0+BH2O*DH2O)
      END IF
      
      RETURN
      END

c************************************************************************
	SUBROUTINE CTN2_vect(SIGVECT,PN2,PH2O,PTOT,T,iWhich,CTN2,iNpts)
C
C
C This routine computes the absorption coefficient in the collision
C nnduced fundamental absorption band of N2 for air in the presence
C of some humidity.using the data that have been read by Subroutine "LECN2"
C
C     The arguments and their units are the following
C  input        
C    SigVECT: wavenumber in units of "1/cm" (reciprocal centimeter)
C    PN2  : partial pressure of N2 in units of "atm" (atmosphere)
C    PH2O : partial pressure of H2O in units of "atm" (atmosphere)
C    PTOT : total pressure in units of "atm" (atmosphere)
C     T   : temperature in units of "K" (Kelvin)
C  iWhich : +1 for N2/N2 and N2/WV and 0 for N2/N2 and -1 for N2/WV
C  output        
C    A   :  absorption coefficient for the considered conditions
C           in units of "1/cm" (reciprocal centimeter). Hence, for
C           an optical path of legth L, the transmission is
C           trans = exp(-CTN2*L) where L has to be in centimeer units 
C
C Important: if the nominal N2 vmr of 0.781 is used to compute
C PN2, then PN2=0.781*(PTOT-PH2O) and NOT PN2=0.781*PTOT
C
C This routine uses a model that is described in the paper 
C "Indirect influence of of humidity on atmospheric emission 
C  and transmission spectra near 4 microns"
C
C Creted by J-M Hartmann, March 2018
C jmhartmann@lmd.polytechnique.fr
C
C
C
C Number of tabulated values and wavenumber step
      PARAMETER (NVAL=901 , StpSig=1.D0)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION CTN2(*),SIGVECT(*)
C      
C Tabulated values of the data for N2-N2 and N2-H2O
C These have been read by routine "LECN2"
C
      COMMON/N2CIA/SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),
     ,   B0H2O(NVAL),BETA0H2O(NVAL)    
C
      DATA T0/273.16D0/
      DATA TREF/296.D0/
C
C
      DO iI = 1,iNpts
        SIGMA    = SIGVECT(iI)
        CTN2(iI) = 0.D0	
        IF (( T.GT.350.D0 ) .OR. (SIGMA.LT.SIGRF(1)) .OR.(SIGMA.GT.SIGRF(NVAL))) THEN
          CTN2(iI) = 0.D0
        ELSE
C Compute the N2-N2 and N2-H2O CIA absorption coefficients
C (Bair and BH2o, respectively) by using the exponential Temperature
C dependence from the tabulated values and a liner interpolation versus
C wavenumber using the two sorrounding points (INF and SUP)
c  The CIA at T is computed from B0*exp[BETA0*(1/Tref-1/T)]
C
          IINF=INT( (SIGMA-SIGRF(1)+0.1D-4)/StpSig ) + 1
          IINF=MIN0(IINF,NVAL-1)
          ISUP=IINF+1
          D1ST=(1.D0/TREF)-(1.D0/T)
          BINFair=B0air(IINF)*DEXP(BETA0air(IINF)*D1ST)
          BSUPair=B0air(ISUP)*DEXP(BETA0air(ISUP)*D1ST)
          Bair=BINFair+(BSUPair-BINFair)*(SIGMA-SIGRF(IINF))/StpSig
          BINFH2O=B0H2O(IINF)*DEXP(BETA0H2O(IINF)*D1ST)
          BSUPH2O=B0H2O(ISUP)*DEXP(BETA0H2O(ISUP)*D1ST)
          BH2O=BINFH2O+(BSUPH2O-BINFH2O)*(SIGMA-SIGRF(IINF))/StpSig
C
C Then correct Bair by introducing the contribution of the N2-O2 CIA
C
          Bair=Bair*(0.79 + 0.21*(1.294D0-0.4545D0*T/TREF))
C
C Switch from pressures (in atm) to densities (in amagat)
C and compute CIA by combining dry air (N2+O2) and H2O
C contributions
C
          DTOT=PTOT*(T0/T)
          DN2=PN2*(T0/T)
          DH2O=PH2O*(T0/T)

          IF (iWhich .EQ. 1) THEN
            !! both N2/N2 and N2/WV
            CTN2(iI)=DN2*(Bair*(DTOT-DH2O)+BH2O*DH2O)
          ELSEIF (iWhich .EQ. 0) THEN
            !! only N2/N2         
            CTN2(iI)=DN2*(Bair*(DTOT-DH2O)+0.0D0*BH2O*DH2O)
          ELSEIF (iWhich .EQ. -1) THEN
            !! only N2/WV                  
             CTN2(iI)=DN2*(Bair*(DTOT-DH2O)*0.0D0+BH2O*DH2O)
          END IF
        END IF
      END DO
      
      RETURN
      END

c************************************************************************
