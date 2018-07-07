c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

ccccccccc this is NEW N2 and O2 continuum
ccccccccc  copied from LBLRTMv5.10, but modified so that it uses 
ccccccccc  arbitrary N2/O2 mixing ratios

c************************************************************************
ccccccccccc include the common data blocks
c************************************************************************
      include 'brandnew_o2ctt.f'
c************************************************************************
      include 'brandnew_n2ct0.f'
c************************************************************************

ccccccccc this is NEW OXYGEN continuum
ccccccccc  copied from LBLRTMv5.10
ccccccccc it is indpt of the 79/20 N2/O2 air mixing ratio
       SUBROUTINE CALCONOXYNEW(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY,T,P, 
     $    PARTP,AMT,whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
C 
       INTEGER nptabs
       REAL*8 VB7,VT7,DV7,NPT7 
       REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)  
       REAL*8 O2C213(59),O2C253(59),O2C293(59) 
 
       COMMON /CO2CT/  VB7,VT7,DV7,NPT7,O2C213,O2C253,O2C293 
       COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB 
 
       INTEGER IDGAS, NFREQ, NLAY, whichlayer 
       REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*) 
 
       INTEGER IPT,ILAY 
 
      REAL*8 vb,vt,dv 
      INTEGER npt 
 
      REAL*8 ALL,WO2,XLOSMT,RHOFAC,PAVE,TAVE,XO2CN,P0,T0 
 
      REAL*8 V1C,V2C,DVC,tau_fac 
      INTEGER NPTC,J 
 
      REAL*8 C0(5050),V1,raSFr(5050),gamma 
 
       REAL*8 C2, AVOG, TS, PS 
       DATA C2/1.4387863/ 
       DATA AVOG/6.022045E+26/ 
       DATA TS/296.0/ 
       DATA PS/1.0/ 
 
      DATA XLOSMT / 2.68675E+19 /  
      DATA P0/1013.0/ 
      DATA T0/296.0/ 

  
C-----------------------------------------------------------------------

       DVABS = 1.          
       DVABS = FSTEP
       V1ABS = INT(Freq(1))     
       V1=Freq(1) 
       IF (V1.LT.0.) V1ABS = V1ABS-1. 
       V1ABS = V1ABS-3.*DVABS         
       V2ABS = INT(Freq(NFREQ)+3.*DVABS+0.5)   
       NPTABS = (V2ABS-V1ABS)/DVABS+1.5 
       XO2CN=1.0

       IF (IDGAS .EQ. 7)THEN
         VB = VB7
         VT = VT7
         DV = DV7
         NPT = NPT7
         ENDIF

       DO IPT=1,NFREQ
         CON(IPT)=0.0
         END DO

       IF (IDGAS .EQ. 7) THEN
         DO ILAY=whichlayer,whichlayer
           ALL = AMT(ILAY)*AVOG
           WO2 = ALL
           TAVE=t(iLAY)
           PAVE=partp(ILAY)*p0  

           !! has to be partp since  pV = N k T ==> N/V = p/(k T)
           !! rhofac = (Pave/P0)*(273./Tave) == # of molecules cm-3

           gamma = partp(ilay)/p(ilay)

C     ********    O2 OXYGEN COLLISION INDUCED FUNDAMENTAL  ***********   

c          Only calculate if from 1340 to 1850. cm-1 
           IF (((freq(NFREQ).gt.1340.0).AND.(freq(1).lt.1850.))) THEN
             
             rhofac = (Pave/P0)*(273./Tave) 
             tau_fac = WO2 * 1.e-20 * rhofac  
c 
c            Wk(7) === WO2 is the oxygen column amount in units of molec/cm2 
c            rhofac is in units of amagats (air) 
c            The temperature correction is done in the subroutine o2_ver_1: 
c  
             call o2_ver_1_SSM(raSFr,v1c,v2c,dvc,nptc,c0,tave,gamma) 
             !now spline them onto the finer grid 
             CALL XSPL(raSFr,c0,nptc,Freq,con,nfreq) 

c note by sergio : I have removed radiation field in o2_ver_1 so
c           c0 are the oxygen absorption coefficients at temperature tave  
c              - these absorption coefficients are in units of 
c                   [(cm^2/molec) 10^20)]/(amagat)  
c              - cm-1 in the denominator arises through the removal 
c                   of the radiation field 
c              - for this case, an amagat is interpreted as one 
c                   loshmidt of air (273K) 

            DO  J = 1, NFREQ
              CON(J) = tau_fac * con(J) * XO2CN
              IF ((freq(j) .gt. 1850.0) .or. (freq(j) .lt. 1340.0) .or. 
     $          (con(j) .lt. 0.0)) then
                con(j)=0.0
                end if
              END DO
            endif 
           END DO
         END IF

       RETURN
       END

c************************************************************************
      subroutine o2_ver_1_SSM(raF,v1c,v2c,dvc,nptc,c,T,gamma) 

      REAL*8 V1C,V2C,DVC,T
      REAL*8 C(*),raF(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050),gamma    
      INTEGER NPTABS,NPTS
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      
      REAL*8 V1S,V2S,DVS,xo2(103),xo2t(103) 
      COMMON /o2_f  / V1S,V2S,DVS,NPTS,xo2,xo2t 

      INTEGER I1,I2,I,J
      REAL*8 xktfac,factor,To,xlosmt,vj
 
c     Oxygen Collision Induced Fundamental 
 
c     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann, 
c     and Ch. Boulet 
c        Infrared collision-induced absorption by O2 near 6.4 microns for 
c        atmospheric applications: measurements and emprirical modeling,  
c         Appl. Optics, 36, 563-569, (1997). 
 
      DATA To/ 296./, xlosmt/ 2.68675e+19/ 
 
      xktfac = (1./To)-(1./T) 

c the partp/p is to make everything come out OK in the main routine      
c gamma = partp/p for the current layer
      factor = (1.e+20 /xlosmt) / gamma

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

      do j=1,nptc 
         i = i1+j 
         C(J) = 0.
         raF(j)=V1C+DVC*FLOAT(J-1) 
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         VJ = V1C+DVC*FLOAT(J-1)
c        the radiation field is removed with 1/vj  
c        c(j) = factor * xo2(i)* exp(xo2t(i)*xktfac) / vj 
c sergio takes this away
         c(j) = factor * xo2(i)* exp(xo2t(i)*xktfac)
 10      CONTINUE
         END DO
      
      return 
      end 

c************************************************************************

       SUBROUTINE CALCONNITNEW(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY,T,P, 
     $    PARTP,AMT,whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
C 
       INTEGER NPT22,NPTABS
       REAL*8 VB22,VT22,DV22
       REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)  
 
       COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB 
       COMMON /CN2C0/ VB22,VT22,DV22,NPT22 
 
       INTEGER IDGAS, NFREQ, NLAY, whichlayer 
       REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*) 
 
       INTEGER IPT,ILAY 
 
       REAL*8 TFAC(kMaxLayer) 
 
      REAL*8 vb,vt,dv 
      INTEGER npt 
 
      REAL*8 ALL,WN2,WXN2,XLOSMT,RHOFAC,WTOT,PAVE,TAVE 
      REAL*8 XFRCN,XN2CN,P0,T0 
 
      REAL*8 V1C,V2C,DVC,tau_fac 
      INTEGER NPTC,J 
 
      REAL*8 C0(5050),C1(5050),V1,raSFr(5050),gamma 
      REAL*8 CON0(MaxLen),CON1(MaxLen),conself,confor,efficiency 
 
       REAL*8 C2, AVOG, TS, PS 
       DATA C2/1.4387863/ 
       DATA AVOG/6.022045E+26/ 
       DATA TS/296.0/ 
       DATA PS/1.0/ 
 
      DATA XLOSMT / 2.68675E+19 /  
      DATA P0/1013.0/ 
      DATA T0/296.0/ 
  
C-----------------------------------------------------------------------

       DVABS = 1.          
       DVABS = FSTEP
       V1ABS = INT(Freq(1))     
       V1=Freq(1) 
       IF (V1.LT.0.) V1ABS = V1ABS-1. 
       V1ABS = V1ABS-3.*DVABS         
       V2ABS = INT(Freq(NFREQ)+3.*DVABS+0.5)   
       NPTABS = (V2ABS-V1ABS)/DVABS+1.5 

       XN2CN=1.0
       XFRCN=1.0

       IF (IDGAS .EQ. 22)THEN
         VB = VB22
         VT = VT22
         DV = DV22
         NPT = NPT22
         ENDIF

       DO IPT=1,NFREQ
         CON(IPT)=0.0
         END DO

       IF (IDGAS .EQ. 22) THEN
         DO ILAY=whichlayer,whichlayer
           ALL = AMT(ILAY)*AVOG

C          The following puts WXN2 units in 1./(CM AMAGAT) 
c          Wk(22) is the nitrogen column amount in units of molec/cm2 
c remember DATA AVOG/6.022045E+26/ === molecules per kilomole
c we use GENLN2 amounts === kilomoles/cm2 
c so all=amt(ilay)*avog is in molecules/cm2
           WN2=ALL

           WXN2 = WN2/XLOSMT

           !check to see if next two lines are right
           WTOT=WN2*1e-20            
           PAVE=partp(ILAY)*p0  
           !! has to be partp since  pV = N k T ==> N/V = p/(k T)
           !! rhofac = (Pave/P0)*(273./Tave) == # of molecules cm-3
           !!referenced to that at STP (the Loschmidt number)
           !!p/kT == 101300/1.38e-23/273 = 2.6889e25 m-3 = 2.6889e19 cm-3

           TAVE=t(iLAY)
           TFAC(iLay) = (TAVE-T0)/(220.-T0) 

           gamma = partp(ilay)/p(ilay)

C     ******** NITROGEN COLLISION INDUCED PURE ROTATION BAND  ******** 
C           Model used: 
C           Borysow, A, and L. Frommhold, "Collision-induced 
C             rototranslational absorption spectra of N2-N2 
C            pairs for temperatures from 50 to 300 K", The 
C            Astrophysical Journal, 311, 1043-1057, 1986. 
           IF ((freq(NFREQ).gt.-10.0).AND.(freq(1).lt.350.)) THEN

C            RHOFAC units are AMAGATS 
             RHOFAC = ((WN2*1.e-20)/WTOT)*(PAVE/P0)*(273./TAVE) 

             CALL N2R296 (raSFr,V1C,V2C,DVC,NPTC,C0) 
             CALL N2R220 (raSFr,V1C,V2C,DVC,NPTC,C1) 

           CALL XSPL(raSFr,c0,nptc,Freq,con0,nfreq) 
           CALL XSPL(raSFr,c1,nptc,Freq,con1,nfreq) 

             DO J = 1, NFREQ
               IF (CON0(J).GT.0. .AND. CON1(J).GT.0.) THEN
                 CON(J) = (WXN2*RHOFAC*CON0(J)*
     $                     (CON1(J)/CON0(J))**TFAC(iLay)) * XN2CN 
                 END IF
               END DO
              

C        ********    NITROGEN COLLISION INDUCED FUNDAMENTAL ********
c        version_1 of the Nitrogen Collision Induced Fundamental 
C        Only calculate if V2 > 2085. cm-1 and V1 <  2670. cm-1 
         elseif ((freq(NFREQ).gt.2085.0).and.(freq(1).lt.2670.)) then 
 
           rhofac = (Pave/P0)*(273./Tave) 
           tau_fac = Wn2 * 1.e-20 * rhofac

c           Wn2 is in units of molec/cm2 
c           rhofac is in units of amagats (air) 
c 
c           The temperature correction is done in subroutine n2_ver_1: 
c this will be a very coarse output vector, from v1c to v2c at spacing dvc
c the coarse freqs are in raSFr, with the continuum coeffs in c0
c there are nptc non zero points in these vectors
           call n2_ver_1_pureN2(raSFr,v1c,v2c,dvc,nptc,c0,tave,gamma) 
           !now spline them onto the finer grid 
           CALL XSPL(raSFr,c0,nptc,freq,con,nfreq) 

c          c0 are the PURE nitrogen absorption coefficients at  
c          temperature tave; radiation field has been removed so that  
c              - these absorption coefficients are in units of 
c                   [(cm^2/molec) 10^20)]/( amagat)  
c              - cm-1 in the denominator arises through the removal 
c                   of the radiation field 
c              - for this case, an amagat is interpreted as one 
c                   loshmidt of air (273K) 
c       which is number of molecules at STP no=2.678*10^19 cm-3

c this is the efficiency of n2_o2 vs n2_n2 collisions 

          efficiency = (1.294-0.4545*tave/296.0)               

           DO  J = 1, NFREQ 

             !self contribution to the continuum
             conself = tau_fac * con(J) * (partp(ilay)/p(ilay)) * XN2CN 

             !foreign contribution to the continuum
             confor=tau_fac * con(J) *((p(ilay)-partp(ilay))/p(ilay)) 
     $               * XFRCN * efficiency 

             CON(J)=conself+confor

             IF ((freq(j) .gt. 2670.0) .or. (freq(j) .lt. 2085.0) .or. 
     $          (con(j) .lt. 0.0)) then
               con(j)=0.0
               end if

             END DO

            endif 
         END DO               !end loop over layers
       ENDIF

       RETURN
       END

c************************************************************************
c this computes the PURE n2 abs
      subroutine n2_ver_1_pureN2 (raF,v1c,v2c,dvc,nptc,c,T,gamma)

      REAL*8 V1C,V2C,DVC,T
      REAL*8 C(*),raF(*),gamma               
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,xn2(118),xn2t(118) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /N2_F/ V1S,V2S,DVS,NPTS,xn2,xn2t                         
      
      INTEGER I1,I2,I,J
      REAL*8 xktfac,a1,a2,factor,vj

c     Nitrogen Collision Induced Fundamental

c     Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M. Hartman
c        Infrared collision-induced absorption by N2 near 4.3 microns for
c        atmospheric applications: measurements and emprirical modeling, 
c         Appl. Optics, 35, 5911-5917, (1996).
c
      REAL*8 To,xlosmt,vmr_n2
      DATA  To/ 296./, xlosmt/ 2.68675e+19/, vmr_n2/ 0.78 /

      xktfac = (1./To)-(1./T)
      
      a1 = 0.8387
      a2 = 0.0754


c the partp/p is to make everything come out OK in the main routine      
c gamma = partp/p for the current layer
      factor = (1.e+20 /xlosmt) / gamma

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
 
      do j=1,nptc
         i = i1+j
         C(J) = 0.
         VJ = V1C+DVC*FLOAT(J-1)
         raF(j)=vj
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         VJ = V1C+DVC*FLOAT(J-1)
         raF(j)=vj
c        in orig code the radiation field is removed with 1/vj 
c         c(j) = factor * xn2(i)* exp(xn2t(i)*xktfac) / vj
c sergio takes this away
         c(j) = factor * xn2(i)* exp(xn2t(i)*xktfac)
 10      continue
         end do

      return
      end
c************************************************************************
      SUBROUTINE N2R296 (raF,V1C,V2C,DVC,NPTC,C)
 
C     Model used:
C      Borysow, A, and L. Frommhold, "Collision-induced
C         rototranslational absorption spectra of N2-N2
C         pairs for temperatures from 50 to 300 K", The
C         Astrophysical Journal, 311, 1043-1057, 1986.
 
      REAL*8 V1C,V2C,DVC
      REAL*8 C(*),raF(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(73) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /N2RT0/ V1S,V2S,DVS,NPTS,S                           
      
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
 
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
 
      DO 10 J = 1, NPTC
         I = I1+J
         raF(J) = V1C+DVC*FLOAT(J-1)
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
 10       CONTINUE
 
      RETURN
 
      END
 
c************************************************************************
      SUBROUTINE N2R220 (raF,V1C,V2C,DVC,NPTC,C)
C
C     Model used:
C      Borysow, A, and L. Frommhold, "Collision-induced
C         rototranslational absorption spectra of N2-N2
C         pairs for temperatures from 50 to 300 K", The
C         Astrophysical Journal, 311, 1043-1057, 1986.
C

      REAL*8 V1C,V2C,DVC
      REAL*8 C(*),raF(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      REAL*8 V1S,V2S,DVS,S(73) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /N2RT1/ V1S,V2S,DVS,NPTS,S                           
      
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
 
c*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2 
 
      DO 10 J = 1, NPTC
         I = I1+J
         raF(J) = V1C+DVC*FLOAT(J-1)
         C(J) = 0.
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
         C(J) = S(I)
 10       CONTINUE
 
      RETURN
 
      END
 
c************************************************************************
