c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

ccccccccc this is NEW OXYGEN continuum
ccccccccc  copied from LBLRTMv5.10
ccccccccc it is tied to the 79/20 N2/O2 air mixing ratio

c************************************************************************
ccccccccccc include the common data blocks
c************************************************************************
      include 'brandnew_o2ctt.f'
c************************************************************************
      include 'brandnew_n2ct0.f'
c************************************************************************

ccccccccc it is tied to the 79/20 N2/O2 air mixing ratio
       SUBROUTINE CALCONOXYNEW_LBLRTM(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY,
     $    T,P,PARTP,AMT,whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
       REAL*8 VB1,VT1,DV1,NPT1,VB7,VT7,DV7,NPT7
       REAL*8 V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050) 
       REAL*8 VB7,VT7,DV7,NPT7,O2C213(59),O2C253(59),O2C293(59)

       COMMON /CO2CT/  VB7,VT7,DV7,NPT7,O2C213,O2C253,O2C293
       COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB

       REAL*8 X

       INTEGER IDGAS, NFREQ, NLAY, whichlayer
       REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*)

C      Variables for new water continuum adjustment terms
       INTEGER JFAC
       REAL*8 ALPHS2, BETAS, V0S, FACTRS, HWSQF, BETAF, V0F, FACTRF,
     $    V0F2, HWSQF2, BETA2, SFAC, alpha2, scor, FSCAL

       INTEGER IPT,ILAY

       REAL*8 C2, AVOG, TS, PS, PO
       DATA C2/1.4387863/
       DATA AVOG/6.022045E+26/
       DATA TS/296.0/
       DATA PS/1.0/

       REAL*8 TFAC(kMaxLayer),XAMT(kMaxLayer)

      REAL*8 vb,vt,dv,vs2,vf2,vf6,vf4
      INTEGER iL,npt,iF

      REAL*8 ALL,FF,WO2,WXN2,XLOSMT,RHOFAC,WTOT,PAVE,TAVE,XO2CN,P0,T0

      DATA XLOSMT / 2.68675E+19 / 
      DATA P0/1013.0/
      DATA T0/296.0/

      REAL*8 V1C,V2C,DVC,VJ,tau_fac
      INTEGER NPTC,J

      REAL*8 C0(5050),V1,raSSFREQ(5050)
  
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


C     ********    O2 OXYGEN COLLISION INDUCED FUNDAMENTAL  ***********   
c 
c     version_1 of the Oxygen Collision Induced Fundamental  
c     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann, 
c        and Ch. Boulet, 
c        Infrared collision-induced absorption by O2 near 6.4 microns for 
c        atmospheric applications: measurements and emprirical modeling,  
c        Appl. Optics, 35, 5911-5917, (1996). 
 
c          Only calculate if V2 > 1340. cm-1 and V1 <  1850. cm-1 
           IF (((freq(NFREQ).gt.1340.0).AND.(freq(1).lt.1850.))) THEN
             
             rhofac = (Pave/P0)*(273./Tave) 
             tau_fac = WO2 * 1.e-20 * rhofac  
c 
c            Wk(7) === WO2 is the oxygen column amount in units of molec/cm2 
c            rhofac is in units of amagats (air) 
c  
c            The temperature correction is done in the subroutine o2_ver_1: 
c  
             call o2_ver_1 (raSSFREQ,v1c,v2c,dvc,nptc,c0,tave) 
             !now spline them onto the finer grid 
             CALL XSPL(raSSFreq,c0,nptc,Freq,con,nfreq) 

c           c0 are the oxygen absorption coefficients at temperature tave  
c              - these absorption coefficients are in units of 
c                   [(cm^2/molec) 10^20)]/(cm-1  amagat)  
c              - cm-1 in the denominator arises through the removal 
c                   of the radiation field 
c              - for this case, an amagat is interpreted as one 
c                   loshmidt of air (273K) 
ccccccc note by sergio : I have removed radiation field in n2_ver_1 so
cccccccc     these absorption coefficients are in units of 
cccc                   [(cm^2/molec) 10^20)]/(amagat)  

            DO  J = 1, NFREQ
              IF ((freq(j) .gt. 1850.0) .or. (freq(j) .lt. 1340.0)) then
                con(j)=0.0
              ELSE
                CON(J) = tau_fac * con(J) * XO2CN
                end if
              IF (CON(J) .lt. 0) CON(J)= 0.0
              END DO
            endif 
           END DO
         END IF

       RETURN
       END

c************************************************************************
      subroutine o2_ver_1 (raF,v1c,v2c,dvc,nptc,c,T) 

      REAL*8 V1C,V2C,DVC,T
      REAL*8 C(*),raF(*)                                                   
      INTEGER NPTC
    
      REAL*8 V1ABS,V2ABS,DVABS,ABSRB(5050)             
      INTEGER NPTABS,NPTS
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      
      REAL*8 V1S,V2S,DVS,xo2(103),xo2t(103) 
      COMMON /o2_f  / V1S,V2S,DVS,NPTS,xo2,xo2t 

      INTEGER I1,I2,I,J
      REAL*8 xktfac,factor,To,xlosmt,vj
 
c     Oxygen Collision Induced Fundamental 
 
c     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann, 
c                                                         and Ch. Boulet 
c        Infrared collision-induced absorption by O2 near 6.4 microns for 
c        atmospheric applications: measurements and emprirical modeling,  
c         Appl. Optics, 35, 5911-5917, (1996). 
 
      DATA To/ 296./, xlosmt/ 2.68675e+19/ 
 
      xktfac = (1./To)-(1./T) 
      
c     correct formulation for consistency with LBLRTM:  
      factor = (1.e+20 /xlosmt)  

c     A factor of 0.21, the mixing ration of oxygen, in the Thibault et al. 
c     formulation is not included here.  This factor is in the column amount. 
c  sergio puts this back in
      factor = (1.e+20 /xlosmt) / 0.21 

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

ccccccccc this is NEW n2 continuum
ccccccccc  copied from LBLRTMv5.10
ccccccccc it is tied to the 79/20 N2/O2 air mixing ratio
       SUBROUTINE CALCONNITNEW_LBLRTM(CON,IDGAS,NFREQ,FREQ,FSTEP,
     $    NLAY,T,P,PARTP,AMT,whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
       REAL*8 VB1,VT1,DV1,NPT1,VB22,VT22,DV22,NPT22
       REAL*8 V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050) 

       COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
       COMMON /CN2C0/ VB22,VT22,DV22,NPT22

       REAL*8 X

       INTEGER IDGAS, NFREQ, NLAY, whichlayer
       REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*)

C      Variables for new water continuum adjustment terms
       INTEGER JFAC
       REAL*8 ALPHS2, BETAS, V0S, FACTRS, HWSQF, BETAF, V0F, FACTRF,
     $    V0F2, HWSQF2, BETA2, SFAC, alpha2, scor, FSCAL

       INTEGER IPT,ILAY

       REAL*8 C2, AVOG, TS, PS, PO
       DATA C2/1.4387863/
       DATA AVOG/6.022045E+26/
       DATA TS/296.0/
       DATA PS/1.0/

       REAL*8 TFAC(kMaxLayer),XAMT(kMaxLayer)

      REAL*8 vb,vt,dv,vs2,vf2,vf6,vf4
      INTEGER iL,npt,iF

      REAL*8 ALL,FF,WN2,WXN2,XLOSMT,RHOFAC,WTOT,PAVE,TAVE,XN2CN,P0,T0

      DATA XLOSMT / 2.68675E+19 / 
      DATA P0/1013.0/
      DATA T0/296.0/

      REAL*8 V1C,V2C,DVC,VJ,tau_fac
      INTEGER NPTC,J

      REAL*8 C0(5050),C1(5050),V1,raSSFREQ(5050)
      REAL*8 CON0(MaxLen),CON1(MaxLen)
  
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


C     ******** NITROGEN COLLISION INDUCED PURE ROTATION BAND  ******** 
C           Model used: 
C           Borysow, A, and L. Frommhold, "Collision-induced 
C             rototranslational absorption spectra of N2-N2 
C            pairs for temperatures from 50 to 300 K", The 
C            Astrophysical Journal, 311, 1043-1057, 1986. 
           IF ((freq(NFREQ).gt.-10.0).AND.(freq(1).lt.350.)) THEN

C            RHOFAC units are AMAGATS 
             RHOFAC = ((WN2*1.e-20)/WTOT)*(PAVE/P0)*(273./TAVE) 

             CALL N2R296 (raSSFReq,V1C,V2C,DVC,NPTC,C0) 
             CALL N2R220 (raSSFReq,V1C,V2C,DVC,NPTC,C1) 

           CALL XSPL(raSSFreq,c0,nptc,Freq,con0,nfreq) 
           CALL XSPL(raSSFreq,c1,nptc,Freq,con1,nfreq) 

             DO J = 1, NFREQ
               IF (CON0(J).GT.0. .AND. CON1(J).GT.0.) THEN
                 CON(J) = (WXN2*RHOFAC*CON0(J)*
     $                     (CON1(J)/CON0(J))**TFAC(iLay)) * XN2CN 
                 END IF
               END DO
              

C        ********    NITROGEN COLLISION INDUCED FUNDAMENTAL ********
c        version_1 of the Nitrogen Collision Induced Fundamental 
c 
c        Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and  
c        J._M. Hartmann, Infrared collision-induced absorption by  
c        N2 near 4.3 microns for atmospheric applications:  
c        Measurements and emprirical modeling, Appl. Optics, 35,  
c        5911-5917, (1996). 
C        Only calculate if V2 > 2085. cm-1 and V1 <  2670. cm-1 

         elseif ((freq(NFREQ).gt.2085.0).and.(freq(1).lt.2670.)) then 
 
           rhofac = (Pave/P0)*(273./Tave) 
           tau_fac = Wn2 * 1.e-20 * rhofac

c           Wn2 is in units of molec/cm2 
c           rhofac is in units of amagats (air) 
c 
c           The temperature correction is done in subroutine n2_ver_1: 
c this will be a very coarse output vector, from v1c to v2c at spacing dvc
c the coarse freqs are in raSSFreq, with the continuum coeffs in c0
c there are nptc non zero points in these vectors
           call n2_ver_1 (raSSFreq,v1c,v2c,dvc,nptc,c0,tave) 
           !now spline them onto the finer grid 
           CALL XSPL(raSSFreq,c0,nptc,freq,con,nfreq) 
c           CALL XLINEAR(raSSFreq,c0,nptc,freq,con,nfreq) 

c          c0 are the nitrogen absorption coefficients at  
c          temperature tave  
c              - these absorption coefficients are in units of 
c                   [(cm^2/molec) 10^20)]/(cm-1  amagat)  
c              - cm-1 in the denominator arises through the removal 
c                   of the radiation field 
c              - for this case, an amagat is interpreted as one 
c                   loshmidt of air (273K) 
c       which is number of molecules at STP no=2.678*10^19 cm-3
ccccccc note by sergio : I have removed radiation field in n2_ver_1 so
cccccccc     these absorption coefficients are in units of 
cccc                   [(cm^2/molec) 10^20)]/(amagat)  

           DO  J = 1, NFREQ 
             IF ((freq(j) .gt. 2670.0) .or. (freq(j) .lt. 2085.0)) then
               con(j)=0.0
             ELSE
               CON(J) = tau_fac * con(J) * XN2CN
               end if
             IF (CON(J) .lt. 0) CON(J)= 0.0
             END DO

            endif 
         END DO               !end loop over layers
       ENDIF

       RETURN
       END

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
 
      subroutine n2_ver_1 (raF,v1c,v2c,dvc,nptc,c,T)

      REAL*8 V1C,V2C,DVC,T
      REAL*8 C(*),raF(*)                                                   
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
 
c     correct formulation for consistency with LBLRTM:
c     Lafferty et al. reference  assumes that the
c     column amount is that for air 
       factor = (1.e+20 /xlosmt) * (1./vmr_n2) * (a1-a2*(T/To))
c sergio puts this factor squared!
       factor = (1.e+20 /xlosmt) * (a1-a2*(T/To)) /(vmr_n2*vmr_n2)
                            
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
