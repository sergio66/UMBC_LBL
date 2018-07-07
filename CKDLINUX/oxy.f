c************************************************************************
      include 'brandnew_o2ctt.f'
c************************************************************************
ccccccccc this is NEW OXYGEN continuum
ccccccccc  copied from LBLRTMv5.10
ccccccccc it is indpt of the 79/20 N2/O2 air mixing ratio
       SUBROUTINE CALCONOXYNEW(CON,IDGAS,NFREQ,FREQ,FSTEP,NLAY,T,P, 
     $    PARTP,AMT,whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
       INTEGER NPTABS
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






