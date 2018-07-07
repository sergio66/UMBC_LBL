      include '/home/sergio/SPECTRA/FORTRANLINUX/max.inc'      

C	Input/Output Parameters of Routine (Arguments or Common)
C	---------------------------------
C	T	: Temperature in Kelvin (Input).
C	amM1	: Molar mass of the absorber in g/mol(Input).
C	sg0	: Unperturbed line position in cm-1 (Input).
C       GamD	: Doppler HWHM in cm-1 (Input)
C	Gam0	: Speed-averaged line-width in cm-1 (Input). 	
C	Gam2	: Speed dependence of the line-width in cm-1 (Input).
C	anuVC	: Velocity-changing frequency in cm-1 (Input).
C	eta		: Correlation parameter, No unit (Input).
C	Shift0	: Speed-averaged line-shift in cm-1 (Input).
C	Shift2	: Speed dependence of the line-shift in cm-1 (Input)	 
C	raSg	: Current WaveNumber of the Computation in cm-1 (Input).
c       iNumPts : Number of points to loop over

      double precision sg0,GamD,mass,T,stren0
      double precision Gam0,Gam2,anuVC,eta,Shift0,Shift2
      double precision Sig0,SigF,dSig,raSg(MaxLen)
      double complex raVT(MaxLen)
      integer iI
      
      print *,'input ([linecenter] cm-1  lineStrength cm2/mol)'
      read *,sg0,stren0
      print *,'input ([mass(amu) T(K)])'
      read *,mass,T
      GamD = 3.581163D-07*dsqrt(T/mass)*sg0/2.0
      print *,'input ([SpeedAvgWidth] [SpeedDependenceWidth]) cm-1'
      read *,Gam0,Gam2
      print *,'input ([SpeedAvgShift] [SpeedDependanceShift]) cm-1'
      read *,Shift0,Shift2
      print *,'input ([VelChangeFreq cm-1] and [CorrelationParam]) '
      read *,anuVC,eta

      print *,'input Sig0,SigF,dSig '
      read *,Sig0,SigF,dSig
      print *,'Sig0,SigF,dSig before : ',Sig0,SigF,dSig

c remember 5 point boxcar
      Sig0 = Sig0 - 2*dSig
      SigF = SigF - 5*dSig + 2*dSig
      iNumPts = ceiling(SigF-Sig0)/dSig+1
      print *,'Sig0,SigF,dSig after : ',Sig0,SigF,dSig
      print *,'iNumPts = ',iNumPts
      
      IF (iNumPts .GT. MaxLen) THEN
        print *,'Error iNumPts > MaxLen ',iNumPts,MaxLen
	STOP
      END IF
      DO iI = 1,iNumPts
        raSg(iI) = Sig0 + (iI-1)*dSig
      END DO
      
      CALL loop_pCqSDHC_plus(sg0,stren0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,iNumPts,
     &    raSg,raVt)
      CALL loop_pCqSDHC(sg0,stren0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,iNumPts,
     &    raSg,raVt)

      DO iI = 1,iNumPts
        write(*,200) iI,raSg(iI),real(raVt(iI)),imag(raVt(iI))
      END DO
      
 200  FORMAT(I10,3(' ',D16.8))
      end
      
