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

      double precision dasg0(MaxPQR),dastren0(MaxPQR)
      double precision daGam0(MaxPQR),daGam2(MaxPQR)
      double precision daanuVC(MaxPQR),daeta(MaxPQR),daShift0(MaxPQR),daShift2(MaxPQR)

      double precision Sig0,SigF,dSig,daSg(MaxLen)
      double complex daVt(MaxLen),daVtall(MaxLen)
      integer iI,iJ,iAll

      print *,'input numlines   mass(amu) T(K) : '
      read *,iAll,mass,T
      GamD = 3.581163D-07*dsqrt(T/mass)*sg0/2.0
      
      print *,'input Sig0,SigF,dSig '
      read *,Sig0,SigF,dSig

      print *,'input the following line parameters, separately line by line ...'
      print *,'  [linecenter] cm-1  lineStrength cm2/mol'
      print *,'  [Lorentz SpeedAvgWidth] [SpeedDependenceWidth] cm-1'
      print *,'  [SpeedAvgShift] [SpeedDependanceShift] cm-1'
      print *,'  [VelChangeFreq cm-1] and [CorrelationParam] '

      do iI = 1,iAll
        read *,sg0,stren0,Gam0,Gam2,Shift0,Shift2,anuVC,eta
	dasg0(iI) = sg0
	dastren0(iI) = stren0
	daGam0(iI) = Gam0
	daGam2(iI) = Gam2
	daShift0(iI) = Shift0
	daShift2(iI) = Shift2
	daanuVC(iI) = anuVC
	daeta(iI) = eta
      end do
      
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
        daSg(iI) = Sig0 + (iI-1)*dSig
      END DO


c      CALL loop_pCqSDHC_plus(sg0,stren0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,iNumPts,
c     &    daSg,daVt)
      do iI = 1,MaxLen
        daVtall(iI) = cmplx(0.0d0,0.0d0)
      end do
      do iI = 1,iAll
	sg0 = dasg0(iI)
	stren0 = dastren0(iI)
	Gam0 = daGam0(iI)
	Gam2 = daGam2(iI)
	Shift0 = daShift0(iI)
	Shift2 = daShift2(iI)
	anuVC = daanuVC(iI)
	eta = daeta(iI)      
        CALL loop_pCqSDHC(sg0,stren0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,iNumPts,
     &                    daSg,daVt)
        do iJ = 1,NumPts
           daVtall(iI) = daVtall(iI) + daVt(iI)
        end do       
      END DO
      
      DO iI = 1,iNumPts
        write(*,200) iI,daSg(iI),real(daVtall(iI)),imag(daVtall(iI))
      END DO
      
 200  FORMAT(I10,3(' ',D16.8))
      end
      
