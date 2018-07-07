
	subroutine 	qSDHC(sg0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,
     &sg,LS_qSDHC_R,LS_qSDHC_I)
C-------------------------------------------------
C	"qSDHC": quadratic-Speed-Dependent Hard-Collision
C	Subroutine to Compute the complex normalized spectral shape of an 
C	isolated line by the qSDHC model
C
C	Input/Output Parameters of Routine (Arguments or Common)
C	---------------------------------
C	T	    : Temperature in Kelvin (Input).
C	amM1	: Molar mass of the absorber in g/mol(Input).
C	sg0		: Unperturbed line position in cm-1 (Input).
C     GamD	: Doppler HWHM in cm-1 (Input)
C	Gam0	: Speed-averaged line-width in cm-1 (Input). 	
C	Gam2	: Speed dependence of the line-width in cm-1 (Input).
C	anuVC	: Velocity-changing frequency in cm-1 (Input).
C	Shift0	: Speed-averaged line-shift in cm-1 (Input).
C	Shift2	: Speed dependence of the line-shift in cm-1 (Input)	 
C	sg		: Current WaveNumber of the Computation in cm-1 (Input).
C
C	Output Quantities (through Common Statements)
C	-----------------
C	LS_qSDHC_R: Real part of the normalized spectral shape (cm)
C	LS_qSDHC_I: Imaginary part of the normalized spectral shape (cm)
C
C	Called Routines: 'CPF'	(Complex Probability Function)
C	---------------  'CPF3'	(Complex Probability Function for the region 3)
C
C	Called By: Main Program
C	---------
C
C     Double Precision Version
C
C-------------------------------------------------
	implicit none
	 double precision sg0,GamD
	 double precision Gam0,Gam2,anuVC,Shift0,Shift2
	 double precision sg
	 double precision pi,rpi,cte
	 double precision xz1,xz2,yz1,yz2,xXb,yXb
	 double precision wr1,wi1,wr2,wi2,wrb,wib
	 double precision SZ1,SZ2,DSZ,SZmx,SZmn
	 double precision LS_qSDHC_R,LS_qSDHC_I
	double complex c0,c2,c0t,c2t
	double complex X,Y,iz,Z1,Z2
	double complex Aterm,LS_qSDHC
C
C-------------------------------------------------
C
	cte=dsqrt(dlog(2.D0))/GamD
	pi=4.d0*datan(1.d0)
	rpi=dsqrt(pi)
	iz=dcmplx(0.d0,1.d0)
c Calculating the different parameters 
	c0=dcmplx(Gam0,-Shift0)
	c2=dcmplx(Gam2,-Shift2)
	c0t=(c0-1.5d0*c2)+anuVC
	c2t=c2
	Y=1.d0/((2.d0*cte*C2t))**2			
C
	
	X=(iz*(sg-sg0)+c0t)/c2t
c	
	if (cdabs(C2t).eq.0.d0) go to 110
	if (cdabs(X).le.3.d-8*cdabs(Y)) go to 120
	if (cdabs(Y).le.1.d-15*cdabs(X)) go to 140
c calculating Z1 and Z2
	Z1=cdsqrt(X+Y)-cdsqrt(Y)
	Z2=Z1+2.d0*cdsqrt(Y)
c calculating the real and imaginary parts of Z1 and Z2
	xZ1=-dimag(Z1)
	yZ1=dreal(Z1)
	xZ2=-dimag(Z2)
	yZ2=dreal(Z2)
c check if Z1 and Z2 are close to each other
	SZ1=dsqrt(xZ1*xZ1+yZ1*yZ1)
	SZ2=dsqrt(xZ2*xZ2+yZ2*yZ2)
	DSZ=dabs(SZ1-SZ2)
	SZmx=dmax1(SZ1,SZ2)
	SZmn=dmin1(SZ1,SZ2)
c when Z1 and Z2 are close to each other, ensure that they are in 
c the same interval of CPF 
	if (DSZ.le.1.d0.and.SZmx.gt.8.d0.and.SZmn.le.8.d0) then
	Call CPF3 ( xZ1, yZ1, WR1, WI1 ) 
	Call CPF3 ( xZ2, yZ2, WR2, WI2 ) 
	else	
	Call CPF ( xZ1, yZ1, WR1, WI1 ) 
	Call CPF ( xZ2, yZ2, WR2, WI2 ) 
	endif
c calculating the A term of the profile
	Aterm=rpi*cte*(dcmplx(wr1,wi1)-dcmplx(wr2,wi2))
	go to 10
c when C2t=0
110   continue
	Z1=(iz*(sg-sg0)+C0t)*cte
	xZ1=-dimag(Z1)
	yZ1=dreal(Z1)
	Call CPF ( xZ1, yZ1, WR1, WI1 )
	Aterm=rpi*cte*dcmplx(WR1,WI1)
	go to 10
c when abs(Y) is much larger than abs(X)
120   continue
	Z1=(iz*(sg-sg0)+C0t)*cte
	Z2=cdsqrt(X+Y)+cdsqrt(Y)
	xZ1=-dimag(z1)
	yZ1=dreal(z1)
	xZ2=-dimag(z2)
	yZ2=dreal(z2)
	Call CPF ( xZ1, yZ1, WR1, WI1 )
	Call CPF ( xZ2, yZ2, WR2, WI2 ) 
	Aterm=rpi*cte*(dcmplx(WR1,WI1)-dcmplx(WR2,WI2))
	go to 10
c when abs(X) is much larger than abs(Y)
140   continue
	if (cdabs(cdsqrt(X)).le.4.d3) then
	  xXb=-dimag(cdsqrt(X))
	  yXb=dreal(cdsqrt(X))
	  Call CPF ( xXb, yXb, WRb, WIb ) 
	  Aterm=(2.d0*rpi/C2t)*(1.d0/rpi-cdsqrt(X)*dcmplx(WRb,WIb))
cc and when abs(X) is much larger than 1
	else
	  Aterm=(1.d0/C2t)*(1.d0/X-1.5d0/(X**2))
	endif
c
10    continue
c
	LS_qSDHC=(1.d0/pi)*Aterm/
     ,(1.d0-anuVC*Aterm)

	LS_qSDHC_R=dreal(LS_qSDHC)
	LS_qSDHC_I=dimag(LS_qSDHC)

   
      Return
      End Subroutine qSDHC