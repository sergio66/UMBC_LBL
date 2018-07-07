      Subroutine CPF3(X,Y,WR,WI)
C-------------------------------------------------
C "CPF": Complex Probability Function
C .........................................................
C         .       Subroutine to Compute the Complex       .
C         .        Probability Function W(z=X+iY)         .
C         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
C         .    Which Appears when Convoluting a Complex   .
C         .     Lorentzian Profile by a Gaussian Shape    .
C         .................................................
C
C             WR : Real Part of W(z)
C             WI : Imaginary Part of W(z)
C
C This Routine takes into account the region 3 only, i.e. when sqrt(x**2+y**2)>8. 
C
C Accessed Files:  None
C --------------
C
C Called Routines: None                               
C ---------------                                 
C
C Called By: 'pCqSDHC'
C ---------
C
C Double Precision Version
C 
C-------------------------------------------------
C      
      Implicit None
        Integer I
	double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision TT(15),pipwoeronehalf
C      
	Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
	data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,
     ,        9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
	data pipwoeronehalf/0.564189583547756d0/

C Region 3
	zm1=zone/dcmplx(x,y)
	zm2=zm1*zm1
	zsum=zone
	zterm=zone
	do i=1,15
	zterm=zterm*zm2*tt(i)
	zsum=zsum+zterm
	end do
	zsum=zsum*zi*zm1*pipwoeronehalf
	wr=dreal(zsum)
	wi=dimag(zsum)
	return
      End