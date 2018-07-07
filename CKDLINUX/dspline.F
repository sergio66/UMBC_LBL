c************************************************************************
       SUBROUTINE xlinear(XA,YA,N,XOUT,YOUT,NOUT)

       include '/home/sergio/SPECTRA/FORTRANLINUX/max.inc'

C real*8 version
C      -----------------------------------------------------------------
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
C      XA  : I  : DOUB arr : x array(N) in increasing order  IN
C      YA  : I  : DOUB arr : y array(N)                      IN
C      N   : I  : INT      : number of points in arrays
C      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
C      YOUT   : O  : DOUB ARR     : y points from spline interpolation
C      NOUT   : I  : INT          : number of points at which to spline
C      -----------------------------------------------------------------
C
C      Parameters
       REAL*8 XA(*),YA(*),XOUT(*),YOUT(*)
       INTEGER N,NOUT

       INTEGER I

       IF (NOUT .GT. Maxlen) THEN
         print *,'in xlinear, can only have Maxlen points! '
         print *,' use smaller input arrays!!!'         
         STOP
         END IF
       
       DO I=1,NOUT
         CALL LINEAR(XA,YA,N,XOUT(I),YOUT(I)) 
         END DO

       RETURN
       END

c************************************************************************
       SUBROUTINE xspl(XA,YA,N,XOUT,YOUT,NOUT)

       include '/home/sergio/SPECTRA/FORTRANLINUX/max.inc'

C real*8 version
C      -----------------------------------------------------------------
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
C      XA  : I  : DOUB arr : x array(N) in increasing order  IN
C      YA  : I  : DOUB arr : y array(N)                      IN
C      N   : I  : INT      : number of points in arrays
C      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
C      YOUT   : O  : DOUB ARR     : y points from spline interpolation
C      NOUT   : I  : INT          : number of points at which to spline
C      -----------------------------------------------------------------
C
C      Parameters
       REAL*8 XA(*),YA(*),XOUT(*),YOUT(*)
       INTEGER N,NOUT

       REAL*8 Y2A(MaxLen),dyp1,dyp2

       INTEGER I

       IF (NOUT .GT. Maxlen) THEN
         print *,'in spline, can only have Maxlen points! '
         print *,' use smaller input arrays!!!'         
         STOP
         END IF
       
       dyp1=1.0e16
       dyp2=1.0e16

       CALL dSPLY2(XA,YA,N,dyp1,dyp2,Y2A)
       DO I=1,NOUT
         CALL dSPLIN(XA,YA,Y2A,N,XOUT(I),YOUT(I)) 
         END DO

       RETURN
       END

c************************************************************************
C  double precision
       SUBROUTINE LINEAR(XA,YA,N,X,Y)

C linear interpolation
C double precision version
C      -----------------------------------------------------------------
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       REAL*8 XA(*),YA(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       REAL*8 A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of points X falls (bisect loop)
       KLO=1
       KHI=N
 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          WRITE(*,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF

       A=(XA(KHI) - X)/H
       B=YA(KHI)-YA(KLO)

       Y=YA(KHI)-A*B

       RETURN
       END

c************************************************************************
       SUBROUTINE dSPLIN(XA,YA,Y2A,N,X,Y)

       include '/home/sergio/SPECTRA/FORTRANLINUX/max.inc'

C real*8 version
C      -----------------------------------------------------------------
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      Y2A : I  : DOUB arr : 2nd derivative of points
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       REAL*8 XA(*),YA(*),Y2A(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       REAL*8 A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of pints X falls (bisect loop)
       KLO=1
       KHI=N
 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          print *,KLO,KHI,XA(KLO),XA(KHI)
 1010     FORMAT('ERROR! dSPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF
C
       A=(XA(KHI) - X)/H
       B=(X - XA(KLO))/H
C
       Y=A*YA(KLO) + B*YA(KHI) + ( Y2A(KLO)*(A**3 - A) +
     $    Y2A(KHI)*(B**3 - B) )*(H**2)/6.0
C
       RETURN
       END

c************************************************************************
C
       SUBROUTINE dSPLY2(XA,YA,N,DYP1,DYPN,Y2A)

       include '/home/sergio/SPECTRA/FORTRANLINUX/max.inc'

C
C real*8 version
C      -----------------------------------------------------------------
C      Calc 2nd derivative as preperation for SPLINT spline routine
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      Y2A : O  : DOUB arr : 2nd derivative array(N)
C      -----------------------------------------------------------------
C
C      Parameters
       REAL*8 XA(*),YA(*),Y2A(*),DYP1,DYPN
       INTEGER N
C
C      Local Variables
       INTEGER I,K
       REAL*8 P,QN,SIG,UN,WORK(MaxLen)
C      -----------------------------------------------------------------
C
C      Lower boundary
       IF (DYP1 .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          Y2A(1)=0.0
          WORK(1)=0.0
       ELSE
C         Set to a specific first derivative
          Y2A(1)=-0.5
          WORK(1)=( 3.0/(XA(2) - XA(1)) )*( (YA(2) - YA(1))/
     $       (XA(2) - XA(1)) - DYP1)
       ENDIF
C
C      Decomposition loop of the tridiagonal algorithm
       DO I=2,N-1
c this is from the progas code
c          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I))
          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
          P=SIG*Y2A(I-1) + 2.0
          Y2A(I)=(SIG - 1.0)/P
          WORK(I)=(YA(I+1) - YA(I))/(XA(I+1) - XA(I)) -
     $       (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
          WORK(I)=( 6.0*WORK(I)/(XA(I+1) - XA(I-1)) -
     $       SIG*WORK(I-1) )/P
       ENDDO
C
C      Upper boundary
       IF (DYPN .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          QN=0.0
          UN=0.0
       ELSE
C         Set to a specific first derivative
          QN=0.5
          UN=( 3.0/(XA(N) - XA(N-1)) )*( DYPN -
     $       (YA(N) - YA(N-1))/(XA(N) - XA(N-1)) )
       ENDIF
       Y2A(N)=(UN - QN*WORK(N-1))/(QN*Y2A(N-1) + 1.0)
C
C      Assign the other 2nd derivatives using the back-substitution
C      loop of the tridiagonal algorithm
       DO K=N-1,1,-1
          Y2A(K)=Y2A(K)*Y2A(K+1) + WORK(K)
       ENDDO
C
       RETURN
       END

c************************************************************************
