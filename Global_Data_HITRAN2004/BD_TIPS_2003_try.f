c compile with 
c   f77 -o BD_TIPS_2003_try.x -N109 -W -O2 BD_TIPS_2003_try.f

      implicit DOUBLE PRECISION (a-h,o-z)
c      IMPLICIT NONE

      INTEGER iMol,iISO,iG
      REAL*8  rTemp,rQT

      print *,'Enter iMol iISO  iG  rT : '
      read *,iMol,iISO,iG,rTemp

      CALL BD_TIPS_2003(iMol,rTemp,iISO,iG,rQT)
      print *,rQT
    
      END

c************************************************************************
      include 'BD_TIPS_2003_try_include1.f'
c      include 'BD_TIPS_2003_try_include2.f'

