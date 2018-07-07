c compile with 
c f77 -o BD_TIPS_2003_allisotopes.x -N109 -W -O3 BD_TIPS_2003_allisotopes.f

c      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT NONE

      INTEGER iMol,iG,iI,iNumISO
      REAL*8  rTemp,rQT

      print *,'%% Enter iMol iNumISO  rT : '
      read *,iMol,iNumISO,rTemp

      iG = 1   !!!! looking at the HITRAN code, this nuclear degeneracy 
               !!!! is automatically included in the source code

      rQT = 100.0d0
      do iI = 1,iNumISO
        iG = 1     !!! some obscure thing makes this reset to "0" after each
                   !!! call to BD_TIPS_2003
        CALL BD_TIPS_2003(iMol,rTemp,iI,iG,rQT)
        print *,rQT
        end do
  
      END

c************************************************************************
      include 'BD_TIPS_2003_try_include1.f'
c      include 'BD_TIPS_2003_try_include2.f'

