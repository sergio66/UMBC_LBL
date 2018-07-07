c compile with 
c f77 -o TIPS_2012_allisotopes.x -N109 -W -O3 TIPS_2012_allisotopes.f or
c /asl/opt/absoft/absoft10.0/bin/af77 -o TIPS_2012_allisotopes.x -N109 -W -O3 TIPS_2012_allisotopes.f

c      implicit DOUBLE PRECISION (a-h,o-z)
      IMPLICIT NONE

      INTEGER iMol,iG,iI,iNumISO
      DOUBLE PRECISION rTemp,rQT
      REAL*8 rrQT

      print *,'%% Enter iMol iNumISO  rT : '
      read *,iMol,iNumISO,rTemp

      iG = 1   !!!! looking at the HITRAN code, this nuclear degeneracy 
               !!!! is automatically included in the source code

      !!!!!!! >>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<<<<<<< !!!!!!!
      !!!!!!! note we modified TIPS_2012_try_include1.f to account for
      !!!!!!! gasIDs 39.40.41.42 as those did not have qtips fcns
      !!!!!!! while mass.dat has qtips(296) for those molecules
      !!!!!!! so now qtips(T) = qtips(296) for those molecules
      rQT = 100.0d0
      do iI = 1,iNumISO
        iG = 1     !!! some obscure thing makes this reset to "0" after each
                   !!! call to BD_TIPS_2003
        CALL TIPS_2012(iMol,rTemp,iI,iG*1.d0,rQT)
        rrQT = rQT
        print *,rrQT
        end do
      !!!!!!! note we modified TIPS_2012_try_include1.f to account for
      !!!!!!! gasIDs 39.40.41.42 as those did not have qtips fcns
      !!!!!!! while mass.dat has qtips(296) for those molecules
      !!!!!!! so now qtips(T) = qtips(296) for those molecules
      !!!!!!! >>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<<<<<<< !!!!!!!
  
      END

c************************************************************************
      include 'TIPS_2012_try_include1.f'

