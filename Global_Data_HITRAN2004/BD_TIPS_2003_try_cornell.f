c compile with 
c   f77 -o BD_TIPS_2003_try.x -N109 -W BD_TIPS_2003_try.f

      IMPLICIT NONE

      INTEGER iMol,iISO,iG
      REAL*8  rTemp,rQT

      print *,'Enter iMol iISO  iG  rT : '
      read *,iMol,iISO,iG,rTemp

      CALL BD_TIPS_2003(iMol,rTemp,iISO,iG,rQT)
      print *,rQT
    
      END

c************************************************************************
c this is from /asl/data/hitran/HITRAN04/Global_Data/TIPS/BD_TIPS_2003.for
c this is from http://www.astro.cornell.edu/share/gierasch/akm/linebyline.f

