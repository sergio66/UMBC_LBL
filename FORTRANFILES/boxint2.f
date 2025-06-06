      subroutine boxint2(z,y,nbox,zlen)
c     subroutine boxint2(z,y,nbox,zlen) does a box car integration
c     note that the array sizes are fixed
      
      real*8 z(zlen),y(nbox*zlen)
      integer nbox,zlen
    
      integer jj,ii,is,ie

      do 30 jj=1,zlen
        is = 1   + (jj-1)*nbox
        ie = nbox+ (jj-1)*nbox
        z(jj)=0.0
        do 40 ii=is,ie
          z(jj)=z(jj)+y(ii)
 40       CONTINUE
 30     CONTINUE

      if (nbox .gt. 1) then
        do 50 jj=1,zlen
          z(jj)=z(jj)/nbox
 50       CONTINUE
        end if

      return
      end
