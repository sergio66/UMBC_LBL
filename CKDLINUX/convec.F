C     computational subroutine
      subroutine convec(x,y,z,nx,ny)
      complex*16 x(*), y(*), z(*)
c      mwsize nx, ny
      integer nx, ny

C     Initialize the output array
      do 10 i=1,nx+ny-1
         z(i) = (0.0,0.0)
 10   continue

      do 30 i=1,nx
         do 20 j=1,ny
            z(i+j-1) = z(i+j-1) + x(i) * y(j)
 20      continue
 30   continue
      return
      end
