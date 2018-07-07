      subroutine loopco2(outvect,ziso,mass_iso,brd,strength, 
     $        centerfreq,calcfreq,tempr,numlines,sizewave,numisotopes,
     $        chichi,linetype,pfor,pco2)

      include 'max.inc'

      integer numlines,sizewave,numisotopes,chichi,linetype
      real*8 ziso(MaxLen),mass_iso(MaxLen),brd(MaxLen)
      real*8 strength(MaxLen),outvect(MaxLen)
      real*8 centerfreq(MaxLen),calcfreq(MaxLen),tempr,pfor,pco2

      real*8 y(MaxLen),yi(MaxLen),z(MaxLen),mass
      integer ii,kk

c      IF (n_in .gt. MaxLen) THEN
c        print *,'in loopco2.f, n_in .gt. MaxLen',n_in,MaxLen
c        STOP
c        END IF

      do ii=1,sizewave
        outvect(ii)=0.0
        y(ii)=0.0
        yi(ii)=0.0
        end do

c      print *,'computing lineshape for ',numlines,' lines'    
c      print *,'chichi= ',chichi,' linetype = ',linetype
      do kk=1,numlines
        !compute the line shape
        mass=mass_iso(idint(ziso(kk)))

        if (linetype .gt. 0) then
          call vhh2RI(y,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
        else
          call lorentz_vhh(y,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
          end if

        if (chichi .lt. 0) then
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+y(ii)*strength(kk)
            end do
        else
          call cousin(z,calcfreq,centerfreq(kk),brd(kk),tempr,
     $                pfor,pco2,sizewave)                       
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+y(ii)*z(ii)*strength(kk)
            end do
          end if
          
        end do

      return
      end

c ***********************************************************************
      include 'vhh2RI.F'
      include 'lorentz_vhh.F'
      include 'vhh.F'
      include 'vhh1.F'
      include 'cousin1.F'
