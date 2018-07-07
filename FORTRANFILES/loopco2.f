      subroutine loop(outvect,ziso,mass_iso,brd,strength, 
     $        centerfreq,calcfreq,tempr,numlines,sizewave,numisotopes,
     $        chichi,linetype,pfor,pco2)

      integer numlines,sizewave,numisotopes,chichi,linetype
      real*8 ziso(numlines),mass_iso(numisotopes),brd(numlines)
      real*8 strength(numlines),outvect(sizewave)
      real*8 centerfreq(numlines),calcfreq(sizewave),tempr,pfor,pco2

      real*8 y(sizewave),yi(sizewave),z(sizewave),mass
      integer ii,kk

      do ii=1,sizewave
        outvect(ii)=0.0
        y(ii)=0.0
        yi(ii)=0.0
        end do

c      print *,'computing lineshape for ',numlines,' lines'    
c      print *,'chichi= ',chichi,' linetype = ',linetype
      do kk=1,numlines
        !compute the line shape
        mass=mass_iso(ziso(kk))

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
      include 'vhh2RI.f'
      include 'lorentz_vhh.f'
      include 'vhh.f'
      include 'vhh1.f'
      include 'cousin1.f'
