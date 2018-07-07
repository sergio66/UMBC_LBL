      subroutine loop(outvect,ziso,mass_iso,brd,strength, 
     $   centerfreq,calcfreq,tempr,numlines,sizewave,numisotopes,lvg)

      integer numlines,sizewave,numisotopes,lvg
      real*8 ziso(numlines),mass_iso(numisotopes),brd(numlines)
      real*8 strength(numlines),outvect(sizewave)
      real*8 centerfreq(numlines),calcfreq(sizewave),tempr

      real*8 yr(sizewave),yi(sizewave),mass
      integer ii,kk

      do ii=1,sizewave
        outvect(ii)=0.0
        yr(ii)=0.0
        yi(ii)=0.0
        end do

c lorentz >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (lvg .eq. -1) then    
c        print *,'computing Lorentz lineshape for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(ziso(kk))
          call lorentz(yr,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
            end do          
          end do

c van huber >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. 0) then    
c        print *,'computing VanHuber lineshape for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(ziso(kk))
          call vhh2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
            end do          
          end do

c voigt >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. +1) then    
c        print *,'computing Voigt lineshape for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(ziso(kk))
c this is the complex voigter currently being used by Genln2
          call voigt2R(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
            end do          
          end do

        endif

      return
      end

c ***********************************************************************
c  vhh2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),sizewave) 
c     is a FAST Genln2 complex voigter, rewritten so its does VHH
c  vhh1(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),sizewave) 
c     is a slower Genln2 complex voigter, rewritten so its does VHH
c  if (brd(kk) .gt. 1e-2) then
c    vhh(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),sizewave)
c       is a quite fast real voigter, using polyomial approx. only good
c       for "wider" lines   
c DITTO
c          if (brd(kk) .gt. 1e-2) then
c this is the fast one .........
c            call voigt(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
c     $             sizewave)
c          else
c this is the FASTTTT!!!! complex voigter currently being used by Genln2
c            call voigt2R(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
c     $             sizewave)
c            end if
c 
c ***********************************************************************
      include 'lorentz.f'
      include 'vhh2RI.f'
      include 'voigt2R.f'

c      include 'vhh1.f'
c      include 'vhh.f'
c      include 'voigt.f'
