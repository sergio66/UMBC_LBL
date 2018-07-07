      subroutine loop_basement(outvect,ziso,mass_iso,brd,strength, 
     $   centerfreq,calcfreq,tempr,numlines,sizewave,numisotopes,lvg)

      include 'max.inc'

c same as loop.f except so that it is consistent with CKD defn of continnum,
c it has the basement term subtracted off
c see /salsify/packages/Genln2/Genln2/vlocal_newshp.f
      integer numlines,sizewave,numisotopes,lvg
      real*8 ziso(kMaxBandLines)
      real*8 mass_iso(kMaxIsotopes)
      real*8 brd(kMaxBandLines)
      real*8 strength(kMaxBandLines),outvect(MaxLen)
      real*8 centerfreq(kMaxBandLines),calcfreq(MaxLen),tempr

      real*8 yr(MaxLen),yi(MaxLen),mass,br,bi
      integer ii,kk
      integer explicit_turn_off

c      1/PI = 0.31830988
c      625.0 = 25.0**2
c       basment = 0.31830988*STRPAR*WIDPAR/(625.0 + WIDPAR**2)
      
      IF (numlines .gt. kMaxBandLines) THEN
        print *,'in loop_loc.f, numlines  .gt. kMaxBandLines'
        STOP
        END IF

      do ii=1,sizewave
        outvect(ii)=0.0
        yr(ii)=0.0
        yi(ii)=0.0
        end do

      explicit_turn_off=1   !this makes sure k(local) = 0 for |v-v0| > 25

      if (explicit_turn_off .gt. 0) then
        do kk=1,numlines
          mass=mass_iso(idint(ziso(kk)))
c compute the basement term
c this is just the wavenumber 25 cm-1 away from line center
          call vhh2RI(br,bi,centerfreq(kk)-25,centerfreq(kk),
     $             tempr,mass,brd(kk),1)
          
          if (lvg .eq. -1) then    !Lorentz lineshape
            call lorentz_base2(yr,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave,br)        
          elseif (lvg .eq. 0) then    !VanHuber lineshape
            call vhh2RI_base2(yr,yi,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave,br)
          elseif (lvg .eq. +1) then    !Voigt lineshape
            call voigt2R_base2(yr,yi,calcfreq,centerfreq(kk),tempr,
     $             mass,brd(kk),sizewave,br)
            end if
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
            end do          
          end do

      elseif (explicit_turn_off .lt. 0) then
        do kk=1,numlines
          mass=mass_iso(idint(ziso(kk)))
          if (lvg .eq. -1) then    !Lorentz lineshape
            call lorentz_basement(yr,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave)
          elseif (lvg .eq. 0) then    !VanHuber lineshape
            call vhh2RI_basement(yr,yi,calcfreq,centerfreq(kk),tempr,
     $             mass,brd(kk),sizewave)
          elseif (lvg .eq. +1) then    !Voigt lineshape
            call voigt2R_basement(yr,yi,calcfreq,centerfreq(kk),tempr,
     $             mass,brd(kk),sizewave)
            end if
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
            end do          
          end do
        end if

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
      include 'lorentz.F'
      include 'vhh2RI.F'
      include 'voigt2R.F'
