      subroutine loop(outvect,ziso,mass_iso,brd,strength, 
     $   centerfreq,calcfreq,tempr,numlines,sizewave,numisotopes,lvg)

c      = -3   voigt                  voigt1RI        HUMLICEK = used in HITTRAN LM 2018 
c rlvg = -2 : vhh without pedestal   vhh2RI_basement VOIVEC = GENLN2
c        -1 : lerentz
c         0 : vanhuber               vhh2RI          VOIVEC = GENLN2
c         1 : voigt                  voigt2RI        VOIVEC = GENLN2
c         2 : speed dependent
c      = +3   vanhuber               vhh1RI           HUMLICEK = used in HITTRAN LM 2018

c      include 'max.inc' 
      include '/home/sergio/SPECTRA/FORTRANFILES/max.inc'

      integer numlines,sizewave,numisotopes,lvg
      real*8 ziso(kMaxBandLines)
      real*8 mass_iso(kMaxIsotopes)
      real*8 brd(kMaxBandLines)
      real*8 strength(kMaxBandLines),outvect(MaxLen)
      real*8 centerfreq(kMaxBandLines),calcfreq(MaxLen),tempr

      real*8 yr(MaxLen),yi(MaxLen),mass
      integer ii,kk

c      IF (sizewave .gt. MaxLen) THEN
c        print *,'in loop.f, n_in .gt. MaxLen',n_in,MaxLen
c        STOP
c      END IF

      IF (numlines .gt. kMaxBandLines) THEN
        print *,'in loop.f, numlines   > kMaxBandLines',numlines,kMaxBandLines
        STOP
      END IF

      do ii=1,sizewave
        outvect(ii)=0.0
        yr(ii)=0.0
        yi(ii)=0.0
      end do

c voigt humlicek >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (lvg .eq. -3) then    
c        print *,'computing voigt humlicek for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
c this is the complex voigter from Humlicek
          call voigt1RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
          end do          
        end do

c GENLN2 voigt without basement >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. -2) then    
c        print *,'computing GENLN2 voivec VHH without pedestal for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
          call vhh2RI_basement(yr,yi,calcfreq,centerfreq(kk),tempr,
     $             mass,brd(kk),sizewave)
c          call vhh2RI_pedestal(yr,yi,calcfreq,centerfreq(kk),tempr,
c     $             mass,brd(kk),sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
          end do          
        end do

c lorentz >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. -1) then    
c        print *,'computing Lorentz lineshape for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
          call lorentz(yr,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
          end do          
        end do

c GENLN2 voivec van huber >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. 0) then
c        print *,'computing GENLN2 voivec VHH lineshape for ',numlines,' lines'          
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
c          print *,kk,ziso(kk),idint(ziso(kk)),numisotopes,mass
          call vhh2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
     $             sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
          end do          
        end do

c GENLN2 voivec voigt >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. +1) then    
c        print *,'computing Genln2 voivec Voigt lineshape for ',numlines,' lines'    
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
c this is the complex voigter currently being used by Genln2 (voivec)
          call voigt2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
          end do          
        end do

c speed dependent >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. +2) then    
c        print *,'computing SpeedDependent Voigt lineshape for ',numlines,' lines'    
        print *,'ERROR cannot do SpeedDependent Voigt lineshape for ',numlines,' lines'            
        stop
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
c this is the complex voigter currently being used by Genln2 (voivec)
c          call SDvhh2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
c     $             sizewave)
          do ii=1,sizewave
            outvect(ii)=outvect(ii)+yr(ii)*strength(kk)
          end do          
        end do

c van huber Humlicek >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (lvg .eq. +3) then
c        print *,'computing van huber humlicek for ',numlines,' lines'          
c this uses Humlicek        
        do kk=1,numlines
          !compute the line shape
          mass=mass_iso(idint(ziso(kk)))
c          print *,kk,ziso(kk),idint(ziso(kk)),numisotopes,mass
          call vhh1RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
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
c            call voigt2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,brd(kk),
c     $             sizewave)
c          end if
c 
c ***********************************************************************
      include 'lorentz.F'

c RI = real+ imag, R = real only
      include 'vhh2RI.F'
      include 'voigt2RI.F'

c RI = real+ imag, R = real only
      include 'vhh1RI.F'
      include 'voigt1RI.F'

c      include 'SpeedDependentVoigt/SDvhh2RI.F'
c      include 'SpeedDependentVoigt/cpf3.F'
c      include 'SpeedDependentVoigt/cpf.F'
c      include 'SpeedDependentVoigt/pCqSDHC.F'
c      include 'SpeedDependentVoigt/q_loop_SDV.F'
c      include 'SpeedDependentVoigt/loop_pCqSDHC.F'
c      include 'SpeedDependentVoigt/qsdhc.F'
c      include 'SpeedDependentVoigt/qSDV.F'

ccc      include 'vhh1.f'
ccc      include 'vhh.f'
ccc      include 'voigt.f'
