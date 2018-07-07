      subroutine loop_chibasement2(outvect,ziso,mass_iso,brd,strength, 
     $   centerfreq,calcfreq,tempr,numlines,sizewave,numisotopes,lvg,
     $   pair,ph2o,doband)

c this uses the SERGIO chi functions

c same as loop.f except so that it is consistent with CKD defn of continnum,
c it has the basement term subtracted off
c note that whether we are using lorentz, voigt or VHH, the local term that is
c subtrcatetd off is the lorenntz term
c see /salsify/packages/Genln2/Genln2/vlocal_newshp.f

c also this function multiplies in the local chi fcn so as to make the 
c lineshape super lorentzian near the wings
c see /salsify/packages/Genln2/Genln2/localwat_newshp.f
C 
C  PURPOSE        Compute H2O Voigt lineshape modified for 
C                 super-Lorentzian near wings using chi function, 
C                 and subtract the "basement" term. 
C 
C  VERSION        1.0 Scott Hannon, July 1996 
C 
C  DESCRIPTION    This routine uses interpolations of a pre-tabulated 
C                 chi factor. The chi factor consists of two parts: a 
C                 self broadened part (which uses Tony Clough's self 
C                 continuum model), and a foreign broadened part (which 
C                 uses a spline of David Tobin's results). The so 
C                 called basement term is the value at 25 cm-1 from 
C                 line center. 


      include 'max.inc' 

      integer numlines,sizewave,numisotopes,lvg,doband
      real*8 ziso(kMaxBandLines)
      real*8 mass_iso(kMaxIsotopes)
      real*8 brd(kMaxBandLines)
      real*8 strength(kMaxBandLines),outvect(MaxLen)
      real*8 centerfreq(kMaxBandLines),calcfreq(MaxLen),tempr
      real*8 pair,ph2o

      real*8 yr(MaxLen),yi(MaxLen),mass,br,bi
      integer ip,kk,j
     
      real*8 prs,prf,chifac,base,fdif,fchi,df,chif,chis
      real*8 sergioCHIS,sergioCHIF
      
      real*8 xwc(11),S_chi(11),F_chi(11)        !!!the stored chi functions
      real*8 slope, yint
      integer iIOUN,ierr,iS,iFOUND
      character*80 caFName

C      Self and foreign pressure broadening ratios 
       PRS=PH2O/PAIR 
       PRF=(PAIR-PH2O)/PAIR 
C      Note that PRF + PRS = 1.0 

c      print *,'doband = ',doband

      iIOUN = 11
      if (doband .EQ. 1250) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1250.dat'
      elseif (doband .EQ. 1350) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1350.dat'
      elseif (doband .EQ. 1450) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1450.dat'
      elseif (doband .EQ. 1550) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1550.dat'
      elseif (doband .EQ. 1650) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1650.dat'
      elseif (doband .EQ. 1750) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1750.dat'
      elseif (doband .EQ. 1850) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1850.dat'
      elseif (doband .EQ. 1950) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself1950.dat'
      elseif (doband .EQ. 2050) then
        caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself2050.dat'
        end if
      OPEN(UNIT=iIOUN,FILE=caFName,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
        print *, 'error ',IERR,'opening file ', caFName
        STOP
        END IF
      ip = 1
      do ip = 1,11
        read(iIOUN,*) xwc(ip),S_chi(ip)
        end do
      close(iIOUN)

      if (prf .gt. 1.0D-8) then
        iIOUN = 11
        if (doband .EQ. 1250) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1250.dat'
        elseif (doband .EQ. 1350) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1350.dat'
        elseif (doband .EQ. 1450) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1450.dat'
        elseif (doband .EQ. 1550) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1550.dat'
        elseif (doband .EQ. 1650) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1650.dat'
        elseif (doband .EQ. 1750) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1750.dat'
        elseif (doband .EQ. 1850) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1850.dat'
        elseif (doband .EQ. 1950) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn1950.dat'
        elseif (doband .EQ. 2050) then
          caFName = '/home/sergio/WATER/CONTINUUM/sergiochiforn2050.dat'
          end if
        OPEN(UNIT=iIOUN,FILE=caFName,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          print *, 'error ',IERR,'opening file ', caFName
          STOP
          END IF
        ip = 1
        do ip = 1,11
          read(iIOUN,*) xwc(ip),F_chi(ip)
          end do
        close(iIOUN)
        end if

      IF (sizewave .gt. MaxLen) THEN
        print *,'in loop_chi.f, sizewave .gt. MaxLen',sizewave,MaxLen
        STOP
        END IF

      IF (numlines .gt. kMaxBandLines) THEN
        print *,'in loop_chi.f, numlines  .gt. kMaxBandLines'
        STOP
        END IF

c      1/PI = 0.31830988
c      625.0 = 25.0**2
c      basment = 0.31830988*STRPAR*WIDPAR/(625.0 + WIDPAR**2)
      
      do ip=1,sizewave
        outvect(ip)=0.0
        yr(ip)=0.0
        yi(ip)=0.0
        end do
          
      do kk=1,numlines
C       calculate the usual lineshape (lorentz,voigt or vhh) 
        mass=mass_iso(idint(ziso(kk)))
        if (lvg .eq. -1) then       !Lorentz lineshape
          call lorentz(yr,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave)
        elseif (lvg .eq. 0) then    !VanHuber lineshape
          call vhh2RI(yr,yi,calcfreq,centerfreq(kk),tempr,mass,
     $             brd(kk),sizewave)
        elseif (lvg .eq. +1) then   !Voigt lineshape
          call voigt2R(yr,yi,calcfreq,centerfreq(kk),tempr,
     $             mass,brd(kk),sizewave)
          endif
          
        CHIS = 1.0
        FDIF = 25.0 - 1.0e-5
        iS = 2
 555    CONTINUE
        IF (FDIF .LE. XWC(iS)) THEN
          iFOUND = +1
          SLOPE = (S_chi(iS) - S_chi(iS-1))/(xwc(iS) - xwc(iS-1))
          YINT = S_chi(iS) - slope*xwc(iS)
          CHIS = slope * fdif + YINT 
        ELSEIF ((iS .LT. 11) .AND. (iFound .LT. 0)) THEN
          iS = iS + 1
          GOTO 555
          END IF
        sergioCHIS = CHIS

        CHIF = 0.0
        if (prf .gt. 1.0e-8) then
          CHIF = 1.0
          FDIF = 25.0 - 1.0e-5
          iS = 2
 666      CONTINUE
          IF (FDIF .LE. XWC(iS)) THEN
            iFOUND = +1
            SLOPE = (F_chi(iS) - F_chi(iS-1))/(xwc(iS) - xwc(iS-1))
            YINT = F_chi(iS) - slope*xwc(iS)
            CHIF = slope * fdif + YINT 
          ELSEIF ((iS .LT. 11) .AND. (iFound .LT. 0)) THEN
            iS = iS + 1
            GOTO 666
            END IF
        sergioCHIF = CHIF
        end if

C       calc the basement chi at 25 cm-1 away 
C       I25WX is the half index (1001 of 2001 points, which are 0-50 cm-1 away
C       from line center)
ccc        CHIFAC=PRF*WFX(I25WX) + PRS*WSX(I25WX) 
        CHIFAC = PRF*sergioCHIF + PRS*sergioCHIS

c then compute the basement lorentz term
c this is just the wavenumber 25 cm-1 away from line center
        call vhh2RI(br,bi,centerfreq(kk)-25,centerfreq(kk),
     $             tempr,mass,brd(kk),1)

c combine basement lorentz term with basement chi term, at 25cm-1
        base=chifac*br

C       Apply the chi factor to the Voigt result and subtract the 
C       basement term. 
        DO IP=1,sizewave
          CHIS = 1.0
          CHIF = 1.0
          FDIF = DABS(calcfreq(IP) - centerfreq(kk)) 
          iFOUND = -1
          iS = 2
          SLOPE = 0.0
          YINT = 1.0
          IF (FDIF .LE. 25.0) THEN
            iS = 2
 777        CONTINUE
            IF (FDIF .LE. XWC(iS)) THEN
              iFOUND = +1
              SLOPE = (S_chi(iS) - S_chi(iS-1))/(xwc(iS) - xwc(iS-1))
              YINT = S_chi(iS) - slope*xwc(iS)
              CHIS = slope * fdif + YINT 
           ELSEIF ((iS .LT. 11) .AND. (iFound .LT. 0)) THEN
              iS = iS + 1
              GOTO 777
              END IF
            END IF

          IF ((FDIF .LE. 25.0) .AND. (prf .gt. 1.0d-8)) THEN
            iFOUND = -1
            iS = 2
            SLOPE = 0.0
            YINT = 1.0
 888        CONTINUE
            IF (FDIF .LE. XWC(iS)) THEN
              iFOUND = +1
              SLOPE = (F_chi(iS) - F_chi(iS-1))/(xwc(iS) - xwc(iS-1))
              YINT = F_chi(iS) - slope*xwc(iS)
              CHIF = slope * fdif + YINT 
          ELSEIF ((iS .LT. 11) .AND. (iFound .LT. 0)) THEN
              iS = iS + 1
              GOTO 888
              END IF
            END IF

C          CHIF = WFX(J) + (WFX(J+1) - WFX(J) )*DF 
C          CHIS = WSX(J) + (WSX(J+1) - WSX(J) )*DF 
          CHIFAC=PRF*CHIF + PRS*CHIS 
c
c         NOTE THIS EXPLICITLY SAYS THAT K(LOCAL) = 0  IF |V-VO| > 25
          if (abs(calcfreq(ip)-centerfreq(kk)) .le. 25.0) then
            outvect(IP)=outvect(ip)+(CHIFAC*YR(IP)-BASE)*strength(kk)
            end if

          ENDDO 

        end do

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
