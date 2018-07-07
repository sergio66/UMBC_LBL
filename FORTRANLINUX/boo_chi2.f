
      implicit none 

      real*8 xwc(11),S_chi(11),F_chi(11),chis      !!!the stored chi functions
      real*8 slope, yint, fdif,df
      integer iIOUN,ierr,iS,ip,ii,iFound
      character*80 caFName

      iIOUN = 11
      caFName = '/home/sergio/WATER/CONTINUUM/sergiochiself.dat'
      OPEN(UNIT=iIOUN,FILE=caFName,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
        print *, 'error ',IERR,'opening file ', caFName
        STOP
        END IF
      ip = 1
      do ip = 1,11
        read(iIOUN,*) xwc(ip),S_chi(ip)
        print *, xwc(ip),S_chi(ip)
        end do
      close(iIOUN)

      caFName = '/home/sergio/WATER/CONTINUUM/chichi.dat'
      OPEN(UNIT=iIOUN,FILE=caFName,STATUS='UNKNOWN',FORM='FORMATTED',
     $    IOSTAT=IERR)
      DO iI = 1, 20000
        CHIS = 1.0
        FDIF = 0.0025 * (iI-1)
        iFOUND = -1
        iS = 2
        SLOPE = 0.0
        YINT = 1.0
        IF (FDIF .LE. 25.0) THEN
          iS = 2
 777      CONTINUE
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
        write(iIOUN,*) fdif,chis,iS,slope,yint
        END DO
      close(iIOUN)

      stop
      end
