c compile with 
c   f77 -o make_mass_HITRAN2MATLAB.x -N109 -W make_mass_HITRAN2MATLAB.f

      IMPLICIT NONE

      INTEGER kNgas             !!number of gases
      INTEGER kMaxISO           !!max number of isotopes
      PARAMETER (kNgas = 40,kMaxISO = 10)

      INTEGER iWriteOutMax

      CHARACTER*80  fgasids     !!gasids for reference
      CHARACTER*80  fpartsum    !!molar mass and molecular abundance
      CHARACTER*80  fiso        !!number of isotopes per molecule
      CHARACTER*80  fout        !!dumps out massXXX.dat for load mass.dat in
                                !!runY.m suite of code

      INTEGER iUN1,iUN2,iUN3,iUN4,iERR
      INTEGER iI,iJ,iK

      INTEGER iNmols
     
      INTEGER iaNISO(kNGAS)           !!number of isotopes
      INTEGER iaaISO(kNgas,kMaxISO)   !!isotope identifiers
      REAL    raaABN(kNgas,kMaxISO)   !!relative abundance of isotopes
      REAL    raaMASS(kNgas,kMaxISO)  !!molar mass of isotope
      REAL    raaG(kNgas,kMaxISO)     !!(nuclear) degeneracy of isotopes
      REAL    raaPF296(kNgas,kMaxISO) !!total internal partition sums at 296K

      CHARACTER*80 caStr
      CHARACTER*6  caGas,caIso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c this is hard coded (right now) for HITRAN 2004
c we want to make mass04X.dat using HITRAN2004 data
      fgasids  = '/home/sergio/SPECTRA/gasids'
c copy /home/sergio/SPECTRA/NEW_HITRAN_MISC/h2k_ISO.DAT to the file below, and
c add on info for molecule 39, by looking at fpartsum
      fiso     = '/home/sergio/SPECTRA/NEW_HITRAN_MISC/h04_ISO.DAT'
      fpartsum = '/asl/data/hitran/HITRAN04/Global_Data/molparam.txt' 
      fout     = '/home/sergio/SPECTRA/mass04.dat'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c this is hard coded (right now) for HITRAN 2000
c we want to make mass00X.dat using HITRAN2000 data
      fgasids  = '/home/sergio/SPECTRA/gasids'
c      fiso     = '/asl/data/hitran/HITRAN2k/ISO.DAT'  !!!missing comma gas 28
      fiso     = '/home/sergio/SPECTRA/NEW_HITRAN_MISC/h2k_ISO.DAT'
c      fpartsum = '/asl/data/hitran/HITRAN2k/PartSum/MOLPARAM.txt' !!missing g
c damn DOS characters
c tr -d '\015\032' < h2k_MOLPARAM.txt > h2k_MOLPARAM.txt.new   OVERKILL
c tr -d '\032' < h2k_MOLPARAM.txt > h2k_MOLPARAM.txt.new
c dos2unix < /asl/data/hitran/HITRAN2k/PartSum/MOLPARAM.txt > h2k_MOLPARAM.txt
c awk '{ sub("\r$", ""); print }' dosfile.txt > unixfile.txt
      fpartsum  = '/home/sergio/SPECTRA/NEW_HITRAN_MISC/h2k_MOLPARAM.txt'
      fout     = '/home/sergio/SPECTRA/mass00.dat'
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      iWriteOutMax = 32   !! this is the number of gases H92,H96,H98 had, so
                          !! for now, stick to that number

      iUN1 = 11
      iUN2 = 12
      iUN3 = 13
      iUN4 = 14

      CALL read_fgasids(iUN1,fgasids)
      CALL read_fiso(iUN2,fiso,    iNmols,iaNISO,iaaISO)
      CALL read_fpartsum(iUN3,fpartsum,iNmols,iaNISO,iaaISO,
     $                   raaMASS,raaABN,raaG,raaPF296)
      CALL write_fout(iUN4,fout,iWriteOutMax,iNmols,iaNISO,
     $                raaMASS,raaABN,raaG,raaPF296)

      END

c************************************************************************
c this subroutine reads in fgasids
      SUBROUTINE read_fgasids(iUN1,fgasids)

      IMPLICIT NONE

      CHARACTER*80  fgasids     !!gasids for reference
      INTEGER       iUN1

      INTEGER iERR
       
      OPEN(UNIT=iUN1,FILE=fgasids, STATUS='OLD',FORM='FORMATTED',IOSTAT=IERR) 
      IF (iERR .NE. 0) THEN
        write(6,*) 'Error number ',IERR,' with fgasids file ',fgasids
        write(6,*) 'probably means file does not exist; check again!'
        STOP
        END IF      
      CLOSE(iUN1)

      RETURN
      END

c************************************************************************
c this subroutine reads in fiso
      SUBROUTINE read_fiso(iUN2,fiso,iNmols,iaNISO,iaaISO)

      IMPLICIT NONE

      INTEGER kNgas             !!number of gases
      INTEGER kMaxISO           !!max number of isotopes
      PARAMETER (kNgas = 40,kMaxISO = 10)

c input vars
      CHARACTER*80  fiso        !!number of isotopes per molecule
      INTEGER iUN2

c output vars
      INTEGER iNmols
      INTEGER iaNISO(kNGAS)          !!number of isotopes
      INTEGER iaaISO(kNgas,kMaxISO)  !!isotope identifiers

c local vars
      INTEGER iI,iJ,iK,iERR,iSum,iaISO(kMaxISO),iSpieces
      CHARACTER*80 caStr,caStr1
      CHARACTER*6  caGas,caIso

      iSum = 0
      OPEN(UNIT=iUN2,FILE=fiso,    STATUS='OLD',FORM='FORMATTED',IOSTAT=IERR) 
      IF (iERR .NE. 0) THEN
        write(6,*) 'Error number ',IERR,' with fiso file ',fiso
        write(6,*) 'probably means file does not exist; check again!'
        STOP
        END IF   

      READ(iUN2,*) caStr
      READ(iUN2,*) iNmols,iSpieces
      print *,'file says there are iNmols,iSpieces = ',iNmols,iSpieces

c now read the molecules and number of isotopes info
      READ(iUN2,*) caStr
      DO iI = 1,iNmols
        READ(iUN2,50) caStr
        DO iJ = 1,6
          caGas(iJ:iJ) = ' '
          END DO
        iK = 1
 40     CONTINUE
        IF (caStr(iK:iK) .NE. ',') THEN
          caGas(iK:iK) = caStr(iK:iK)
          iK = iK + 1
          GOTO 40
          END IF
        caIso(1:6) = caStr(iK+1:iK+5)
        READ(caIso,60) iaNISO(iI)
        iSUm = iSum + iaNISO(iI)
        print *,iI,caGas,iaNISO(iI),iSum
        END DO

      IF (iSum .NE. iSpieces) THEN
        print *,'Error iSum .NE. iSpieces'
        STOP
        END IF

c now read the molecules and isotopes identifiers info
      READ(iUN2,*) caStr
      DO iI = 1,iNmols
        DO iJ = 1,80
          caStr1(iJ:iJ) = ' '
          END DO
        READ(iUN2,50) caStr

        !!! look for the ' '
        iK = 1
 80     CONTINUE
        IF (caStr(iK:iK) .NE. ' ') THEN
          iK = iK + 1
          GOTO 80
          END IF
        !!! look for the first non' '
        iK = iK
 85     CONTINUE
        IF (caStr(iK:iK) .EQ. ' ') THEN
          iK = iK + 1
          GOTO 85
          END IF
        DO iJ = 1,80-iK+1
          caStr1(iJ:iJ) = caStr(iJ+iK-1:iJ+iK-1)
          END DO
        read(caStr1,*) (iaISO(iJ),iJ = 1,iaNISO(iI))
        DO iJ = 1,iaNISO(iI)
          iaaISO(iI,iJ) = iaISO(iJ)
          END DO
c        print *,iI,(iaISO(iJ),iJ = 1,iaNISO(iI))

       END DO

      CLOSE(iUN2)
 50   FORMAT(A80)
 60   FORMAT(I5)

      RETURN
      END

c************************************************************************
c this reads in molar mass and gas abundance
      SUBROUTINE read_fpartsum(iUN3,fpartsum,iNmols,iaNISO,iaaISO,
     $                         raaMASS,raaABN,raaG,raaPF296)

      IMPLICIT NONE

      INTEGER kNgas             !!number of gases
      INTEGER kMaxISO           !!max number of isotopes
      PARAMETER (kNgas = 40,kMaxISO = 10)

c input vars
      CHARACTER*80  fpartsum    !!molar mass and molecular abundance
      INTEGER iUN3
      INTEGER iNmols
      INTEGER iaNISO(kNGAS)          !!number of isotopes
      INTEGER iaaISO(kNgas,kMaxISO)  !!isotope identifiers

c output vars
      REAL    raaABN(kNgas,kMaxISO)   !!relative abundance of isotopes
      REAL    raaMASS(kNgas,kMaxISO)  !!molar mass of isotopes
      REAL    raaG(kNgas,kMaxISO)     !!(nuclear) degeneracy of isotopes
      REAL    raaPF296(kNgas,kMaxISO) !!total internal partition sums at 296K

c local vars
      INTEGER iI,iJ,iK,iERR,iISO
      REAL    rAbn,rPF,rG,rMM
      CHARACTER*80 caStr

      OPEN(UNIT=iUN3,FILE=fpartsum,STATUS='OLD',FORM='FORMATTED',IOSTAT=IERR) 
      IF (iERR .NE. 0) THEN
        write(6,*) 'Error number ',IERR,' with fpartsum file ',fpartsum
        write(6,*) 'probably means file does not exist; check again!'
        STOP
        END IF      

      READ(iUN3,*) caStr   !! read the header
      DO iI = 1,iNmols
        READ(iUN3,*) caStr
        DO iJ = 1,iaNISO(iI)
          read(iUN3,*) iISO,rAbn,rPF,rG,rMM
          IF (iISO .NE. iaaISO(iI,iJ)) THEN
            print *,'j = ',iJ,'isotopes do not match for gas ',iI
            STOP
          ELSE
            raaABN(iI,iJ)   = rAbn
            raaMASS(iI,iJ)  = rMM
            raaG(iI,iJ)     = rG
            raaPF296(iI,iJ) = rPF
            END IF
          END DO
        END DO 
      CLOSE(iUN3)

      RETURN
      END

c************************************************************************
      SUBROUTINE write_fout(iUN4,fout,iWriteOutMax,iNmols,iaNISO,
     $                      raaMASS,raaABN,raaG,raaPF296)

      IMPLICIT NONE

      INTEGER kNgas             !!number of gases
      INTEGER kMaxISO           !!max number of isotopes
      PARAMETER (kNgas = 40,kMaxISO = 10)

c input vars
      CHARACTER*80  fout        !!dumps out massXXX.dat for load mass.dat in
                                !!runY.m suite of code
      INTEGER iUN4,iNmols,iWriteOutMax
      INTEGER iaNISO(kNGAS)           !!number of isotopes
      REAL    raaABN(kNgas,kMaxISO)   !!relative abundance of isotopes
      REAL    raaMASS(kNgas,kMaxISO)  !!molar mass of isotopes
      REAL    raaG(kNgas,kMaxISO)     !!(nuclear) degeneracy of isotopes
      REAL    raaPF296(kNgas,kMaxISO) !!total internal partition sums at 296K

c local vars
      INTEGER iI,iJ,iK,iERR
      REAL rX,rY

      rX = -1.0
      rY =  0.0
     
      OPEN(UNIT=iUN4,FILE=fout,STATUS='NEW',FORM='FORMATTED',IOSTAT=IERR) 
      IF (iERR .NE. 0) THEN
        write(6,*) 'Error number ',IERR,' with fout file ',fout
        write(6,*) 'probably means output file already exists; check again!'
        STOP
        END IF      

      IF (iWriteOutMax .EQ. iNmols) THEN
        print *,'write out info for all the ',iNmols,' gases we found ...'
      ELSEIF (iWriteOutMax .LT. iNmols) THEN
        print *,'write out info for first ',iWriteOutMax,' gases we found ...'
      ELSEIF (iWriteOutMax .GT. iNmols) THEN
        print *,'Error cannot write out info for ',iWriteOutMax,' gases !!!'
        STOP
        END IF

      DO iI = 1,iWriteOutMax
        write(iUN4,50) iI,iaNISO(iI),rX,rY,rY
        END DO
      DO iI = 1,iWriteOutMax
        DO iJ = 1,iaNISO(iI)        
          write(iUN4,60) iI,raaMASS(iI,iJ),raaABN(iI,iJ),
     $                      INT(raaG(iI,iJ)),raaPF296(iI,iJ)
          END DO
        END DO

 50   FORMAT(I13,' ',I11 ,' ',3(E14.7,' '))
 60   FORMAT(I13,' ',F9.5,' ',E14.7,' ',I5,' ',E14.7)
      CLOSE(iUN4)

      RETURN
      END

c************************************************************************
