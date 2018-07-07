C***********************************************************************
C
C  PROGRAM        HITRD  SUBROUTINE     
C
C  PURPOSE        ROUTINE FOR READING HITRAN UNFORMATTED LINE DATA BASE
C
C  VERSION        3.0   D.P. EDWARDS   02/05/90 
C
C  DESCRIPTION    SUBROUTINE FOR READING THE LINE DATA FILE.
C                 THIS ROUTINE CAN BE USED TO READ THE UNFORMATTED 
C                 DIRECT ACCESS HITRAN FILE USING VAX SPECIFIC FORTRAN 
C                 WHEN USING A VAX OR USING STANDARD FORTRAN77 FOR 
C                 OTHER MACHINES. WHICH EVER IS USED, THE CODE SECTIONS
C                 RELEVANT TO THE OTHER OPTION SHOULD BE COMMENTED OUT 
C                 BELOW.
C
C***********************************************************************
C  QUANTITES AND STORAGE READ FROM THE HITRAN LINE DATA BASE
C***********************************************************************
C  VARIABLE   SOURCE      TYPE  DESCRIPTION
C  ---------------------------------------------------------------------
C  LSTAT      GENLN2       Int  Status of transition information.
C  IGAS       HITRAN       Int  Molecule number.
C  ISO        HITRAN       Int  Isotope number 
C                               (1=most abundant, 2=second,). 
C  WNUM       HITRAN(mod.) DPr  Line frequency [cm-1] in double pecision.
C  STREN      HITRAN(mod.) Re   Line strength  [cm-1./(kg.moles.cm-2)] 
C                               @ 296K. (HITRAN value * 6.022E26 
C                               to avoid underflows.)
C  TPROB      HITRAN       Re   Transition probability [Debyes2].
C  ABROAD     HITRAN       Re   Air-broad halfwidth  (HWHM) [cm-1/atm] 
C                                @ 296K.
C  SBROAD     HITRAN       Re   Self-broad halfwidth (HWHM) [cm-1/atm] 
C                               @ 296K.
C  ELS        HITRAN       Re   Lower-state energy [cm-1].
C  ABCOEF     HITRAN       Re   Coefficient of temperature dependance 
C                               of air-broadened halfwidth.
C  TSP        HITRAN       Re   Transition shift due to pressure 
C                               (presently empty, some coupling coeffs.
C                               inserted).
C  IUSGQ      HITRAN       Int  Upper state global quanta index.
C  ILSGQ      HITRAN       Int  Lower state global quanta index.
C  ---------------------------------------------------------------------
C             14 words total           
C             IBUF1 Total =  56  bytes for a 4 byte/word machine
C                         = 112  bytes for a 8 byte/word machine
C  ---------------------------------------------------------------------
C  USLQ       HITRAN       A*9  Upper state local quanta. 
C  BSLQ       HITRAN       A*9  Lower state local quanta. 
C  AI         HITRAN(mod.) A*3  Accuracy indices for frequency, 
C                               intensity and halfwidth 
C                               (original HITRAN written as 3I1).
C  REF        HITRAN(mod.) A*6  Indices for lookup of references for 
C                               frequency,intensity and halfwidth, 
C                               not yet used. 
C                               (original HITRAN written as 3I2).
C  IFWDPT     GENLN2       I*4  Forward pointer on data line.
C  BLANK                   A*1
C  ---------------------------------------------------------------------
C             IBUF2 Total =  32 bytes/ 8 words for a 4 byte/word machine 
C                         =  36 bytes for a 8 byte/word machine
C                            BUT the CRAY pads the characters with an
C                            extra 4 bytes so that they terminate on a
C                            word boundary.
C  ---------------------------------------------------------------------
C
C  SUBROUTINES    BLCKDT
C
C  I/O UNITS      LUNX   I/P   INPUT LINE DATA FILE   
C                 IOD    O/P   FORMATTED SUMMARY OUTPUT FILE
C
C  SYSTEM SPECIFIC                 VAX EQUIVALENCE OF BYTE DATA
C***********************************************************************
C
       SUBROUTINE HITRD(KFND)
C-----------------------------------------------------------------------
	integer mxgas, IOD, LUN, IG, I, ICC
       PARAMETER (MXGAS=15)
       INTEGER   KFND
C
C  COMMON BLOCK CONTAINING I/O UNITS
C

       COMMON/COMUNT/ IOD,LUN

C
C  COMMON BLOCKS COMHT1 AND COMHT2 CONTAIN THE VARIABLES READ IN ONE 
C  RECORD OF THE LINE DATA BASE. COMHT1 CONTAINS NUMERIC VARABLES AND 
C  COMHT2 CHARACTER.
C
       COMMON /COMHT1/ IDUM,LSTAT,IGAS,ISO,WNUM,STREN,TPROB,ABROAD,
     +                SBROAD,ELS,ABCOEF,TSP,IUSGQ,ILSGQ
       COMMON /COMHT2/ USLQ,BSLQ,AI,REF
C
       INTEGER          IDUM,LSTAT,IGAS,ISO,IUSGQ,ILSGQ,IFWDPT
       real*4             STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP
       DOUBLE PRECISION WNUM
       CHARACTER*3      AI
       CHARACTER*6      REF
       CHARACTER*9      USLQ,BSLQ 
C
C  COMMON BLOCK COMHIT IS USED TO COMMUNICATE BETWEEN HITINI AND HITRD;
C  IT IS NOT INTENDED TO BE READ OUTSIDE THESE ROUTINES.
C
       COMMON /COMHIT/ LUNX,NOGAS,IFWD,IGSNUM,IDGX,ICHX,IFL,
     +                 IREC1,IREC2,IREC
C
       INTEGER   LUNX,NOGAS,IFWD(MXGAS),
     +           IGSNUM,IDGX(MXGAS),ICHX(MXGAS),
     +           IFL,IREC1,IREC2,IREC
C
       INTEGER   IPTRS(14),II1(9),II2(5)
       LOGICAL   IFOUND
C
C$$$$$$ VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       BYTE        IBUF1(56),IBUF2(32)
C
C       EQUIVALENCE (IBUF1(1),LSTAT),(IBUF2(1),USLQ)    
C       EQUIVALENCE (IBUF2(28),IFWDPT)                    
C       EQUIVALENCE (IBUF1(21),II1),(IBUF2(1),II2)      
C$$$$$$ END OF VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
       SAVE 
C-----------------------------------------------------------------------
C
C       DATA REFM/'     5'/
C
C  GET A RECORD NO EITHER AS THE NEXT RECORD IN THE FORWARD POINTER LIST
C  OR THE NEXT RECORD IF THE FORWARD POINTER LIST IS NOT COMPLETE.
C
       IGAS = 0
   10  IF (IFL .NE. 0 .OR. IGSNUM .EQ. 0) THEN
         IREC = IREC + 1
       ELSE
         IREC = IFWD(1)
         DO 20 IG=2,IGSNUM
           IF (IREC .GT. IFWD(IG))
     +     IREC = IFWD(IG)
   20    CONTINUE
       ENDIF
C
C  IF WE HAVE DROPPED OFF THE END OF THE FILE
C  SET THE 'NO MORE LINES' FLAG AND RETURN
C
       IF (IREC .GT. IREC2) THEN
         write(6,*) 'END OF FILE'
         IGAS = -1
         RETURN
       ENDIF
C
C  START OF READ LOOP 30 
C
   30  CONTINUE
C
C$$$$$$ VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       READ(LUNX,REC=IREC)IBUF1,IBUF2 
C$$$$$$ END OF VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$ STANDARD SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       READ(LUNX,REC=IREC) LSTAT
       IF (ABS(LSTAT) .GE. 10) THEN
         READ(LUNX,REC=IREC)
     +   LSTAT,IGAS,ISO,WNUM,
     +   STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP,IUSGQ,ILSGQ,
     +   USLQ,BSLQ,AI,REF,IFWDPT

C
C  READ FORWARD POINTER BLOCK RECORD
C
       ELSE 
         READ(LUNX,REC=IREC)
     +   LSTAT,IGAS,ISO,WNUM,(IPTRS(I),I=1,14)
C     
       ENDIF
C$$$$$$ END OF STANDARD SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  PROCESS A FORWARD POINTER BLOCK IF ONE HAS BEEN READ, LSTAT=-7
C
       IF (LSTAT .EQ. -7) THEN
C
C  IGNORE RECORD IF ALL FWD PTS HAVE BEEN FOUND
C
C         WRITE(6,*) 'Forward ptr block'
         IF (IGSNUM.EQ.0 .OR. IFL.EQ.0) THEN 
           IREC = IREC + 1
           GOTO 30
C
C  UPDATE ANY FWD PTS NOT YET FOUND
C
         ELSE
C
C
C$$$$$$ VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C           DO 50 I=1,14    
C             IF (I .LE. 9) IPTRS(I) = II1(I) 
C             IF (I d.GE.10) IPTRS(I) = II2(I-9) 
C   50      CONTINUE  
C$$$$$$ END OF VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
           DO 70 IG=1,IGSNUM
             IF (IDGX(IG) .GE. IGAS .AND.
     +       IDGX(IG) .LT. IGAS+14) THEN
               IF (IFWD(IG) .EQ. 0 .AND. IFL .NE. 0) 
     +         IFL = IFL - 1
               IFWD(IG) = IPTRS(IDGX(IG)-IGAS+1) + IREC
             ENDIF
   70      CONTINUE
C
C  GO AND GET ANOTHER RECORD, WHICH MAY NOT BE NEXT IN NUMERICAL
C  SEQUENCE IF THE FWD PTS LIST HAS JUST BEEN FILLED 
C
           GOTO 10
         ENDIF
C
C  CHECK IF A LINE MIXING RECORD IS TO FOLLOW
C
cc       ELSEIF (LSTAT .GT. 10 .AND. REF .EQ. REFM) THEN
cc         IF (LMXFLG) THEN
C
cc           READ(INS,1000) LSTATM,IGASM,ISOM,WNUMM,IUSGQM,ILSGQM,USLQM,
cc     +     BSLQM,COMIX
cc 1000      FORMAT(I4,I2,I1,F12.6,2I3,2A9,3(1PE16.8))
cc         
cc          IF ((IUSGQM .NE. IUSGQ) .OR. (ILSGQM .NE. ILSGQ) .OR.
cc     +     (USLQM .NE. USLQ) .OR. (BSLQM .NE. BSLQ)) THEN
cc             WRITE(IOD,1005) IGAS,ISO,WNUM
cc 1005        FORMAT(
cc     +       1X,'!!!!!!!!!! SUBROUTINE HITRD !!!!!!!!!!',/,
cc     +       1X,'ERROR: Line mixing calculation for gas',/,
cc     +       1X,'ID =',I3,' isotope =',I2               ,/,
cc     +       1X,'Mis-match between line and mixing     ',/,
cc     +       1X,'data records on line database at      ',/,
cc     +       1X,'wavenumber = ',F12.6                  )    
cc             call mexerrmsgtxt('HITRD')
cc           ENDIF
cc         ENDIF
C
C  TEST IF THIS IS THE LAST RECORD ON THE FILE, LSTAT=-2
C  SET THE 'NO MORE LINES' FLAG AND RETURN
C
       ELSEIF (LSTAT .EQ. -2) THEN
         write(6,*) 'no more lines'
         IGAS = -1
         RETURN
       ENDIF
C
C  CHECK THAT GAS ID IS WITHIN RANGE 
C
       IF (IGAS .GT. NOGAS) THEN
         WRITE(IOD,1030)IGAS,IREC,NOGAS
 1030    FORMAT(
     +   1X,'!!!!!!!!!! SUBROUTINE HITRD !!!!!!!!!!',/,
     +   1X,'ERROR: gas ID =',I4,' in record',I5    ,/,
     +   1X,'This exceeds the stated number in line',/,
     +   1X,'file header record =',I4               )
         call mexerrmsgtxt('HITRD')
       ENDIF
C
C  SEE IF GAS IS A WANTED GAS
C
       IFOUND = .FALSE.
       DO 80 IG=1,IGSNUM
         IF (IGAS .EQ. IDGX(IG)) THEN
           IFOUND = .TRUE.
           KFND = IG
           GOTO 90
         ENDIF
   80  CONTINUE
C
   90  IF (IFOUND) THEN
C
C  A WANTED GAS, DECREMENT THE 'GASES TO BE FOUND FLAG' IF GAS NOT 
C  PREVIOUSLY ENCOUNTERED
C
         IF (IFL .NE. 0 .AND. IFWD(IG) .EQ. 0) 
     +   IFL = IFL - 1
C
C  LOG THE FORWARD POINTER FOR THE NEXT CALL 
C
         IFWD(IG) = IFWDPT + IREC
C
C  CHECK THAT THE OCCURRANCE OF THIS LINE IS THE ONE THAT IS REQUIRED:-
C  FOR ICHX=ICC, THE MOST RECENT LINES UP TO VERSION ICC WILL BE CHOSEN
C  I.E. LINES WITH LSTAT BETWEEN 10 (ORIGINAL HITRAN) AND ICC, AND ANY
C  -ICC LINES THAT HAVE BEEN SUPERCEEDED BY MORE RECENT VERSIONS. 
C
         ICC = ICHX(KFND)
         IF ((LSTAT .GE. 10 .AND. LSTAT .LE. ICC) .OR.
     +   (LSTAT .EQ. -ICC)) THEN 
C
C  THIS IS THE RECORD WE WANT SO RETURN
C
           RETURN
         ELSE
C
C  THIS IS NOT THE OCCURRANCE WE WANT SO START AGAIN
C
           GOTO 10
         ENDIF
C
C  ELSE PART FOR IFOUND =.FALSE.
C
       ELSE
C
C  IF IFL IS ZERO, GAS FOUND SHOULD BE VALID (I.E. IFOUND IS TRUE)
C  OTHERWISE GENERATE AN 'ERROR IN DATABASE' MESSAGE.
C
         IF(IFL .EQ. 0)THEN
           WRITE(IOD,1040)IREC,IGAS 
 1040      FORMAT(
     +     1X,'!!!!!!!!!! SUBROUTINE HITRD !!!!!!!!!!',/,
     +     1X,'ERROR: Line file record =',I7          ,/,
     +     1X,'Gas ID =',I6,' found and not expected.')
           call mexerrmsgtxt('HITRD')
         ENDIF
C
       ENDIF
C
C  NEXT 
       IREC = IREC + 1
       GOTO 30
C
C  END OF  LOOP 30 
C
       END




