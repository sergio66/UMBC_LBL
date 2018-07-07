#include <fintrf.h>	
#include "mlxin.h"
	SUBROUTINE mexFunction(NLHS, PLHS, NRHS, PRHS)
	POINTER PLHS(*), PRHS(*)
	INTEGER NLHS NRHS
	INTEGER M, N, SIZE, XP, YP
	INTEGER MXGETM, MXGETN 
	POINTER MXGETPR, MXCREATEFULL

	pointer LINCT,WNUMHX, ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB
	POINTER ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,ZILSGQ
	POINTER ZUSLQ,ZBSLQ,ZAI,ZREF

	POINTER WNUMLO,WNUMHI,idg,ich,str
	CHARACTER*80 INFILE
	PLHS(1) = MXCREATEFULL(1,1,0)
	PLHS(2) = MXCREATEFULL(1,1,0)
	PLHS(3) = MXCREATEFULL(1,MXLIN,0)
	PLHS(4) = MXCREATEFULL(1,1,0)
	PLHS(5) = MXCREATEFULL(1,MXLIN,0)
	PLHS(6) = MXCREATEFULL(1,MXLIN,0)
	PLHS(7) = MXCREATEFULL(1,MXLIN,0)
	PLHS(8) = MXCREATEFULL(1,MXLIN,0)
	PLHS(9) = MXCREATEFULL(1,MXLIN,0)
	PLHS(10) = MXCREATEFULL(1,MXLIN,0)
	PLHS(11) = MXCREATEFULL(1,MXLIN,0)
	PLHS(12) = MXCREATEFULL(1,MXLIN,0)
	PLHS(13) = MXCREATEFULL(1,MXLIN,0)
	PLHS(14) = MXCREATEFULL(1,MXLIN,0)
	PLHS(15) = MXCREATEFULL(1,MXLIN,0)
	PLHS(16) = MXCREATEFULL(MXLIN,9,0)
	PLHS(17) = MXCREATEFULL(MXLIN,9,0)
	PLHS(18) = MXCREATEFULL(MXLIN,3,0)
	PLHS(19) = MXCREATEFULL(MXLIN,6,0)

	LINCT = MXGETPR(PLHS(1))
	WNUMHX = MXGETPR(PLHS(2))
	ZLSTAT = MXGETPR(PLHS(3))
	ZIGAS = MXGETPR(PLHS(4))
	ZISO = MXGETPR(PLHS(5))
	ZWNUM = MXGETPR(PLHS(6))
	ZSTREN = MXGETPR(PLHS(7))
	ZTPROB = MXGETPR(PLHS(8))
	ZABROAD = MXGETPR(PLHS(9))
	ZSBROAD = MXGETPR(PLHS(10))
	ZELS = MXGETPR(PLHS(11))
	ZABCOEF = MXGETPR(PLHS(12))
	ZTSP = MXGETPR(PLHS(13))
	ZIUSGQ = MXGETPR(PLHS(14))
	ZILSGQ = MXGETPR(PLHS(15))
	ZUSLQ = MXGETPR(PLHS(16))
	ZBSLQ = MXGETPR(PLHS(17))
	ZAI = MXGETPR(PLHS(18))
	ZREF = MXGETPR(PLHS(19))

	WNUMLO = MXGETPR(PRHS(1))
	WNUMHI = MXGETPR(PRHS(2))
        m= mxgetm(prhs(3))
        n= mxgetn(prhs(3))
        size=max(m,n)
	print*,size
        status=mxgetstring(prhs(3),INFILE,size)
c	INFILE = MXGETPR(PRHS(3))
	idg = MXGETPR(PRHS(4))
	ich = MXGETPR(PRHS(5))
	str = MXGETPR(PRHS(6))

      call hittomat(%val(WNUMLO), %val(WNUMHI),INFILE,%val(idg),
     +       %val(ich),%val(str),
c outputs next
     +       %val(LINCT),%val(WNUMHX),%val(ZLSTAT),%val(ZIGAS),
     +       %val(ZISO),%val(ZWNUM),%val(ZSTREN),%val(ZTPROB),
     +       %val(ZABROAD),%val(ZSBROAD),%val(ZELS),%val(ZABCOEF),
     +       %val(ZTSP),%val(ZIUSGQ),
     +       %val(ZILSGQ),%val(ZUSLQ),%val(ZBSLQ),%val(ZAI),
     +       %val(ZREF))
	END
      subroutine hittomat(WNUMLO, WNUMHI,INFILE,idgo,icho,str,
     +       LINCT2,WNUMHX,ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB,
     +       ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,
     +       ZILSGQ,ZUSLQ,ZBSLQ,ZAI,ZREF)

C***********************************************************************
C
C  Program       HITTOMAT
C
C  Purpose       Program copys lines from HITRAN unformatted data base
C                and puts them in MATLAB .MAT format
C
C***********************************************************************
#include "mlxin.h"

      PARAMETER (MXGAS=32)

C
C  COMMON BLOCK CONTAINING I/O UNITS
C
      COMMON /COMUNT/ IOD,LUN

      INTEGER IOD
C
C  COMMON BLOCKS COMHT1 AND COMHT2 CONTAIN THE VARIABLES READ IN ONE
C  RECORD OF THE LINE DATA BASE. HITC1 CONTAINS NUMERIC VARIABLES AND
C  HITC2 CHARACTER.
C
       COMMON /COMHT1/ IDUM,LSTAT,IGAS,ISO,WNUM,STREN,TPROB,ABROAD,
     +                SBROAD,ELS,ABCOEF,TSP,IUSGQ,ILSGQ
       COMMON /COMHT2/ USLQ,BSLQ,AI,REF
C
       INTEGER          IDUM,LSTAT,IGAS,ISO,IUSGQ,ILSGQ
       REAL*8             STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP
       DOUBLE PRECISION WNUM
       CHARACTER*3      AI
       CHARACTER*6      REF
       CHARACTER*9      USLQ,BSLQ
C
       REAL*8       WNUMLO,WNUMHI,GSTR(MXGAS),STR, linct2
       INTEGER    KFND,ISUM,I,IDG(MXGAS),ICH(MXGAS),LINCT
	real*8   idgo, icho
       CHARACTER  TITLB*48,INFILE*80,OUTFILE*80
C
C  DECLARE ARRAYS TO HOLD LINE DATA FOR MATLAB
C
       REAL*8 WNUMHX,ZLSTAT(MXLIN),ZIGAS,ZISO(MXLIN),
     +                  ZWNUM(MXLIN),ZSTREN(MXLIN),ZTPROB(MXLIN),
     +                  ZABROAD(MXLIN),ZSBROAD(MXLIN),ZELS(MXLIN),
     +                  ZABCOEF(MXLIN),ZTSP(MXLIN),ZIUSGQ(MXLIN),
     +                  ZILSGQ(MXLIN),ZUSLQ(MXLIN,9),ZBSLQ(MXLIN,9),
     +                  ZAI(MXLIN,3),ZREF(MXLIN,6)
C
C  Set i/o unit numbers
C

	idg(1)=idgo
	ich(1)=icho
        LUN = 12
C
C  Read minimum frequency
C

c       WRITE(6,*) 'Enter minimum frequency:'
c       READ(5,*) WNUMLO

C
C  Read maximum frequency
C

c       WRITE(6,*) 'Enter maximum frequency:'
c       READ(5,*) WNUMHI

C
C  Read number of gases
C
C     WRITE(6,*) 'Enter number of gases to select:'
C     READ(5,*) ISUM
       ISUM=1
C
C  Read in gas id numbers
c       DO 10, I=1,ISUM
c         WRITE(6,*) 'Enter id of gas:'
c         READ(5,*) IDG(I)
c         WRITE(6,*) 'Enter line version #:'
c         READ(5,*) ICH(I)
c         WRITE(6,*) 'Enter minimum strength:'
c         READ(5,*) STR
         GSTR(I) = STR * 6.022E26
c10     CONTINUE
C
C  Read in hitran linebase filename
C

C       WRITE(6,*) 'Enter name of hitran linebase file:'
C       READ(5,120) INFILE
C120    FORMAT (A)

C
C  Read in hitran linebase filename
C

C       WRITE(6,*) 
C     +  'Enter name of MATLAB line data output file (xxx.mat) :'
C       READ(5,120) OUTFILE

C
C  Initialize line counter and maximum strength
C
       LINCT=0
       MXSTR = 0.0
C
C  Change wnum high limit to double precision
C
       WNUMHX = WNUMHI
C
C  Initialize hitran file
C
       CALL INIT(INFILE,WNUMLO,LUN,ISUM,IDG,ICH,TITLB)
C
C  Read entry from hitran file
C
200    CALL HITRD(KFND)
      IF (IGAS .GT. 0 .AND. WNUM .LE. WNUMHX) THEN
         IF (STREN .GT. GSTR(KFND)) THEN
C
C  Save valid entry data for gas id
C
           LINCT=LINCT+1
           ZLSTAT(LINCT)=LSTAT
           ZISO(LINCT)=ISO
           ZWNUM(LINCT)=WNUM
           ZSTREN(LINCT) = STREN/6.022E26
           ZTPROB(LINCT)=TPROB
           ZABROAD(LINCT)=ABROAD
           ZSBROAD(LINCT)=SBROAD
           ZELS(LINCT)=ELS
           ZABCOEF(LINCT)=ABCOEF
           ZTSP(LINCT)=TSP
           ZIUSGQ(LINCT)=IUSGQ
           ZILSGQ(LINCT)=ILSGQ
           DO 300 I=1,9
              ZUSLQ(LINCT,I)=ICHAR(USLQ(I:I))
              ZBSLQ(LINCT,I)=ICHAR(BSLQ(I:I))
300        CONTINUE
           DO 310 I=1,3
              ZAI(LINCT,I)=ICHAR(AI(I:I))
310        CONTINUE
           DO 320 I=1,6
              ZREF(LINCT,I)=ICHAR(REF(I:I))
320        CONTINUE
         END IF
C
C  Read another entry from file
C
         GOTO 200
       END IF
C
C      No need to keep all gas id's since this program reads in only one
C
       ZIGAS=IDG(1)
C
C  Print # of lines selected
C
       WRITE(6,399) LINCT
 399   FORMAT('Number of lines selected: ',I6)
C
C  Write out selected data
C
C

c        CALL wrout(LINCT,WNUMHX,ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB,
c     +       ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,
c     +       ZILSGQ,ZUSLQ,ZBSLQ,ZAI,ZREF)
C
	linct2=linct
C       STOP
       END

C***********************************************************************
C
C  PROGRAM        INIT  SUBROUTINE
C
C  PURPOSE        INITIALISING ROUTINE FOR READING LINE DATA FILES
C
C  VERSION        3.0   D.P. EDWARDS, D. ROISIN   02/05/90 
C
C  DESCRIPTION    SUBROUTINE FOR INITIALIZING THE LINE DATA FILE.
C                 THIS ROUTINE CAN BE USED TO INITIALIZE THE 
C                 UNFORMATTED DIRECT ACCESS HITRAN FILE USING VAX 
C                 SPECIFIC FORTRAN WHEN USING A VAX OR USING STANDARD 
C                 FORTRAN77 FOR OTHER MACHINES. WHICH EVER IS USED, 
C                 THE CODE SECTIONS RELEVANT TO THE OTHER OPTION 
C                 SHOULD BE COMMENTED OUT BELOW.
C
C  ARGUMENTS      INFILE  C*80 I/P HITRAN LINEBASE FILENAME
C                 WNUMLO  I*4  I/P LOWEST WAVENUMBER OF INTEREST [cm-1]
C                 LUN     I*4  I/P STREAM FOR READING LINEBASE
C                 ISUM    I*4  I/P NUMBER OF GASES TO BE READ
C                 IDG     I*4  I/P ID'S OF REQUIRED GASES
C                 ICH     I*4  I/P LINE VERSIONS OF REQUIRED GASES
C                 TITLB   C*48 O/P TITLE OF LINEBASE  
C
C  I/O UNITS      LUN    I/P   INPUT LINE DATA FILE   
C                 IOD    O/P   FORMATTED SUMMARY OUTPUT FILE
C
C  SYSTEM SPECIFIC                 VAX OPEN STATEMENT: READONLY, RECL
C                                  CRAY OPEN STATEMENT: RECL
C                                  VAX BYTE DATA EQUIVALENCE STATEMENTS
C***********************************************************************
C
       SUBROUTINE INIT
     + (INFILE,WNUMLO,LUN,ISUM,IDG,ICH,TITLB)
C-----------------------------------------------------------------------
       PARAMETER (MXGAS=15)
C
C  COMMON BLOCK CONTAINING I/O UNITS
C
       COMMON /COMUNT/ IOD

       INTEGER IOD
C
C  COMMON BLOCKS HITC1 AND HITC2 CONTAIN THE VARIABLES READ IN ONE 
C  RECORD OF THE LINE DATA BASE. HITC1 CONTAINS NUMERIC VARABLES AND 
C  HITC2 CHARACTER.
C
       COMMON /COMHT1/ IDUM,LSTAT,IGAS,ISO,WNUM,STREN,TPROB,ABROAD,
     +                SBROAD,ELS,ABCOEF,TSP,IUSGQ,ILSGQ
       COMMON /COMHT2/ USLQ,BSLQ,AI,REF
C
       INTEGER          IDUM,LSTAT,IGAS,ISO,IUSGQ,ILSGQ,IFWDPT
       REAL             STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP
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
       INTEGER     LUNX,NOGAS,IFWD(MXGAS),
     +             IGSNUM,IDGX(MXGAS),ICHX(MXGAS),
     +             IFL,IREC1,IREC2,IREC
C
       INTEGER          IDG(MXGAS),ICH(MXGAS)
       DOUBLE PRECISION WNUMLX
       LOGICAL          INITIAL
       CHARACTER        TITLB*48,INFILE*80,FILMIX*80,WORDIN*70
       INITIAL = .TRUE.
C
C$$$$$$ VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       
C       BYTE        IBUF1(56),IBUF2(32)                  
C       EQUIVALENCE (IBUF1(1),LSTAT),(IBUF2(1),USLQ)    
C       EQUIVALENCE (IBUF2(28),IFWDPT)                  
C      
C$$$$$$ END OF VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C       SAVE          LLS COMMENTED OUT
C-----------------------------------------------------------------------
C       
C  SET UP COMMON BLOCK VARAIABLES
C
       LUNX = LUN
       IGSNUM = ISUM
       DO 10 IG=1,IGSNUM
         IDGX(IG) = IDG(IG)
         ICHX(IG) = ICH(IG)
   10  CONTINUE
C
C  CHECK INITIALISATION STATUS AND OPEN LINEBASE FILE
C
       IF (INITIAL) THEN
C
C
C$$$$$$ STANDARD SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C  FOR USE WITH HITLIN MACH=1 OPTION
C  4 BYTE/WORD MACHINES, RECL IN WORDS
C
         OPEN(LUNX,FILE=INFILE,STATUS='OLD',
     +   IOSTAT=IOST,ACCESS='DIRECT',RECL=22)
C$$$$$$ END OF STANDARD SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$ VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C  FOR USE WITH HITLIN MACH=2 OPTION
C  4 BYTE/WORD VAX, RECL IN WORDS
C
C         OPEN(LUNX,FILE=INFILE,STATUS='OLD',
C     +   IOSTAT=IOST,ACCESS='DIRECT',RECL=22,READONLY)             
C$$$$$$ END OF VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$ CRAY SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C  FOR USE WITH HITLIN MACH=3 OPTION
C  8 BYTE/WORD CRAY, RECL IN BYTES
C
C         OPEN(LUNX,FILE=INFILE,STATUS='OLD',
C     +   IOSTAT=IOST,ACCESS='DIRECT',RECL=152)
C$$$$$$ END OF CRAY SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  CHECK STATUS OF .BIN FILE AFTER OPENING
C
         IF (IOST.NE.0) THEN
           WRITE(IOD,1000)IOST
 1000      FORMAT(
     +     1X,'!!!!!!!!!! SUBROUTINE HITINI !!!!!!!!!!',/,
     +     1X,'ERROR: IOST =',I5                       ,/,
     +     1X,'TRYING TO OPEN LINE DATA FILE          ')
           call mexerrmsgtxt('HITINI')
         ENDIF
C
C  GET FROM FILE HEADER RECORD 
C
         READ(LUNX,REC=1) LSTAT,NOGAS,
     +   IREC1,IREC2,TITLB
       ENDIF
C
C  DO SEARCH IN DOUBLE PRECISION SINCE WNUM IS IN DP
C       
       WNUMLX = WNUMLO
C
C  BINARY SEARCH FOR 1ST RECORD: FOLLOWS CODE IN THE IBM SLUP ROUTINE:
C  K IS FIRST RECORD, L IS LAST
C
       K = IREC1
       L = IREC2
   30  IREC = (K + L)/2
C
C  ONLY NEED WNUM FOR SEARCH
C
C$$$$$$ VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       READ(LUNX,REC=IREC) (IBUF1(I),I=1,20) 
C$$$$$$ END OF VAX SPECIFIC SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$ STANDARD SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       READ(LUNX,REC=IREC) LSTAT,IGAS,ISO,WNUM
C$$$$$$ END OF STANDARD SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
       IF (WNUM .LT. WNUMLX) THEN
C
C  REQUIRED WN HIGHER, SO TRY HALF WAY UP TO CURRENT L
C
         K = IREC
       ELSE
C
C  REQUIRED WN LOWER OR EQUAL, SO TRY HALF WAY DOWN TO CURRENT K
C
         L = IREC
       ENDIF
C
C  IF K & L DIFFER BY 1, HAVE FOUND REQUIRED LOCATION WHERE THE K 
C  WAVENUMBER IS LOWER THAN WNUMLO, AND THE L WAVENUMBER IS GREATER 
C  THAN OR EQUAL TO WNUMLO. CHOOSE RECORD L (NOTE K<L). IF THERE IS 
C  MORE THAN 1 RECORD AT EXACTLY WAVENUMBER WNUMLO, WILL FINISH 
C  POINTING TO FIRST.  FORWARD POINTER BLOCK ARE LABELLED WITH 
C  WAVENUMBER OF NEXT LINE TO EXPLOIT THIS.
C
       IF (K+1 .NE. L) GO TO 30
C
C  NOTE HITRD INCREMENTS IREC BEFORE READING SO WILL START ON L
C
       IREC = K
       IF (IREC .EQ. IREC1) IREC=IREC-1
C
C  SET UP ARRAY OF FORWARD POINTERS TO BE EMPTY:
C
       DO 40 IG=1,IGSNUM
         IFWD(IG) = 0
   40  CONTINUE
C
       IFL = IGSNUM

       RETURN
       END



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
       REAL             STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP
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
C             IF (I .GE.10) IPTRS(I) = II2(I-9) 
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






