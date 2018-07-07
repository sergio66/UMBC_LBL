      subroutine mexFunction(nlhs,plhs,nrhs,prhs)

c      gateway function for matlab mex file.
c      use hittomat.m to run this file.

      include 'fintrf.h'

c      pointers to arguments
      POINTER PLHS(*), PRHS(*)
c      number of arguments
      INTEGER NLHS, NRHS
      INTEGER M, N, SIZE, XP, YP,linct2
      INTEGER MXGETM, MXGETN, status, mxgetstring
      POINTER MXGETPR, MXCREATEFULL

      pointer LINCT,WNUMHX, ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB
      POINTER ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,ZILSGQ
      POINTER ZUSLQ,ZBSLQ,ZAI,ZREF, GASID

      POINTER WNUMLO,WNUMHI,idg,ich,str, isum
      CHARACTER*80 INFILE

c      get pointers to the input arguments
      WNUMLO = MXGETPR(PRHS(1))
      WNUMHI = MXGETPR(PRHS(2))
        m= mxgetm(prhs(3))
        n= mxgetn(prhs(3))
        size=max(m,n)
c      can't use %val for strings.
        status=mxgetstring(prhs(3),INFILE,size)
      idg = MXGETPR(PRHS(4))
      ich = MXGETPR(PRHS(5))
      str = MXGETPR(PRHS(6))
      isum= MXGETPR(PRHS(7))
c      run hittomat once to see how much space we need.
      call hittomatd(%val(WNUMLO), %val(WNUMHI),INFILE,%val(idg),
     +       %val(ich),%val(str),linct2, %val(isum))

c       Note: be default fortran sends a ponter to your data when calling
c       a routine. By using %val Fortran sends a copy of the data itself
c       and not a pointer to your routine. When the data comes into the 
c      the routine, Fortran assumes that it is just a pointer. In this way
c      you can use pointers in Fortran.

      print*,'linct2',linct2
c      use matlab to allocate space for output arguments

      PLHS(1) = MXCREATEFULL(1,1,0)
      PLHS(2) = MXCREATEFULL(1,1,0)
      PLHS(3) = MXCREATEFULL(1,linct2,0)
      PLHS(4) = MXCREATEFULL(1,1,0)
      PLHS(5) = MXCREATEFULL(1,LINCT2,0)
      PLHS(6) = MXCREATEFULL(1,LINCT2,0)
      PLHS(7) = MXCREATEFULL(1,LINCT2,0)
      PLHS(8) = MXCREATEFULL(1,LINCT2,0)
      PLHS(9) = MXCREATEFULL(1,LINCT2,0)
      PLHS(10) = MXCREATEFULL(1,LINCT2,0)
      PLHS(11) = MXCREATEFULL(1,LINCT2,0)
      PLHS(12) = MXCREATEFULL(1,LINCT2,0)
      PLHS(13) = MXCREATEFULL(1,LINCT2,0)
      PLHS(14) = MXCREATEFULL(1,LINCT2,0)
      PLHS(15) = MXCREATEFULL(1,LINCT2,0)

      PLHS(16) = MXCREATEFULL(9,LINCT2,0)
      PLHS(17) = MXCREATEFULL(9,LINCT2,0)
      PLHS(18) = MXCREATEFULL(3,LINCT2,0)
      PLHS(19) = MXCREATEFULL(6,LINCT2,0)
      PLHS(20) = MXCREATEFULL(1,LINCT2,0)

c      PLHS(16) = MXCREATEFULL(1,9*LINCT2,0)
c      PLHS(17) = MXCREATEFULL(1,9*LINCT2,0)
c      PLHS(18) = MXCREATEFULL(1,3*LINCT2,0)
c      PLHS(19) = MXCREATEFULL(1,6*LINCT2,0)

c      find the pointers for the arrays we just created.

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
      GASID = MXGETPR(PLHS(20))
      call hittomat(%val(WNUMLO), %val(WNUMHI),INFILE,%val(idg),
     +       %val(ich),%val(str),
c outputs next
     +       %val(LINCT),%val(WNUMHX),%val(ZLSTAT),%val(ZIGAS),
     +       %val(ZISO),%val(ZWNUM),%val(ZSTREN),%val(ZTPROB),
     +       %val(ZABROAD),%val(ZSBROAD),%val(ZELS),%val(ZABCOEF),
     +       %val(ZTSP),%val(ZIUSGQ),
     +       %val(ZILSGQ),%val(ZUSLQ),%val(ZBSLQ),%val(ZAI),
     +       %val(ZREF),linct2,%val(isum),%val(GASID))
      END
      subroutine hittomat(WNUMLO, WNUMHI,INFILE,idgo,icho,str,
     +       LINCT2,WNUMHX,ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB,
     +       ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,
     +       ZILSGQ,ZUSLQ,ZBSLQ,ZAI,ZREF,linct3, isum2,GASID)

C***********************************************************************
C
C  Program       HITTOMAT
C
C  Purpose       Program copys lines from HITRAN unformatted data base
C                and puts them in MATLAB .MAT format
C
C***********************************************************************
#include "mlxin.h"
      integer mxgas, lun
      real*8 mxstr
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
       INTEGER index9,index6,index3,linct3
       real*4             STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP
       DOUBLE PRECISION WNUM
       CHARACTER*3      AI
       CHARACTER*6      REF
       CHARACTER*9      USLQ,BSLQ
C
       real*8       WNUMLO,WNUMHI,GSTR(MXGAS),STR, linct2
       INTEGER    KFND,ISUM,I,IDG(MXGAS),ICH(MXGAS),LINCT
      real*8   idgo(MXGAS), icho(MXGAS)
       CHARACTER  TITLB*48,INFILE*80,OUTFILE*80
C
C  DECLARE ARRAYS TO HOLD LINE DATA FOR MATLAB
C

c      MXLIN is really a dummy argument here. The space has already
c      been allocated by matlab. Fortran doesn't mind if you
c      over write your subscripts as long as you don't write
c      where you shouldn't.
       real*8 WNUMHX,ZLSTAT(MXLIN),ZIGAS,ZISO(MXLIN),
     +                  ZWNUM(MXLIN),ZSTREN(MXLIN),ZTPROB(MXLIN),
     +                  ZABROAD(MXLIN),ZSBROAD(MXLIN),ZELS(MXLIN),
     +                  ZABCOEF(MXLIN),ZTSP(MXLIN),ZIUSGQ(MXLIN),
     +                  ZILSGQ(MXLIN),ZUSLQ(MXLIN),ZBSLQ(MXLIN),
     +                  ZAI(MXLIN),ZREF(MXLIN),isum2, GASID(MXLIN)

C
C  Set i/o unit numbers
C

c      In matlab everything is real*8. some things need to be placed
c       in an integer for Fortran      

      isum=int(isum2)
      do i=1,isum
            idg(i)= int(idgo(i))
            ich(i)=int(icho(i))
      end do
        LUN = 12
      index9=0
      index6=0
      index3=0
C
C  Read minimum frequency
C

c      old read write commands left over from when this was a stand alone
c      binary (hittomat.f)

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
C
C  Read in gas id numbers
       DO 10, I=1,ISUM
c         WRITE(6,*) 'Enter id of gas:'
c         READ(5,*) IDG(I)
c         WRITE(6,*) 'Enter line version #:'
c         READ(5,*) ICH(I)
c         WRITE(6,*) 'Enter minimum strength:'
c         READ(5,*) STR
         GSTR(I) = STR * 6.022E26
10     CONTINUE
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
C  Save valid entry data for gas id. Note for matlab everything must be real*8.
C

           LINCT=LINCT+1
         GASID(LINCT)=dble(IGAS)
           ZLSTAT(LINCT)=dble(LSTAT)
           ZISO(LINCT)=dble(ISO)
           ZWNUM(LINCT)=dble(WNUM)
           ZSTREN(LINCT) = dble(STREN/6.022E26)
           ZTPROB(LINCT)=dble(TPROB)
           ZABROAD(LINCT)=dble(ABROAD)
           ZSBROAD(LINCT)=dble(SBROAD)
           ZELS(LINCT)=dble(ELS)
           ZABCOEF(LINCT)=dble(ABCOEF)
           ZTSP(LINCT)=dble(TSP)
           ZIUSGQ(LINCT)=IUSGQ
           ZILSGQ(LINCT)=ILSGQ
c      note, these strings are turned into integers so they can be easily
c      read into matlab. Once read into matlab they will be converted back
c      to characters
           DO 300 I=1,9
                 index9=index9+1
              ZUSLQ(index9)=dble(ICHAR(USLQ(I:I)))
              ZBSLQ(index9)=dble(ICHAR(BSLQ(I:I)))
 300        CONTINUE
           DO 310 I=1,3
              index3=index3+1
              ZAI(index3)=dble(ICHAR(AI(I:I)))
 310        CONTINUE
           DO 320 I=1,6
                 index6=index6+1
              ZREF(index6)=dble(ICHAR(REF(I:I)))
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

c      convert linct back to a double to be read by matlab.      

      linct2=dble(linct)
C       STOP
       END

      subroutine hittomatd(WNUMLO, WNUMHI,INFILE,idgo,icho,str,
     +       LINCT, isum2)
c this routine is exactly like hittomat except it doesn't save any data/
c it is used to find how much memory we need to allocate. 
C***********************************************************************
C
C  Program       HITTOMAT
C
C  Purpose       Program copys lines from HITRAN unformatted data base
C                and puts them in MATLAB .MAT format
C
C***********************************************************************
#include "mlxin.h"
      integer mxgas, lun
      real*8 mxstr
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
       real*4             STREN,TPROB,ABROAD,SBROAD,ELS,ABCOEF,TSP
       DOUBLE PRECISION WNUM
       CHARACTER*3      AI
       CHARACTER*6      REF
       CHARACTER*9      USLQ,BSLQ
C
       real*8       WNUMLO,WNUMHI,GSTR(MXGAS),STR, linct2
       INTEGER    KFND,ISUM,I,IDG(MXGAS),ICH(MXGAS),LINCT
      real*8   idgo(mxgas), icho(mxgas)
       CHARACTER  TITLB*48,INFILE*80,OUTFILE*80
       real*8 WNUMHX,ZLSTAT(MXLIN),ZIGAS,ZISO(MXLIN),
     +                  ZWNUM(MXLIN),ZSTREN(MXLIN),ZTPROB(MXLIN),
     +                  ZABROAD(MXLIN),ZSBROAD(MXLIN),ZELS(MXLIN),
     +                  ZABCOEF(MXLIN),ZTSP(MXLIN),ZIUSGQ(MXLIN),
     +                  ZILSGQ(MXLIN),ZUSLQ(MXLIN,9),ZBSLQ(MXLIN,9),
     +                  ZAI(MXLIN,3),ZREF(MXLIN,6),isum2

C
C  DECLARE ARRAYS TO HOLD LINE DATA FOR MATLAB
C
C
C  Set i/o unit numbers
C
        isum=int(isum2)
      do i=1,isum
            idg(i)= int(idgo(i))
            ich(i)=int(icho(i))
      end do
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
c     WRITE(6,*) 'Enter number of gases to select:'
c     READ(5,*) ISUM
c       ISUM=1
C
c  Read in gas id numbers
       DO 10, I=1,ISUM
c         WRITE(6,*) 'Enter id of gas:'
c         READ(5,*) IDG(I)
c         WRITE(6,*) 'Enter line version #:'
c         READ(5,*) ICH(I)
c         WRITE(6,*) 'Enter minimum strength:'
c         READ(5,*) STR
         GSTR(I) = STR * 6.022E26
 10     CONTINUE
C
C  Read in hitran linebase filename
C

c       WRITE(6,*) 'Enter name of hitran linebase file:'
c       READ(5,120) INFILE
c120    FORMAT (A)

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
         END IF
C
C  Read another entry from file
C
         GOTO 200
       END IF
C
C      No need to keep all gas id's since this program reads in only one
C
C
C  Print # of lines selected
C
c       WRITE(6,399) LINCT
 399   FORMAT('Number of lines selected: ',I6)
C
C  Write out selected data
C
C

       END




