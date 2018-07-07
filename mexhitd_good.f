#include <fintrf.h>	
	SUBROUTINE mexFunction(NLHS, PLHS, NRHS, PRHS)
	POINTER PLHS(*), PRHS(*)
	INTEGER NLHS, NRHS
	INTEGER M, N, SIZE, XP, YP,linct2
	INTEGER MXGETM, MXGETN, status, mxgetstring
	POINTER MXGETPR, MXCREATEFULL

	pointer LINCT,WNUMHX, ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB
	POINTER ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,ZILSGQ
	POINTER ZUSLQ,ZBSLQ,ZAI,ZREF

	POINTER WNUMLO,WNUMHI,idg,ich,str, isum
	CHARACTER*80 INFILE

	WNUMLO = MXGETPR(PRHS(1))
	WNUMHI = MXGETPR(PRHS(2))
        m= mxgetm(prhs(3))
        n= mxgetn(prhs(3))
        size=max(m,n)
        status=mxgetstring(prhs(3),INFILE,size)
	idg = MXGETPR(PRHS(4))
	ich = MXGETPR(PRHS(5))
	str = MXGETPR(PRHS(6))
	isum= MXGETPR(PRHS(7))
      call hittomatd(%val(WNUMLO), %val(WNUMHI),INFILE,%val(idg),
     +       %val(ich),%val(str),linct2, %val(isum))


	print*,'linct2',linct2

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

	PLHS(16) = MXCREATEFULL(LINCT2,9,0)
	PLHS(17) = MXCREATEFULL(LINCT2,9,0)
	PLHS(18) = MXCREATEFULL(LINCT2,3,0)
	PLHS(19) = MXCREATEFULL(LINCT2,6,0)

c	PLHS(16) = MXCREATEFULL(9*LINCT2,0)
c	PLHS(17) = MXCREATEFULL(9*LINCT2,0)
c	PLHS(18) = MXCREATEFULL(3*LINCT2,0)
c	PLHS(19) = MXCREATEFULL(6*LINCT2,0)

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
      call hittomat(%val(WNUMLO), %val(WNUMHI),INFILE,%val(idg),
     +       %val(ich),%val(str),
c outputs next
     +       %val(LINCT),%val(WNUMHX),%val(ZLSTAT),%val(ZIGAS),
     +       %val(ZISO),%val(ZWNUM),%val(ZSTREN),%val(ZTPROB),
     +       %val(ZABROAD),%val(ZSBROAD),%val(ZELS),%val(ZABCOEF),
     +       %val(ZTSP),%val(ZIUSGQ),
     +       %val(ZILSGQ),%val(ZUSLQ),%val(ZBSLQ),%val(ZAI),
     +       %val(ZREF),linct2,%val(isum))
	END
      subroutine hittomat(WNUMLO, WNUMHI,INFILE,idgo,icho,str,
     +       LINCT2,WNUMHX,ZLSTAT,ZIGAS,ZISO,ZWNUM,ZSTREN,ZTPROB,
     +       ZABROAD,ZSBROAD,ZELS,ZABCOEF,ZTSP,ZIUSGQ,
     +       ZILSGQ,ZUSLQ,ZBSLQ,ZAI,ZREF,linct3, isum2)

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
       real*8 WNUMHX,ZLSTAT(MXLIN),ZIGAS,ZISO(MXLIN),
     +                  ZWNUM(MXLIN),ZSTREN(MXLIN),ZTPROB(MXLIN),
     +                  ZABROAD(MXLIN),ZSBROAD(MXLIN),ZELS(MXLIN),
     +                  ZABCOEF(MXLIN),ZTSP(MXLIN),ZIUSGQ(MXLIN),
     +                  ZILSGQ(MXLIN),ZUSLQ(9,MXLIN),ZBSLQ(9,MXLIN),
     +                  ZAI(3,MXLIN),ZREF(6,MXLIN),isum2

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
C  Save valid entry data for gas id
C
           LINCT=LINCT+1
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
           DO 300 I=1,9
           	index9=index9+1
              ZUSLQ(linct,i)=ICHAR(USLQ(I:I))
              ZBSLQ(linct,i)=ICHAR(BSLQ(I:I))
 300        CONTINUE
           DO 310 I=1,3
              index3=index3+1
              ZAI(linct,i)=ICHAR(AI(I:I))
 310        CONTINUE
           DO 320 I=1,6
           	index6=index6+1
              ZREF(linct,i)=ICHAR(REF(I:I))
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
	linct2=dble(linct)
C       STOP
       END

      subroutine hittomatd(WNUMLO, WNUMHI,INFILE,idgo,icho,str,
     +       LINCT, isum2)

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




