c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

cccccc these are the original versions, copied from GENLN2 in mid 1999

c************************************************************************
ccccccccccc include the common data blocks
c************************************************************************
      include 'old_n2ct0.f'
c************************************************************************
      include 'old_o2ctt.f'
c************************************************************************

c************************************************************************
c                   OLD NITROGEN CONTINUUM
c************************************************************************
c this is OLD N2 cont                pre April 2000
c copied from /salsify/users/sergio/KCARTA/SPECTRA/Genln2/new_contum3.f

       SUBROUTINE CALCONNIT( CON, IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
     $    PARTP, AMT, whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
       REAL*8 VB1,VT1,DV1,NPT1,H2OS96(2001)
       REAL*8 H2OS60(2001), H2OF96(2001),VB2,VT2,DV2,NPT2,CO2F96(2001)
       REAL*8 CO2F60(2001),CO2F30(2001)
       REAL*8 VB7,VT7,DV7,NPT7,O2C213(59),O2C253(59),O2C293(59)
       REAL*8 VB22,VT22,DV22,NPT22,CN2297(111),CN2273(111)
       REAL*8 CN2253(111),CN2233(111),CN2213(111),CN2193(111)

       COMMON /CH2OS0/ VB1,VT1,DV1,NPT1,H2OS96
       COMMON /CH2OS1/ H2OS60
       COMMON /CH2OF0/ H2OF96
       COMMON /CCO2F0/ VB2,VT2,DV2,NPT2,CO2F96
       COMMON /CCO2F1/ CO2F60
       COMMON /CCO2F2/ CO2F30
       COMMON /CO2CT/  VB7,VT7,DV7,NPT7,O2C213,O2C253,O2C293
C
       COMMON /CN2C0/ VB22,VT22,DV22,NPT22,CN2297,CN2273,
     + CN2253,CN2233,CN2213,CN2193

       REAL*8 CFFF,CFFFH,CFFFL,CE,CEH,CEL,CE297,CE273,CE253,CE233,CE213,
     + CE193,TH,TL,CFF297,CFF273,CFF253,CFF233,CFF213,CFF193,CPART,X,
     + C11(kMaxLayer)
C

C      Arguements
       INTEGER IDGAS, NFREQ, NLAY, whichlayer
       REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*)

C      Variables for new water continuum adjustment terms
       INTEGER JFAC
       REAL*8 ALPHS2, BETAS, V0S, FACTRS, HWSQF, BETAF, V0F, FACTRF,
     $    V0F2, HWSQF2, BETA2, SFAC, alpha2, scor, FSCAL

       INTEGER IPT,ILAY

       REAL*8 C2, AVOG, TS, PS
       DATA C2/1.4387863/
       DATA AVOG/6.022045E+26/
       DATA TS/296.0/
       DATA PS/1.0/

c used in calculating d/dq, d/dT
       REAL*8 ra1,ra2,ra3,rQ
       REAL*8 TFAC(kMaxLayer),XAMT(kMaxLayer)

      REAL*8 vb,vt,dv,vs2,vf2,vf6,vf4,csfft
c      REAL*8 csff60(MaxLen),csff96(MaxLen),cfff96(MaxLen)
      REAL*8 csff60,csff96,cfff96
      REAL*8 f1,gr,v0f3,hwsqf3,beta3
      INTEGER iL,mk,npt,iF,iFrLow,iFrHigh

      REAL*8 ALL,FF,A1,A2,A3,A4
  
C-----------------------------------------------------------------------

       IF (IDGAS .EQ. 22)THEN
         VB = VB22
         VT = VT22
         DV = DV22
         NPT = NPT22
         ENDIF

       IF (IDGAS .EQ. 22) THEN
         DO 10 ILAY=whichlayer,whichlayer
           C11(ILAY) = 0.5*C2/T(ILAY)
           ALL = AMT(ILAY)*AVOG
C
           DO 40 IPT=1,NFREQ
             FF=FREQ(IPT)
             IF (FF .GE. VB .AND. FF .LT. VT) THEN
               MK = 1 + INT((FF - VB)/DV) 
               F1 = VB + (MK - 1)*DV

               GR = (FF - F1)/DV
                 CFF297 = CN2297(MK) + GR*(CN2297(MK+1)-CN2297(MK))
                 CFF273 = CN2273(MK) + GR*(CN2273(MK+1)-CN2273(MK))
                 CFF253 = CN2253(MK) + GR*(CN2253(MK+1)-CN2253(MK))
                 CFF233 = CN2233(MK) + GR*(CN2233(MK+1)-CN2233(MK))
                 CFF213 = CN2213(MK) + GR*(CN2213(MK+1)-CN2213(MK))
                 CFF193 = CN2193(MK) + GR*(CN2193(MK+1)-CN2193(MK))
C
C                Temp dependence of collision efficiency O2 vs N2
                 CE297=0.83
                 CE273=0.90
                 CE253=0.93
                 CE233=0.94
                 CE213=0.94
                 CE193=1.00
C
                 IF (T(ILAY) .GE. 273.0) THEN
                    CFFFH=CFF297
                    CFFFL=CFF273
                    CEH=CE297
                    CEL=CE273
                    TH=297.0
                    TL=273.0
                 ELSEIF (T(ILAY) .GE. 253.0 .AND. T(ILAY) .LT.
     +           273.0) THEN
                    CFFFH=CFF273
                    CFFFL=CFF253
                    CEH=CE273
                    CEL=CE253
                    TH=273.0
                    TL=253.0
                 ELSEIF (T(ILAY) .GE. 233.0 .AND. T(ILAY) .LT.
     +           253.0) THEN
                    CFFFH=CFF253
                    CFFFL=CFF233
                    CEH=CE253
                    CEL=CE233
                    TH=253.0
                    TL=233.0
                 ELSEIF (T(ILAY) .GT. 213.0 .AND. T(ILAY) .LT.
     +           233.0) THEN
                    CFFFH=CFF233
                    CFFFL=CFF213
                    CEH=CE233
                    CEL=CE213
                    TH=233.0
                    TL=213.0
                 ELSEIF (T(ILAY) .LE. 213.0) THEN
                    CFFFH=CFF213
                    CFFFL=CFF193
                    CEH=CE213
                    CEL=CE193
                    TH=213.0
                    TL=193.0
                 ENDIF
C
C                Linear interpolations in temperature
                 CFFF=CFFFL+(CFFFH-CFFFL)*(T(ILAY)-TL)/(TH-TL)
                 CE=CEL+(CEH-CEL)*(T(ILAY)-TL)/(TH-TL)
C
C                Note: 273.15^2=74610.922
                 CPART=( PARTP(ILAY)*74610.922/(T(ILAY)*T(ILAY)) )*
     +           CFFF*( PARTP(ILAY)+CE*( P(ILAY)-PARTP(ILAY) ) )
C
C                Note: AVOG/(2.6867E+19*273.15) = 82058.508
                 X=82058.508*( AMT(ILAY)*T(ILAY) )/PARTP(ILAY)
C
                 A4 = CPART*X*1.00E-06
ccc
C                note: this only adjusts the entire con, not self/foreign
cccccccccccc                 if (xgas .eq. 22) a4=a4*xcon1
ccc

               CON(IPT)=a4

               END IF
   40        CONTINUE
   10     CONTINUE
       ENDIF

       RETURN
       END

c************************************************************************
c                   OLD OXYGEN CONTINUUM
c************************************************************************
c this is OLD OXYGEN cont       pre April 2000
c copied from /salsify/users/sergio/KCARTA/SPECTRA/Genln2/new_contum3.f

c this subroutine calculates the continuum for O2(ID=7)
       SUBROUTINE CALCONOXY( CON, IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
     $    PARTP, AMT, whichlayer)

      include '../FORTRANFILES/max.inc'

C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
       REAL*8 VB1,VT1,DV1,NPT1,H2OS96(2001)
       REAL*8 H2OS60(2001), H2OF96(2001),VB2,VT2,DV2,NPT2,CO2F96(2001)
       REAL*8 CO2F60(2001),CO2F30(2001)
       REAL*8 VB7,VT7,DV7,NPT7,O2C213(59),O2C253(59),O2C293(59)
       REAL*8 VB22,VT22,DV22,NPT22,CN2297(111),CN2273(111)
       REAL*8 CN2253(111),CN2233(111),CN2213(111),CN2193(111)

       COMMON /CH2OS0/ VB1,VT1,DV1,NPT1,H2OS96
       COMMON /CH2OS1/ H2OS60
       COMMON /CH2OF0/ H2OF96
       COMMON /CCO2F0/ VB2,VT2,DV2,NPT2,CO2F96
       COMMON /CCO2F1/ CO2F60
       COMMON /CCO2F2/ CO2F30
       COMMON /CO2CT/  VB7,VT7,DV7,NPT7,O2C213,O2C253,O2C293
C
       COMMON /CN2C0/ VB22,VT22,DV22,NPT22,CN2297,CN2273,
     + CN2253,CN2233,CN2213,CN2193

       REAL*8 CFFF,CFFFH,CFFFL,CE,CEH,CEL,CE297,CE273,CE253,CE233,CE213,
     + CE193,TH,TL,CFF297,CFF273,CFF253,CFF233,CFF213,CFF193,CPART,X,
     + C11(kMaxLayer)
C

C      Arguements
       INTEGER IDGAS, NFREQ, NLAY, whichlayer
       REAL*8 FREQ(*), FSTEP, T(*), P(*),PARTP(*),AMT(*),CON(*)

C      Variables for new water continuum adjustment terms
       INTEGER JFAC
       REAL*8 ALPHS2, BETAS, V0S, FACTRS, HWSQF, BETAF, V0F, FACTRF,
     $    V0F2, HWSQF2, BETA2, SFAC, alpha2, scor, FSCAL

       INTEGER IPT,ILAY

       REAL*8 C2, AVOG, TS, PS
       DATA C2/1.4387863/
       DATA AVOG/6.022045E+26/
       DATA TS/296.0/
       DATA PS/1.0/

c used in calculating d/dq, d/dT
       REAL*8 ra1,ra2,ra3,rQ
       REAL*8 TFAC(kMaxLayer),XAMT(kMaxLayer)

      REAL*8 vb,vt,dv,vs2,vf2,vf6,vf4,csfft
c      REAL*8 csff60(MaxLen),csff96(MaxLen),cfff96(MaxLen)
      REAL*8 csff60,csff96,cfff96
      REAL*8 f1,gr,v0f3,hwsqf3,beta3
      INTEGER iL,mk,npt,iF,iFrLow,iFrHigh

      REAL*8 ALL,FF,A1,A2,A3,A4

      REAL*8 CFFF13,CFFF53,CFFF93,CC2,CC1,CC0,CFFFT
C-----------------------------------------------------------------------

       IF (IDGAS .EQ. 7)THEN
         VB = VB7
         VT = VT7
         DV = DV7
         NPT = NPT7
         ENDIF

       IF (IDGAS .EQ. 7) THEN
         DO 10 ILAY=whichlayer,whichlayer
           C11(ILAY) = 0.5*C2/T(ILAY)
           ALL = AMT(ILAY)*AVOG
C
           DO 40 IPT=1,NFREQ
             FF=FREQ(IPT)
             IF (FF .GE. VB .AND. FF .LT. VT) THEN
               MK = 1 + INT((FF - VB)/DV) 
               F1 = VB + (MK - 1)*DV

               GR = (FF - F1)/DV
               CFFF13 = O2C213(MK) + GR*(O2C213(MK+1)-O2C213(MK))
               CFFF53 = O2C253(MK) + GR*(O2C253(MK+1)-O2C253(MK))
               CFFF93 = O2C293(MK) + GR*(O2C293(MK+1)-O2C293(MK))
C
               CC2 = (CFFF93 - 2.0*CFFF53 + CFFF13)/3.2E+03
               CC1 = (CFFF53 - CFFF13 - 1.864E+04*CC2)/4.0E+01
               CC0 = CFFF13 - 2.13E+02*CC1 - 4.5369E+04*CC2
C
               CFFFT = CC0 + CC1*T(ILAY) + CC2*T(ILAY)**2
C
C  COEFFICIENTS ARE FOR TOTAL MASS PATH NOT O2 MASS PATH (TOTAL=O2/0.22) 
C  CONV FACTOR = R*E4/(1.01325E5*0.22) = 3.729672
C
               A1 = 3.729672*P(ILAY)*T(ILAY)*AMT(ILAY)*CFFFT

               a4=a1

c               print *,ipt,a4

               CON(IPT)=a4

               END IF
   40        CONTINUE
   10     CONTINUE
       ENDIF

       RETURN
       END

c************************************************************************
