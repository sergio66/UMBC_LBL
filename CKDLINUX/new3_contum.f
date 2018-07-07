C***********************************************************************
C
C  PROGRAM        CONTUM   SUBROUTINE
C
C  PURPOSE        TO CALCULATE THE CONTINUUM ABSORPTION
C
C  VERSION        3.Y   D.P. EDWARDS   24/02/93
c Scott, June 1995. modified version of "new_contum" for Clough's
c version 2.1 continuum.
C                 Scott Hannon, November 1993
c modified by Scott to allow multiplying the (self/foreign) continuum.
c requires that genln2 and input also be modified (both for common block
c and only input for keyword xcontm). (program previously called xcontum)
c
c Scott, April 1998.  Modified new2_contum for CKD2.3.  Note that for
c correct use, it it necessary to use the *old* h2ost[01] and h2oft0
c file that *include* the "basement" term.
C 
C  DESCRIPTION    THE CONTINUUM ABSORPTION IS CALCULATED FOR H2O, CO2, 
C                 O2 AND N2 AT THE FREQUENCY FF CM-1.
C
C  REFS.          CLOUGH S.A.,KNEIZYS F.X.,DAVIES R.,GAMACHE R., 
C                 TIPPING R. (1980)
C                 Theoretical line shape for H2O vapour: 
C                 Application to the continuum 
C                 Atmospheric Water Vapour, Eds. A.Deepak,T.D.Wilkeson, 
C                 L.H.Ruhnke, Academic Press, New York
C
C                 CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,GALLERY W.O. 
C                 (1981), Atmospheric spectral transmittance and 
C                 radiance: FASCOD1B
C                 SPIE 277 Atmospheric Transmission  152-166
C
C                 CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,ANDERSON G.P., 
C                 SHETTLE E.P. (1987)
C                 Current issues in infrared atmospheric transparency
C                 International meeting on
C                 Atmospheric Transparency for Satellite Applications
C                 15-19 Sept. 1986 Capri, Italy. Ed. G.V. Silvestrini.
C                 CUEN.
C
C                 TIMOFEYEV, Y.M. AND M.V. TONKOV
C                 EFFECT OF THE INDUCED OXYGEN ABSORPTION BAND ON THE
C                 TRANSFORMATION OF RADIATION IN THE 6UM REGION OF THE
C                 EARTH'S ATMOSPHERE, 
C                 IZV.ACAD.SCI. USSR
C                 ATMOS.OCEAN.PHYS., ENGLISH, 14, 437-441, 1978.
C
C  ARGUMENT       NPATH  I*4  I/P  NO. OF PATHS FOR CALCULATION   
C                 IFLAG  I*4  I/P  WIDE MESH CALCULATION SPECIFIER
C                 WABS   R*4  I/P  WIDEPASS ABSORPTION BY 
C                                  PAROBOLIC POINT, WIDE MESH AND PATH 
C                 CABS   R*4  I/P WIDEPASS NONLTE FACTOR ABSORPTION BY
C                                 PARABOLIC POINT, WIDE MESH, AND PATH
C
C  SUBROUTINES    GENDAT, CO2FT0, CO2FT1, CO2FT2, H2OFT0, H2OST0, 
C                 H2OST1, N2CT93,  O2CTT
C
C***********************************************************************
C
       SUBROUTINE CONTUM(NPATH,IFLAG,WABS,CABS)
C-----------------------------------------------------------------------
#include "parray.f"
C
ccc
c common block for continuum change
       common /xcontm/ xgas,xcon1,xcon2
       integer xgas
       real    xcon1,xcon2
ccc
C  COMMON BLOCK CONTAINING PATH GAS PARAMETERS
C
       COMMON /COMPTH/ 
     + IDFIL(MXPTH),IDGAS(MXPTH),IDISO(MXPTH),AMT(MXPTH),
     + T(MXPTH),DELT(MXPTH),P(MXPTH),DELP(MXPTH),PARTP(MXPTH),
     + VEL(MXPTH),LINSHP(MXPTH),LCONT(MXPTH)
C
C  COMMON BLOCK CONTAINING FREQUENCY MESH PARAMETERS
C
       COMMON /COMFRQ/ 
     + FBDY(MXINT),FWIND,FEXC,NPL,NPU,MINT,NDIV,NWME,NFME(MXINT),
     + NWID
C
C  COMMON BLOCKS CONTAINING CONTINUUM VALUES 
C
       COMMON /CH2OS0/ VB1,VT1,DV1,NPT1,H2OS96(2001)
       COMMON /CH2OS1/ H2OS60(2001)
       COMMON /CH2OF0/ H2OF96(2001)
       COMMON /CCO2F0/ VB2,VT2,DV2,NPT2,CO2F96(2001)
       COMMON /CCO2F1/ CO2F60(2001)
       COMMON /CCO2F2/ CO2F30(2001)
       COMMON /CO2CT/  VB7,VT7,DV7,NPT7,O2C213(59),O2C253(59),O2C293(59)
C
       COMMON /CN2C0/ VB22,VT22,DV22,NPT22,CN2297(111),CN2273(111),
     + CN2253(111),CN2233(111),CN2213(111),CN2193(111)
C
C  COMMON BLOCK CONTAINING LINE SHAPE INTEGER CODES
C
       COMMON /COMLIN/ 
     + LVOIGT,LLORTZ,LVANVH,LNEWSH,LWSCON,LDOPLR,LVOSUB,LWFCON,
     + LXSECN,LNOCON,LWTHCN
C
C  COMMON BLOCK COTAINING PHYSICAL CONSTANTS
C
       COMMON/COMCON/ 
     + TS,PS,C1,C2,VLIGHT,AVOG,R2,PI,GRAV,CPAIR,ATMB
C
       INTEGER   IFLAG(MXINT)
       REAL      WABS(3,MXINT,MXPTH),CABS(3,MXINT,MXPTH)
C
       REAL CFFF,CFFFH,CFFFL,CE,CEH,CEL,CE297,CE273,CE253,CE233,CE213,
     + CE193,TH,TL,CFF297,CFF273,CFF253,CFF233,CFF213,CFF193,CPART,X
C
ccc for new continuum adjustment terms
       INTEGER JFAC
       REAL ALPHS2, BETAS, V0S, FACTRS, HWSQF, BETAF, V0F, FACTRF,
     $    V0F2, HWSQF2, BETA2, SFAC, alpha2, scor, FSCAL,
     $    V0F3, HWSQF3, BETA3
       DIMENSION XFAC(0:50)
       DATA (XFAC(I),I=0,50)/
     $    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
     $    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
     $    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
     $    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
     $    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
     $    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
     $    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
     $    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
     $    1.00000,1.00000,1.00000/
ccc

C-----------------------------------------------------------------------
C
       DO 10 IPATH=1,NPATH
        IF (LCONT(IPATH) .EQ. LWTHCN) THEN
         IGS = IDGAS(IPATH)
C
         IF (IGS .EQ. 1)THEN
           VB = VB1
           VT = VT1
           DV = DV1
           NPT = NPT1
ccc terms for new water continuum adjustment factors
                 TFAC =  2.77777E-2*(296.0 - T(IPATH)) 
      ALPHA2 = 200.0**2
      ALPHS2 = 120.0**2
      BETAS = 5.0E-06
      V0S = 1310.0
      FACTRS = 0.15
C
      HWSQF = 330.0**2
      BETAF = 8.0E-11
      V0F = 1130.0
      FACTRF = 0.97
C
      V0F2 = 1900.0
      HWSQF2 = 150.0**2
      BETA2 = 3.0E-06
C
      V0F3 = 1596.0
      HWSQF3 = 150.0**2
      BETA3 = 3.0E-6
C
ccc
         ELSEIF (IGS .EQ. 2) THEN
           VB = VB2
           VT = VT2
           DV = DV2
           NPT = NPT2
         ELSEIF (IGS .EQ. 7) THEN
           VB = VB7
           VT = VT7
           DV = DV7
           NPT = NPT7
         ELSEIF (IGS .EQ. 22) THEN
           VB = VB22
           VT = VT22
           DV = DV22
           NPT = NPT22
         ELSE 
           GOTO 10
         ENDIF
C
         C11 = 0.5*C2/T(IPATH)
         ALL = AMT(IPATH)*AVOG
C
C  INTERPOLATE CONTINUUM DATA POINTS AT WIDEMESH FREQUENCY POINTS
C
         DO 30 IPW=1,NWME
C
C  CHECK IF CALCULATIONS ARE BEING PERFORMED FOR THIS WIDE MESH
C
          IF (IFLAG(IPW) .EQ. 1) THEN
C
           MLOW = NPL + IPW - 1
           DO 40 IP=1,3
             FF = FBDY(MLOW) +
     +       FLOAT(IP-1)*0.5*(FBDY(MLOW+1) - FBDY(MLOW))  
C
             IF (FF .GE. VB .AND. FF .LT. VT) THEN
               MK = 1 + INT((FF - VB)/DV) 
               F1 = VB + (MK - 1)*DV
C
C              ==========================
C              H2O WATER VAPOUR CONTINUUM
C              ==========================
               IF (IGS .EQ. 1) THEN
                 GR = (FF - F1)/DV 
                 CSFF60 = H2OS60(MK) + GR*(H2OS60(MK+1)-H2OS60(MK))       
                 CSFF96 = H2OS96(MK) + GR*(H2OS96(MK+1)-H2OS96(MK))
                 CFFF96 = H2OF96(MK) + GR*(H2OF96(MK+1)-H2OF96(MK))
C
C                ------------------------------------------
C                Water Self Continuum adjustment/correction
C
C                Adjust for temperature
                 CSFFT = CSFF96*(CSFF60/CSFF96)**TFAC
C
C                Make correction
                 vs2=(FF-v0s)**2
C
                 sfac=1.0
                 if (FF .GE. 700.0 .AND. FF .LE. 1200.0) THEN
                    jfac=(FF-700.0)/10.0 + 0.0001
                    sfac=xfac(jfac)
                 endif
c
                 SCOR = sfac*
     $           (1.0 + 2.02*(1.0E+4/(FF**2 + 1.0E-4*FF**4 + 1.0E+4)))*
     $           (1.0 - 0.2333*(ALPHA2/((FF - 1050.0)**2 + ALPHA2)))*
     $           (1.0 - factrs*(alphs2/(vs2 + (betas*vs2**2) + alphs2)))
c
                 CSFFT = CSFFT* SCOR
C
C                ---------------------------------------------
C                Water Foreign Continuum adjustment/correction
C Note: there is a new water foreign con for 0-800 cm-1 that should
C put in here.
C
                 vf2=(FF-v0f)**2
                 vf6=vf2*vf2*vf2
                 fscal=(1.0-factrf*(hwsqf/(vf2 + (betaf*vf6) +
     $              hwsqf)))
c
                 vf2=(FF-v0f2)**2
                 vf4=vf2*vf2
                 fscal=fscal*(1.0-0.6*(hwsqf2/(vf2 + beta2*vf4 +
     $              hwsqf2)))
c
                 vf2=(FF-v0f3)**2
                 vf4=vf2*vf2
                 fscal=fscal*(1.0-0.2*(hwsqf3/(vf2 + beta3*vf4 +
     $              hwsqf3)))
c
                 CFFF96=CFFF96*fscal
C
C                ------------------------------
C                Add foreign and self continuum
ccc
                 if (xgas .eq. 1) then
                    CSFFT = CSFFT*xcon1
                    CFFF96 = CFFF96*xcon2
                 endif
ccc
cc
                 if (linshp(ipath) .EQ. LWSCON) then
c                   self only; foreign to zero
                    CFFF96=0.0
                 elseif (linshp(ipath) .EQ. LWFCON) then
c                   foreign only; set self to zero
                    CSFFT=0.0
                 endif
cc
                 A1 = FF*ALL*TANH(C11*FF)
                 A2 = TS/(T(IPATH)*PS)
                 A3 = 1.0E-20*(PARTP(IPATH)*CSFFT + 
     +                (P(IPATH) - PARTP(IPATH))*CFFF96)
                 A4 = A1*A2*A3
                 WABS(IP,IPW,IPATH)  = WABS(IP,IPW,IPATH) + A4
                 CABS(IP,IPW,IPATH)  = CABS(IP,IPW,IPATH) + A4
C
C              =====================
C              CO2 FOREIGN CONTINUUM
C              =====================
               ELSEIF (IGS .EQ. 2) THEN
C
C  LINEARLY INTERPOLATE IN FREQUENCY AT EACH TEMPERATURE
C
                 GR = (FF - F1)/DV
                 CFFF96 = CO2F96(MK) + GR*(CO2F96(MK+1)-CO2F96(MK))
                 CFFF60 = CO2F60(MK) + GR*(CO2F60(MK+1)-CO2F60(MK))
                 CFFF30 = CO2F30(MK) + GR*(CO2F30(MK+1)-CO2F30(MK))
C
C  LAGRANGIAN INTERPOLATION IN TEMPERATURE
C
                 DT30 = T(IPATH) - 230.
                 DT60 = T(IPATH) - 260.
                 DT96 = T(IPATH) - 296.
                 CFFFT = 5.050505E-4*DT60*DT96*CFFF30
     +           - 9.259259E-4*DT30*DT96*CFFF60
     +           + 4.208754E-4*DT30*DT60*CFFF96
C
                 A4 = AMT(IPATH)*(P(IPATH)/PS)*CFFFT 
ccc
                if (xgas .eq. 2) a4=a4*xcon1
ccc
C
                 WABS(IP,IPW,IPATH)  = WABS(IP,IPW,IPATH) + A4
                 CABS(IP,IPW,IPATH)  = CABS(IP,IPW,IPATH) + A4
C
C              ==============================
C              O2 COLLISION INDUCED CONTINUUM
C              ==============================
               ELSEIF (IGS .EQ. 7) THEN
C
                 GR = (FF - F1)/DV
                 CFFF13 = O2C213(MK) + GR*(O2C213(MK+1)-O2C213(MK))
                 CFFF53 = O2C253(MK) + GR*(O2C253(MK+1)-O2C253(MK))
                 CFFF93 = O2C293(MK) + GR*(O2C293(MK+1)-O2C293(MK))
C
                 CC2 = (CFFF93 - 2.0*CFFF53 + CFFF13)/3.2E+03
                 CC1 = (CFFF53 - CFFF13 - 1.864E+04*CC2)/4.0E+01
                 CC0 = CFFF13 - 2.13E+02*CC1 - 4.5369E+04*CC2
C
                 CFFFT = CC0 + CC1*T(IPATH) + CC2*T(IPATH)**2
C
C  COEFFICIENTS ARE FOR TOTAL MASS PATH NOT O2 MASS PATH (TOTAL=O2/0.22) 
C  CONV FACTOR = R*E4/(1.01325E5*0.22) = 3.729672
C
                 A1 = 3.729672*P(IPATH)*T(IPATH)*AMT(IPATH)*CFFFT
ccc
                 if (xgas .eq. 7) a1=a1*xcon1
ccc
                 WABS(IP,IPW,IPATH)  = WABS(IP,IPW,IPATH) + A1
                 CABS(IP,IPW,IPATH)  = CABS(IP,IPW,IPATH) + A1
C
C              ==============================
C              N2 COLLISION INDUCED CONTINUUM
C              ==============================
               ELSEIF (IGS .EQ. 22) THEN
C
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
                 IF (T(IPATH) .GE. 273.0) THEN
                    CFFFH=CFF297
                    CFFFL=CFF273
                    CEH=CE297
                    CEL=CE273
                    TH=297.0
                    TL=273.0
                 ELSEIF (T(IPATH) .GE. 253.0 .AND. T(IPATH) .LT.
     +           273.0) THEN
                    CFFFH=CFF273
                    CFFFL=CFF253
                    CEH=CE273
                    CEL=CE253
                    TH=273.0
                    TL=253.0
                 ELSEIF (T(IPATH) .GE. 233.0 .AND. T(IPATH) .LT.
     +           253.0) THEN
                    CFFFH=CFF253
                    CFFFL=CFF233
                    CEH=CE253
                    CEL=CE233
                    TH=253.0
                    TL=233.0
                 ELSEIF (T(IPATH) .GT. 213.0 .AND. T(IPATH) .LT.
     +           233.0) THEN
                    CFFFH=CFF233
                    CFFFL=CFF213
                    CEH=CE233
                    CEL=CE213
                    TH=233.0
                    TL=213.0
                 ELSEIF (T(IPATH) .LE. 213.0) THEN
                    CFFFH=CFF213
                    CFFFL=CFF193
                    CEH=CE213
                    CEL=CE193
                    TH=213.0
                    TL=193.0
                 ENDIF
C
C                Linear interpolations in temperature
                 CFFF=CFFFL+(CFFFH-CFFFL)*(T(IPATH)-TL)/(TH-TL)
                 CE=CEL+(CEH-CEL)*(T(IPATH)-TL)/(TH-TL)
C
C                Note: 273.15^2=74610.922
                 CPART=( PARTP(IPATH)*74610.922/(T(IPATH)*T(IPATH)) )*
     +           CFFF*( PARTP(IPATH)+CE*( P(IPATH)-PARTP(IPATH) ) )
C
C                Note: AVOG/(2.6867E+19*273.15) = 82058.508
                 X=82058.508*( AMT(IPATH)*T(IPATH) )/PARTP(IPATH)
C
                 A4 = CPART*X*1.00E-06
ccc
C                note: this only adjusts the entire con, not self/foreign
                 if (xgas .eq. 22) a4=a4*xcon1
ccc
                 WABS(IP,IPW,IPATH)  = WABS(IP,IPW,IPATH) + A4
                 CABS(IP,IPW,IPATH)  = CABS(IP,IPW,IPATH) + A4
               ENDIF
C              *********************************************************
             ENDIF
C
   40      CONTINUE
          ENDIF
   30    CONTINUE 
        ENDIF
   10  CONTINUE
C
       RETURN
       END
