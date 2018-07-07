c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c these water subroutines calculate the continuum for WATER(ID=1)
c and are based entirely in the KCARTA routines, except NO jacobians

c similarly the O2 and N2 subroutines are based ENTIRELY on GENLN2 stuff

c note we are using        REAL CON(*),WHICHLAYER instead of
c                          REAL CON(KMAXLAYER,KMAXPTS)

c also real ------> real*8

c because of memeory problems 
c and so we do        DO 10 ILAY=whichlayer,whichlayer
c instead of          DO 10 ILAY=1,NLAY

c THIS IS JUST FOR THE O2/N2 continuums

c************************************************************************
      include '/home/sergio/SPECTRA/CKDLINUX/dspline.F'
c************************************************************************
c this is the LBLRTMv5.10 code .. set up for arb N2/O2 mix ratio
      include '/home/sergio/SPECTRA/CKDLINUX/oxynit_new_lbl5.1mod.F'
c************************************************************************
c this is the LBLRTMv5.10 code .. set up for N2/O2 = 79/21
c      include 'oxynit_new_lbl5.1orig.F'
c************************************************************************
c this is the original code from GENLN2
c      include 'oxynit_old.F'
c************************************************************************
      
c************************************************************************
c                   WATER CONTINUUM CKD 0,21,23,24
c           this would be the CKD defn (using local lineshape)
c            this would be included from calconwater_loc.F 
c************************************************************************

c************************************************************************
c                   WATER CONTINUUM CKD 0,21,23
c            this is the Lorentz lineshape (from GENLN2)
c            this would be included from calconwater.F 
c************************************************************************
