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

c NOTE THIS IS SAME BASE CODE AS CALCONWATER.F --- THE ONLY DIFFERENCE
c IS THAT THE INCLUDE FILES ARE DIFFERENT!!!!!

c************************************************************************
       include '/home/sergio/SPECTRA/CKDLINUX/dspline.F'
c************************************************************************
c this is MT_CKD2.5
c directly taken from LBLRTM 12.1, and probably from AER MT_CKD2.5 code
      include '/home/sergio/SPECTRA/CKDLINUX/calcon_12_1.F'
